# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Trixi.jl neither creates the auxiliary variables of a `SemidiscretizationHyperbolicParabolic`
# nor passes them to the parabolic fluxes and boundary conditions. For the
# `TurbulentDiffusion3D`, the constructor below adds the auxiliary variables of the hyperbolic
# equations, the upward unit vectors `k` of all nodes, the number of layers of the columns of
# the `P4estMeshCubedSphereTopography`, and the interpolation to the lowest levels of the
# columns to the `cache`. Then, `calc_parabolic_fluxes!` passes `k` and the variables of the
# lowest level of the column to `flux` and stores the lowest levels in `cache.lowest_level`,
# which the `BoundaryConditionSurfaceFluxes` uses.
function Trixi.SemidiscretizationHyperbolicParabolic(mesh::P4estMesh{3},
                                                     equations::Tuple{Any,
                                                                      TurbulentDiffusion3D},
                                                     initial_condition, solver;
                                                     aux_field = nothing, kwargs...)
    semi = invoke(Trixi.SemidiscretizationHyperbolicParabolic,
                  Tuple{Any, Tuple, Any, Any}, mesh, equations, initial_condition,
                  solver;
                  kwargs...)
    equations_hyperbolic = first(equations)
    if aux_field === nothing
        throw(ArgumentError("`TurbulentDiffusion3D` requires an `aux_field` whose first " *
                            "auxiliary variable is the geopotential"))
    end

    nlayers_ = nlayers(mesh)
    check_vertical_ordering(nlayers_, semi)

    (; cache) = semi
    aux_vars = Trixi.init_aux_vars(mesh, equations_hyperbolic, solver, cache, aux_field)
    upward_unit_vectors = compute_upward_unit_vectors(aux_field, equations_hyperbolic,
                                                      solver, cache)
    lowest_level_weights, lowest_level_directions = compute_lowest_levels(last(equations).lowest_level_height,
                                                                          aux_field,
                                                                          nlayers_,
                                                                          equations_hyperbolic,
                                                                          solver, cache)
    # Specific variables, upward unit vector, and height above the surface of the lowest
    # levels, laid out as `[variable, i, j, column]`
    nnodes_ = nnodes(solver)
    lowest_level = zeros(eltype(cache.elements.node_coordinates),
                         nvariables(equations_hyperbolic) + 4, nnodes_, nnodes_,
                         nelements(solver, cache) ÷ nlayers_)
    return setproperties(semi;
                         cache = (; cache..., aux_vars, upward_unit_vectors,
                                  lowest_level_weights, lowest_level_directions,
                                  lowest_level, nlayers = nlayers_))
end

# Lowest levels of the columns at the height `height` above the surface: the values
# `l_k(xi_a)` of the Lagrange polynomials of the nodes `k` of the lowest element at the
# reference coordinate `xi_a` of the lowest level, laid out as `[k, i, j, column]`, and the
# upward unit vectors there, laid out as `[component, i, j, column]`. The height of a point of
# the column above the surface is the difference of the distances to the center, and `xi_a`
# is the root of the polynomial of the heights of the nodes minus `height`.
function compute_lowest_levels(height, aux_field, nlayers_, equations, solver, cache)
    (; node_coordinates) = cache.elements
    (; nodes) = solver.basis
    nnodes_ = nnodes(solver)
    ncolumns = nelements(solver, cache) ÷ nlayers_
    RealT = eltype(node_coordinates)
    lowest_level_weights = zeros(RealT, nnodes_, nnodes_, nnodes_, ncolumns)
    lowest_level_directions = zeros(RealT, 3, nnodes_, nnodes_, ncolumns)
    for column in 1:ncolumns, j in Trixi.eachnode(solver), i in Trixi.eachnode(solver)
        element = (column - 1) * nlayers_ + 1
        x_nodes = [Trixi.get_node_coords(node_coordinates, equations, solver, i, j, k,
                                         element) for k in Trixi.eachnode(solver)]
        heights = norm.(x_nodes) .- norm(first(x_nodes))
        if !(0 < height < last(heights))
            throw(ArgumentError("the height $height of the lowest level must lie in the " *
                                "lowest element of the columns, whose height is " *
                                "$(last(heights))"))
        end
        height_polynomial(xi) = sum(lagrange_polynomials(xi, nodes) .* heights)
        xi = -1 + 2 * height / last(heights)
        for _ in 1:100
            delta = (height_polynomial(xi) - height) /
                    ForwardDiff.derivative(height_polynomial, xi)
            xi -= delta
            abs(delta) <= 10 * eps(RealT) && break
        end
        weights = lagrange_polynomials(xi, nodes)
        lowest_level_weights[:, i, j, column] .= weights
        x_a = sum(weights .* x_nodes)
        lowest_level_directions[:, i, j, column] .= upward_unit_vector(x_a, aux_field,
                                                                       equations)
    end
    return lowest_level_weights, lowest_level_directions
end

# Values of the Lagrange polynomials of the `nodes` at `xi`
function lagrange_polynomials(xi, nodes)
    return SVector(ntuple(length(nodes)) do k
                       prod((xi - nodes[m]) / (nodes[k] - nodes[m])
                            for m in eachindex(nodes) if m != k)
                   end)
end

# Upward unit vectors `k` of all nodes, laid out as `[component, i, j, k, element]`
function compute_upward_unit_vectors(aux_field, equations, solver, cache)
    (; node_coordinates) = cache.elements
    upward_unit_vectors = similar(node_coordinates)
    Trixi.@threaded for element in Trixi.eachelement(solver, cache)
        for k in Trixi.eachnode(solver), j in Trixi.eachnode(solver),
            i in Trixi.eachnode(solver)

            x = Trixi.get_node_coords(node_coordinates, equations, solver, i, j, k,
                                      element)
            k_hat = upward_unit_vector(x, aux_field, equations)
            for d in 1:3
                upward_unit_vectors[d, i, j, k, element] = k_hat[d]
            end
        end
    end
    return upward_unit_vectors
end

# Arrays of the columns passed to the kernels
function columns_parabolic(cache)
    (; upward_unit_vectors, lowest_level_weights, lowest_level_directions, lowest_level,
    nlayers) = cache
    return (; upward_unit_vectors, lowest_level_weights, lowest_level_directions,
            lowest_level, nlayers)
end

function Trixi.calc_parabolic_fluxes!(flux_parabolic, gradients, u_transformed,
                                      mesh::P4estMesh{3},
                                      equations_parabolic::TurbulentDiffusion3D,
                                      dg::DG, cache)
    columns = columns_parabolic(cache)
    # Explicit bounds check, which allows us to assume inbounds access below
    @boundscheck begin
        check_axes_parabolic_fluxes(flux_parabolic, gradients, u_transformed, columns,
                                    mesh, equations_parabolic, dg, cache)
    end

    Trixi.@threaded for element in Trixi.eachelement(dg, cache)
        for k in Trixi.eachnode(dg), j in Trixi.eachnode(dg), i in Trixi.eachnode(dg)
            @inbounds calc_parabolic_fluxes_per_node!(flux_parabolic, gradients,
                                                      u_transformed, columns,
                                                      equations_parabolic, dg,
                                                      i, j, k, element)
        end
    end

    return nothing
end

function Trixi.calc_parabolic_fluxes!(backend::Backend, flux_parabolic, gradients,
                                      u_transformed, mesh::P4estMesh{3},
                                      equations_parabolic::TurbulentDiffusion3D,
                                      dg::DG, cache)
    Trixi.nelements(dg, cache) == 0 && return nothing
    columns = columns_parabolic(cache)
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck begin
        check_axes_parabolic_fluxes(flux_parabolic, gradients, u_transformed, columns,
                                    mesh, equations_parabolic, dg, cache)
    end

    kernel! = calc_parabolic_fluxes_KAkernel!(backend)
    kernel!(flux_parabolic, gradients, u_transformed, columns, equations_parabolic, dg,
            ndrange = (Trixi.element_node_ndrange(mesh, dg)...,
                       Trixi.nelements(dg, cache)))
    return nothing
end

@kernel inbounds=true function calc_parabolic_fluxes_KAkernel!(flux_parabolic,
                                                               gradients,
                                                               u_transformed,
                                                               columns,
                                                               equations_parabolic,
                                                               dg)
    i, j, k, element = @index(Global, NTuple)
    calc_parabolic_fluxes_per_node!(flux_parabolic, gradients, u_transformed, columns,
                                    equations_parabolic, dg, i, j, k, element)
end

Base.@propagate_inbounds function calc_parabolic_fluxes_per_node!(flux_parabolic,
                                                                  gradients,
                                                                  u_transformed,
                                                                  columns,
                                                                  equations_parabolic,
                                                                  dg, i, j, k, element)
    flux_parabolic_x, flux_parabolic_y, flux_parabolic_z = flux_parabolic
    gradients_x, gradients_y, gradients_z = gradients
    (; upward_unit_vectors, lowest_level_weights, lowest_level_directions, lowest_level,
    nlayers) = columns

    w_node = Trixi.get_node_vars(u_transformed, equations_parabolic, dg,
                                 i, j, k, element)
    gradients_node = (Trixi.get_node_vars(gradients_x, equations_parabolic, dg,
                                          i, j, k, element),
                      Trixi.get_node_vars(gradients_y, equations_parabolic, dg,
                                          i, j, k, element),
                      Trixi.get_node_vars(gradients_z, equations_parabolic, dg,
                                          i, j, k, element))
    k_node = get_upward_unit_vector(upward_unit_vectors, i, j, k, element)

    # Lowest level of the column in the lowest element, whose layers are consecutive
    # elements, see `compute_element_index_maps`, with the conservative variables
    # interpolated by the polynomial of the element
    column = (element - 1) ÷ nlayers + 1
    element_lowest = (column - 1) * nlayers + 1
    u_a = zero(w_node)
    for kk in Trixi.eachnode(dg)
        w_kk = Trixi.get_node_vars(u_transformed, equations_parabolic, dg,
                                   i, j, kk, element_lowest)
        u_a = u_a + lowest_level_weights[kk, i, j, column] * specific2cons(w_kk)
    end
    w_a = cons2specific(u_a, equations_parabolic)
    k_a = SVector(lowest_level_directions[1, i, j, column],
                  lowest_level_directions[2, i, j, column],
                  lowest_level_directions[3, i, j, column])
    v_a = SVector(w_a[2], w_a[3], w_a[4])
    wind_speed = norm(v_a - dot(v_a, k_a) * k_a)
    z_a = equations_parabolic.lowest_level_height
    column_variables = SVector(wind_speed, z_a)

    if k == 1 && element == element_lowest
        for v in eachindex(w_a)
            lowest_level[v, i, j, column] = w_a[v]
        end
        nvars = length(w_a)
        for d in 1:3
            lowest_level[nvars + d, i, j, column] = k_a[d]
        end
        lowest_level[nvars + 4, i, j, column] = z_a
    end

    Trixi.set_node_vars!(flux_parabolic_x,
                         flux(w_node, gradients_node, k_node, column_variables, 1,
                              equations_parabolic),
                         equations_parabolic, dg, i, j, k, element)
    Trixi.set_node_vars!(flux_parabolic_y,
                         flux(w_node, gradients_node, k_node, column_variables, 2,
                              equations_parabolic),
                         equations_parabolic, dg, i, j, k, element)
    Trixi.set_node_vars!(flux_parabolic_z,
                         flux(w_node, gradients_node, k_node, column_variables, 3,
                              equations_parabolic),
                         equations_parabolic, dg, i, j, k, element)
    return nothing
end

Base.@propagate_inbounds function get_upward_unit_vector(upward_unit_vectors, i, j, k,
                                                         element)
    return SVector(upward_unit_vectors[1, i, j, k, element],
                   upward_unit_vectors[2, i, j, k, element],
                   upward_unit_vectors[3, i, j, k, element])
end

function check_axes_parabolic_fluxes(flux_parabolic, gradients, u_transformed, columns,
                                     mesh, equations_parabolic, dg, cache)
    Trixi.check_axes(u_transformed, mesh, equations_parabolic, dg, cache)
    for gradient in gradients
        Trixi.check_axes(gradient, mesh, equations_parabolic, dg, cache)
    end
    for flux_parabolic_i in flux_parabolic
        Trixi.check_axes(flux_parabolic_i, mesh, equations_parabolic, dg, cache)
    end
    (; upward_unit_vectors, lowest_level_weights, lowest_level_directions, lowest_level,
    nlayers) = columns
    node_axes = (Trixi.eachnode(dg), Trixi.eachnode(dg), Trixi.eachnode(dg),
                 Trixi.eachelement(dg, cache))
    column_axes = (Trixi.eachnode(dg), Trixi.eachnode(dg),
                   Base.OneTo(Trixi.nelements(dg, cache) ÷ nlayers))
    Trixi.check_axes(upward_unit_vectors, (Base.OneTo(3), node_axes...))
    Trixi.check_axes(lowest_level_weights, (Trixi.eachnode(dg), column_axes...))
    Trixi.check_axes(lowest_level_directions, (Base.OneTo(3), column_axes...))
    Trixi.check_axes(lowest_level,
                     (Base.OneTo(nvariables(equations_parabolic) + 4),
                      Trixi.eachnode(dg),
                      Trixi.eachnode(dg),
                      Base.OneTo(Trixi.nelements(dg, cache) ÷ nlayers)))
    return nothing
end

# Fluxes of the `BoundaryConditionSurfaceFluxes` through the boundary faces of the
# `boundary_condition_indices` from the lowest levels of the columns
function Trixi.calc_boundary_flux!(cache, t,
                                   boundary_condition::BoundaryConditionSurfaceFluxes,
                                   boundary_condition_indices,
                                   operator_type::Trixi.Divergence, mesh::P4estMesh{3},
                                   equations_parabolic::TurbulentDiffusion3D,
                                   surface_integral, dg::DG)
    (; boundaries) = cache
    (; node_coordinates, contravariant_vectors, surface_flux_values) = cache.elements
    (; lowest_level, nlayers) = cache
    index_range = Trixi.eachnode(dg)

    Trixi.@threaded for local_index in eachindex(boundary_condition_indices)
        boundary = boundary_condition_indices[local_index]
        for j in Trixi.eachnode(dg), i in Trixi.eachnode(dg)
            @inbounds calc_surface_flux_per_node!(surface_flux_values, lowest_level,
                                                  nlayers, boundary_condition,
                                                  equations_parabolic, dg, boundary,
                                                  boundaries.neighbor_ids,
                                                  boundaries.node_indices, index_range,
                                                  node_coordinates,
                                                  contravariant_vectors, i, j)
        end
    end
    return nothing
end

function Trixi.calc_boundary_flux_by_type!(backend::Backend, cache, t, flux_parabolic,
                                           BCs::Tuple{BoundaryConditionSurfaceFluxes,
                                                      Vararg{Any}},
                                           BC_indices::Tuple{AbstractVector{Int},
                                                             Vararg{AbstractVector{Int}}},
                                           operator_type::Trixi.Divergence,
                                           mesh::P4estMesh{3},
                                           equations_parabolic::TurbulentDiffusion3D,
                                           dg::DG)
    boundary_condition = first(BCs)
    boundary_condition_indices = first(BC_indices)
    n_boundaries = length(boundary_condition_indices)
    if n_boundaries > 0
        (; neighbor_ids, node_indices) = cache.boundaries
        (; node_coordinates, contravariant_vectors, surface_flux_values) = cache.elements
        (; lowest_level, nlayers) = cache
        kernel! = calc_surface_flux_KAkernel!(backend)
        kernel!(surface_flux_values, lowest_level, nlayers, boundary_condition,
                equations_parabolic, dg, boundary_condition_indices, neighbor_ids,
                node_indices, Trixi.eachnode(dg), node_coordinates,
                contravariant_vectors,
                ndrange = (Trixi.boundary_node_ndrange(mesh, dg)..., n_boundaries))
    end

    Trixi.calc_boundary_flux_by_type!(backend, cache, t, flux_parabolic, Base.tail(BCs),
                                      Base.tail(BC_indices), operator_type, mesh,
                                      equations_parabolic, dg)
    return nothing
end

@kernel inbounds=true function calc_surface_flux_KAkernel!(surface_flux_values,
                                                           lowest_level, nlayers,
                                                           boundary_condition,
                                                           equations_parabolic, dg,
                                                           boundary_condition_indices,
                                                           neighbor_ids, node_indices,
                                                           index_range,
                                                           node_coordinates,
                                                           contravariant_vectors)
    i, j, local_index = @index(Global, NTuple)
    boundary = boundary_condition_indices[local_index]
    calc_surface_flux_per_node!(surface_flux_values, lowest_level, nlayers,
                                boundary_condition, equations_parabolic, dg, boundary,
                                neighbor_ids, node_indices, index_range,
                                node_coordinates,
                                contravariant_vectors, i, j)
end

Base.@propagate_inbounds function calc_surface_flux_per_node!(surface_flux_values,
                                                              lowest_level, nlayers,
                                                              boundary_condition,
                                                              equations_parabolic, dg,
                                                              boundary, neighbor_ids,
                                                              node_indices_arr,
                                                              index_range,
                                                              node_coordinates,
                                                              contravariant_vectors, i,
                                                              j)
    element = neighbor_ids[boundary]
    node_indices = node_indices_arr[boundary]
    direction = Trixi.indices2direction(node_indices)
    n = length(index_range)

    i_node_start, i_node_step_i, i_node_step_j = Trixi.index_to_start_step_3d(node_indices[1],
                                                                              index_range)
    j_node_start, j_node_step_i, j_node_step_j = Trixi.index_to_start_step_3d(node_indices[2],
                                                                              index_range)
    k_node_start, k_node_step_i, k_node_step_j = Trixi.index_to_start_step_3d(node_indices[3],
                                                                              index_range)
    i_node = Trixi.delayed_index_3d(i_node_start, i_node_step_i, i_node_step_j, i, j, n)
    j_node = Trixi.delayed_index_3d(j_node_start, j_node_step_i, j_node_step_j, i, j, n)
    k_node = Trixi.delayed_index_3d(k_node_start, k_node_step_i, k_node_step_j, i, j, n)

    normal_direction = Trixi.get_normal_direction(direction, contravariant_vectors,
                                                  i_node, j_node, k_node, element)
    x = SVector(node_coordinates[1, i_node, j_node, k_node, element],
                node_coordinates[2, i_node, j_node, k_node, element],
                node_coordinates[3, i_node, j_node, k_node, element])

    # The surface is the lower face of the lowest element of the column
    column = (element - 1) ÷ nlayers + 1
    nvars = nvariables(equations_parabolic)
    w_a = SVector(ntuple(v -> lowest_level[v, i_node, j_node, column], Val(nvars)))
    k_a = SVector(lowest_level[nvars + 1, i_node, j_node, column],
                  lowest_level[nvars + 2, i_node, j_node, column],
                  lowest_level[nvars + 3, i_node, j_node, column])
    z_a = lowest_level[nvars + 4, i_node, j_node, column]

    flux_ = boundary_condition(w_a, k_a, z_a, normal_direction, x, equations_parabolic)
    for v in Trixi.eachvariable(equations_parabolic)
        surface_flux_values[v, i, j, direction, element] = flux_[v]
    end
    return nothing
end
end # @muladd
