# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Interface of the approximate Jacobians of the implicit time integration methods

"""
    update_jacobian!(jacobian, u_ode, semi, gamma_dt)

Assemble (and factorize) the approximate Jacobian `jacobian` of the right-hand side of the
semidiscretization `semi` at the solution `u_ode`, such that [`solve_jacobian!`](@ref) can
solve linear systems with the matrix `I / gamma_dt - J`. For `jacobian = nothing`, nothing
is done.
"""
update_jacobian!(::Nothing, u_ode, semi, gamma_dt) = nothing

"""
    solve_jacobian!(b_ode, jacobian, semi, gamma_dt)

Solve the linear system `(I / gamma_dt - J) x = b_ode` in place, where `J` is the
approximate Jacobian assembled by [`update_jacobian!`](@ref) and `gamma_dt = gamma * dt`.
For `jacobian = nothing`, i.e., `J = 0`, this scales `b_ode` by `gamma_dt`.
"""
function solve_jacobian!(b_ode, ::Nothing, semi, gamma_dt)
    b_ode .*= gamma_dt
    return nothing
end

###############################################################################
# Infrastructure shared by the vertically implicit Jacobians of the
# `CompressibleEulerPotentialTemperatureEquationsWithGravity3D` on a
# `P4estMeshCubedSphereTopography`: the vertical layout of the degrees of freedom, the
# metric terms of the vertical operator, the copies between Trixi.jl's layout and the
# columns, and the dense LU of the local Schur complements. The
# Jacobians decouple into independent systems for each vertical column of nodes, which are
# solved by KernelAbstractions.jl kernels.

# Launch keywords of the column-parallel GPU kernels (one thread per vertical column) and
# of the level-parallel kernels (one thread per (column, layer)): the workgroup sizes
# measured fastest on an NVIDIA GH200 and on an AMD MI210.
@inline column_kernel_workgroupsize() = (workgroupsize = 64,)
@inline level_kernel_workgroupsize() = (workgroupsize = 128,)

# The number of layers of elements of a `P4estMeshCubedSphereTopography` is stored in the
# registry of the mesh (see `src/meshes/p4est_cubed_sphere_topography.jl`). For parallel
# meshes, the columns must be kept together on one rank.
function number_of_layers(mesh::Trixi.P4estMesh{3})
    key = UInt64(Trixi.P4est.PointerWrappers.pointer(mesh.p4est))
    if !haskey(MeshContextRegistry.REGISTRY, key)
        throw(ArgumentError("the vertical Jacobians require a " *
                            "`P4estMeshCubedSphereTopography`"))
    end
    ctx = MeshContextRegistry.REGISTRY[key]::ColumnPartitionContext
    if mesh isa Trixi.P4estMeshParallel{3} && !ctx.keep_columns_together
        throw(ArgumentError("the vertical Jacobians require " *
                            "`keep_columns_together = true` for parallel meshes"))
    end
    return Int(ctx.num_layers_per_column)
end

function number_of_layers(mesh::Trixi.AbstractMesh)
    throw(ArgumentError("the vertical Jacobians require a " *
                        "`P4estMeshCubedSphereTopography`"))
end

# Column and layer of every node. The trees of the cubed sphere are ordered such that the
# layers of a column are consecutive, with the lowest layer first.
function compute_vertical_mapping(nnodes, nelements, nlayers)
    vertical_line_indices = zeros(Int, nnodes, nnodes, nnodes, nelements)
    vertical_element_indices = zeros(Int, nnodes, nnodes, nnodes, nelements)

    for element in 1:nelements
        layer = (element - 1) % nlayers + 1
        column_offset = ((element - 1) ÷ nlayers) * nnodes^2
        for k in 1:nnodes, j in 1:nnodes, i in 1:nnodes
            vertical_line_indices[i, j, k, element] = column_offset + (j - 1) * nnodes +
                                                      i
            vertical_element_indices[i, j, k, element] = layer
        end
    end

    return vertical_line_indices, vertical_element_indices
end

# Per-element maps used to copy data between Trixi.jl's layout and the columns: the
# horizontal column index and the layer, one per element instead of one per node.
function compute_element_index_maps(vertical_line_indices, vertical_element_indices,
                                    nnodes)
    nelements = size(vertical_line_indices, 4)
    horizontal_index_of_element = [(vertical_line_indices[1, 1, 1, element] - 1) ÷
                                   nnodes^2 + 1 for element in 1:nelements]
    vertical_index_of_element = [vertical_element_indices[1, 1, 1, element]
                                 for element in 1:nelements]
    for element in 1:nelements, k in 1:nnodes, j in 1:nnodes, i in 1:nnodes
        @assert vertical_line_indices[i, j, k, element] ==
                (horizontal_index_of_element[element] - 1) * nnodes^2 +
                (j - 1) * nnodes + i
        @assert vertical_element_indices[i, j, k, element] ==
                vertical_index_of_element[element]
    end
    return horizontal_index_of_element, vertical_index_of_element
end

# Vertical frame of every node: `vertical_frame[1:3, ...]` holds the unit contravariant
# direction `nu = Ja^3 / n` of the third reference coordinate, oriented outwards, and
# `vertical_frame[4, ...]` the norm `n = |Ja^3|` normalized by its global maximum (a
# global constant cancels in the similarity transformation of the copies below)
function compute_vertical_frame(semi)
    (; solver, equations) = semi
    node_coordinates = Adapt.adapt(Array, semi.cache.elements.node_coordinates)
    contravariant_vectors = Adapt.adapt(Array,
                                        semi.cache.elements.contravariant_vectors)
    vertical_frame = zeros(eltype(node_coordinates), 4,
                           size(node_coordinates)[2:end]...)
    area_max = zero(eltype(node_coordinates))
    for element in axes(node_coordinates, 5)
        for k in eachnode(solver), j in eachnode(solver), i in eachnode(solver)
            Ja3 = get_contravariant_vector(3, contravariant_vectors, i, j, k, element)
            area = norm(Ja3)
            x = get_node_coords(node_coordinates, equations, solver, i, j, k, element)
            direction = sign(dot(Ja3, x)) * Ja3 / area
            vertical_frame[1:3, i, j, k, element] .= direction
            vertical_frame[4, i, j, k, element] = area
            area_max = max(area_max, area)
        end
    end
    vertical_frame[4, :, :, :, :] ./= area_max
    return vertical_frame
end

# Metric factor `n / J` of the vertical operator of every node, with the norm
# `n = |Ja^3|` of the contravariant vector of the third reference coordinate and the
# Jacobian `J`, laid out [column, node, layer]
function compute_metric_factor(ncolumns, nnodes, nlayers, semi,
                               vertical_line_indices,
                               vertical_element_indices)
    contravariant_vectors = Adapt.adapt(Array,
                                        semi.cache.elements.contravariant_vectors)
    inverse_jacobian = Adapt.adapt(Array, semi.cache.elements.inverse_jacobian)
    metric_factor = zeros(eltype(inverse_jacobian), ncolumns, nnodes,
                          nlayers)
    for element in axes(inverse_jacobian, 4), k in 1:nnodes, j in 1:nnodes,
        i in 1:nnodes

        Ja3 = get_contravariant_vector(3, contravariant_vectors, i, j, k, element)
        column = vertical_line_indices[i, j, k, element]
        layer = vertical_element_indices[i, j, k, element]
        metric_factor[column, k, layer] = norm(Ja3) *
                                          inverse_jacobian[i, j, k, element]
    end
    return metric_factor
end

# Consistent pointwise gravity g * (nu . r_hat), the projection of the gravitational
# acceleration on the unit contravariant direction `nu` of the third reference coordinate
# (see `compute_vertical_frame`). On a sphere without orography, it reduces to g; with
# terrain-following levels, nu deviates from the radial direction r_hat.
function compute_gravity(ncolumns, nnodes, nlayers, semi, vertical_line_indices,
                         vertical_element_indices)
    (; solver, equations) = semi
    node_coordinates = Adapt.adapt(Array, semi.cache.elements.node_coordinates)
    contravariant_vectors = Adapt.adapt(Array,
                                        semi.cache.elements.contravariant_vectors)
    gravity = zeros(eltype(node_coordinates), ncolumns, nnodes, nlayers)
    for element in axes(node_coordinates, 5)
        for k in eachnode(solver), j in eachnode(solver), i in eachnode(solver)
            Ja3 = get_contravariant_vector(3, contravariant_vectors, i, j, k, element)
            x = get_node_coords(node_coordinates, equations, solver, i, j, k, element)
            # |nu . r_hat| with the outward orientation of `compute_vertical_frame`
            projection = abs(dot(Ja3, x)) / (norm(Ja3) * norm(x))
            column = vertical_line_indices[i, j, k, element]
            layer = vertical_element_indices[i, j, k, element]
            gravity[column, k, layer] = equations.g * projection
        end
    end
    return gravity
end

# Copy data between Trixi.jl's layout `[v, i, j, k, element]` and the layout of the
# columns `[column, v, node, layer]`, with the variables `1 => rho`, `2 => rho_w` (vertical
# momentum), and `3 => rho_theta`

# Copy the vertical variables of the node `(i, j, k, element)` of `u`, multiplied by
# `scaling`, where the vertical momentum is the component of the momentum in the
# `vertical_direction` `nu` (see `compute_vertical_frame`). The `scaling` is one for the
# state and the normalized norm `n` for the right-hand side of the linear system, i.e.,
# the diagonal similarity transformation `N` of the vertical operator.
Base.@propagate_inbounds function copy_node_to_columns!(u_vertical, u, scaling,
                                                        vertical_frame,
                                                        horizontal_index_of_element,
                                                        vertical_index_of_element,
                                                        nnodes,
                                                        i, j, k, element)
    column = (horizontal_index_of_element[element] - 1) * nnodes * nnodes +
             (j - 1) * nnodes + i
    layer = vertical_index_of_element[element]
    vertical_direction = SVector(vertical_frame[1, i, j, k, element],
                                 vertical_frame[2, i, j, k, element],
                                 vertical_frame[3, i, j, k, element])
    rho_v = SVector(u[2, i, j, k, element], u[3, i, j, k, element],
                    u[4, i, j, k, element])
    u_vertical[column, 1, k, layer] = scaling * u[1, i, j, k, element]
    u_vertical[column, 2, k, layer] = scaling * dot(vertical_direction, rho_v)
    u_vertical[column, 3, k, layer] = scaling * u[5, i, j, k, element]
    return nothing
end

# Copy the solution of the linear system, divided by the normalized norm `n`, to `b`. The
# horizontal components of the momentum are multiplied by `gamma_dt`, i.e., solved with
# `J = 0`, and the component in the `vertical_direction` is replaced by the solution.
Base.@propagate_inbounds function copy_node_from_columns!(b, u_vertical, vertical_frame,
                                                          horizontal_index_of_element,
                                                          vertical_index_of_element,
                                                          nnodes, gamma_dt,
                                                          i, j, k, element)
    column = (horizontal_index_of_element[element] - 1) * nnodes * nnodes +
             (j - 1) * nnodes + i
    layer = vertical_index_of_element[element]
    vertical_direction = SVector(vertical_frame[1, i, j, k, element],
                                 vertical_frame[2, i, j, k, element],
                                 vertical_frame[3, i, j, k, element])
    inv_area = inv(vertical_frame[4, i, j, k, element])
    rho_v = SVector(b[2, i, j, k, element], b[3, i, j, k, element],
                    b[4, i, j, k, element])
    rho_w = inv_area * u_vertical[column, 2, k, layer]
    rho_v = gamma_dt * rho_v +
            (rho_w - gamma_dt * dot(vertical_direction, rho_v)) * vertical_direction
    b[1, i, j, k, element] = inv_area * u_vertical[column, 1, k, layer]
    b[2, i, j, k, element] = rho_v[1]
    b[3, i, j, k, element] = rho_v[2]
    b[4, i, j, k, element] = rho_v[3]
    b[5, i, j, k, element] = inv_area * u_vertical[column, 3, k, layer]
    return nothing
end

@kernel inbounds=true function copy_state_to_columns_KAkernel!(u_vertical, @Const(u),
                                                               @Const(vertical_frame),
                                                               @Const(horizontal_index_of_element),
                                                               @Const(vertical_index_of_element),
                                                               nnodes)
    i, j, k, element = @index(Global, NTuple)
    copy_node_to_columns!(u_vertical, u, one(eltype(u_vertical)), vertical_frame,
                          horizontal_index_of_element, vertical_index_of_element,
                          nnodes, i, j, k, element)
end

@kernel inbounds=true function copy_rhs_to_columns_KAkernel!(u_vertical, @Const(b),
                                                             @Const(vertical_frame),
                                                             @Const(horizontal_index_of_element),
                                                             @Const(vertical_index_of_element),
                                                             nnodes)
    i, j, k, element = @index(Global, NTuple)
    copy_node_to_columns!(u_vertical, b, vertical_frame[4, i, j, k, element],
                          vertical_frame, horizontal_index_of_element,
                          vertical_index_of_element, nnodes, i, j, k, element)
end

@kernel inbounds=true function copy_solution_from_columns_KAkernel!(b,
                                                                    @Const(u_vertical),
                                                                    @Const(vertical_frame),
                                                                    @Const(horizontal_index_of_element),
                                                                    @Const(vertical_index_of_element),
                                                                    nnodes, gamma_dt)
    i, j, k, element = @index(Global, NTuple)
    copy_node_from_columns!(b, u_vertical, vertical_frame, horizontal_index_of_element,
                            vertical_index_of_element, nnodes, gamma_dt, i, j, k,
                            element)
end

# Copy the state `u`, from which the Jacobian is assembled, to the columns
function copy_state_to_columns!(backend::Backend, u_vertical, u, jacobian, semi)
    (; solver) = semi
    (; vertical_frame, horizontal_index_of_element, vertical_index_of_element) = jacobian
    nnodes_ = number_of_nodes(jacobian)
    kernel! = copy_state_to_columns_KAkernel!(backend)
    kernel!(u_vertical, u, vertical_frame, horizontal_index_of_element,
            vertical_index_of_element, nnodes_;
            ndrange = (nnodes_, nnodes_, nnodes_, nelements(solver, semi.cache)))
    return nothing
end

# Copy the right-hand side `b` of the linear system, multiplied by the normalized norm
# `n`, to the columns
function copy_rhs_to_columns!(backend::Backend, u_vertical, b, jacobian, semi)
    (; solver) = semi
    (; vertical_frame, horizontal_index_of_element, vertical_index_of_element) = jacobian
    nnodes_ = number_of_nodes(jacobian)
    kernel! = copy_rhs_to_columns_KAkernel!(backend)
    kernel!(u_vertical, b, vertical_frame, horizontal_index_of_element,
            vertical_index_of_element, nnodes_;
            ndrange = (nnodes_, nnodes_, nnodes_, nelements(solver, semi.cache)))
    return nothing
end

# Copy the solution of the linear system, divided by the normalized norm `n`, from the
# columns to `b`
function copy_solution_from_columns!(backend::Backend, b, u_vertical, jacobian, semi,
                                     gamma_dt)
    (; solver) = semi
    (; vertical_frame, horizontal_index_of_element, vertical_index_of_element) = jacobian
    nnodes_ = number_of_nodes(jacobian)
    kernel! = copy_solution_from_columns_KAkernel!(backend)
    kernel!(b, u_vertical, vertical_frame, horizontal_index_of_element,
            vertical_index_of_element, nnodes_, gamma_dt;
            ndrange = (nnodes_, nnodes_, nnodes_, nelements(solver, semi.cache)))
    return nothing
end

###############################################################################
# Dense LU factorization and solve of the local Schur complement `SA` of one element
# (indexed `[column, i, j, layer]`, `N = nnodes - 2`)
Base.@propagate_inbounds function lu_dense!(layer, column, A, ::Val{N}) where {N}
    for j in 1:(N - 1)
        # One division and (N - j) multiplications instead of (N - j) divisions
        inv_jj = inv(A[column, j, j, layer])
        for i in (j + 1):N
            A[column, i, j, layer] *= inv_jj
        end
        for k in (j + 1):N
            a_jk = A[column, j, k, layer]
            @simd for i in (j + 1):N
                A[column, i, k, layer] -= A[column, i, j, layer] * a_jk
            end
        end
    end
end

Base.@propagate_inbounds function ldiv_dense!(layer, column, A, b, ::Val{N}) where {N}
    # Forward loop
    for k in 1:(N - 1)
        for i in (k + 1):N
            b[i] -= A[column, i, k, layer] * b[k]
        end
    end
    # Backward loop
    for k in N:-1:1
        b[k] /= A[column, k, k, layer]
        for i in 1:(k - 1)
            b[i] -= A[column, i, k, layer] * b[k]
        end
    end
end
end # @muladd

include("jacobian_two_sided_gpu.jl")
