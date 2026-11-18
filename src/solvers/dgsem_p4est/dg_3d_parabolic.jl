# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Trixi.jl neither creates the auxiliary variables of a `SemidiscretizationHyperbolicParabolic`
# nor passes them to the parabolic fluxes. For the `TurbulentDiffusion3D`, the constructor
# below adds the auxiliary variables of the hyperbolic equations and the upward unit vectors
# `k` of all nodes to the `cache`, and `calc_parabolic_fluxes!` passes `k` to `flux`.
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

    (; cache) = semi
    aux_vars = Trixi.init_aux_vars(mesh, equations_hyperbolic, solver, cache, aux_field)
    upward_unit_vectors = compute_upward_unit_vectors(aux_field, equations_hyperbolic,
                                                      solver, cache)
    return setproperties(semi; cache = (; cache..., aux_vars, upward_unit_vectors))
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

function Trixi.calc_parabolic_fluxes!(flux_parabolic, gradients, u_transformed,
                                      mesh::P4estMesh{3},
                                      equations_parabolic::TurbulentDiffusion3D,
                                      dg::DG, cache)
    (; upward_unit_vectors) = cache
    # Explicit bounds check, which allows us to assume inbounds access below
    @boundscheck begin
        check_axes_parabolic_fluxes(flux_parabolic, gradients, u_transformed,
                                    upward_unit_vectors, mesh, equations_parabolic, dg,
                                    cache)
    end

    Trixi.@threaded for element in Trixi.eachelement(dg, cache)
        for k in Trixi.eachnode(dg), j in Trixi.eachnode(dg), i in Trixi.eachnode(dg)
            @inbounds calc_parabolic_fluxes_per_node!(flux_parabolic, gradients,
                                                      u_transformed,
                                                      upward_unit_vectors,
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
    (; upward_unit_vectors) = cache
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck begin
        check_axes_parabolic_fluxes(flux_parabolic, gradients, u_transformed,
                                    upward_unit_vectors, mesh, equations_parabolic, dg,
                                    cache)
    end

    kernel! = calc_parabolic_fluxes_KAkernel!(backend)
    kernel!(flux_parabolic, gradients, u_transformed, upward_unit_vectors,
            equations_parabolic, dg,
            ndrange = (Trixi.element_node_ndrange(mesh, dg)...,
                       Trixi.nelements(dg, cache)))
    return nothing
end

@kernel inbounds=true function calc_parabolic_fluxes_KAkernel!(flux_parabolic,
                                                               gradients,
                                                               u_transformed,
                                                               upward_unit_vectors,
                                                               equations_parabolic,
                                                               dg)
    i, j, k, element = @index(Global, NTuple)
    calc_parabolic_fluxes_per_node!(flux_parabolic, gradients, u_transformed,
                                    upward_unit_vectors, equations_parabolic, dg,
                                    i, j, k, element)
end

Base.@propagate_inbounds function calc_parabolic_fluxes_per_node!(flux_parabolic,
                                                                  gradients,
                                                                  u_transformed,
                                                                  upward_unit_vectors,
                                                                  equations_parabolic,
                                                                  dg, i, j, k, element)
    flux_parabolic_x, flux_parabolic_y, flux_parabolic_z = flux_parabolic
    gradients_x, gradients_y, gradients_z = gradients

    w_node = Trixi.get_node_vars(u_transformed, equations_parabolic, dg,
                                 i, j, k, element)
    gradients_node = (Trixi.get_node_vars(gradients_x, equations_parabolic, dg,
                                          i, j, k, element),
                      Trixi.get_node_vars(gradients_y, equations_parabolic, dg,
                                          i, j, k, element),
                      Trixi.get_node_vars(gradients_z, equations_parabolic, dg,
                                          i, j, k, element))
    k_node = SVector(upward_unit_vectors[1, i, j, k, element],
                     upward_unit_vectors[2, i, j, k, element],
                     upward_unit_vectors[3, i, j, k, element])

    Trixi.set_node_vars!(flux_parabolic_x,
                         flux(w_node, gradients_node, k_node, 1, equations_parabolic),
                         equations_parabolic, dg, i, j, k, element)
    Trixi.set_node_vars!(flux_parabolic_y,
                         flux(w_node, gradients_node, k_node, 2, equations_parabolic),
                         equations_parabolic, dg, i, j, k, element)
    Trixi.set_node_vars!(flux_parabolic_z,
                         flux(w_node, gradients_node, k_node, 3, equations_parabolic),
                         equations_parabolic, dg, i, j, k, element)
    return nothing
end

function check_axes_parabolic_fluxes(flux_parabolic, gradients, u_transformed,
                                     upward_unit_vectors, mesh, equations_parabolic, dg,
                                     cache)
    Trixi.check_axes(u_transformed, mesh, equations_parabolic, dg, cache)
    for gradient in gradients
        Trixi.check_axes(gradient, mesh, equations_parabolic, dg, cache)
    end
    for flux_parabolic_i in flux_parabolic
        Trixi.check_axes(flux_parabolic_i, mesh, equations_parabolic, dg, cache)
    end
    Trixi.check_axes(upward_unit_vectors,
                     (Base.OneTo(3), Trixi.eachnode(dg), Trixi.eachnode(dg),
                      Trixi.eachnode(dg), Trixi.eachelement(dg, cache)))
    return nothing
end
end # @muladd
