@muladd begin
#! format: noindent

function Trixi.rhs!(du, u, t,
                    mesh::Union{P4estMesh{2}, T8codeMesh{2}},
                    equations::AbstractEquations{3},
                    boundary_conditions, source_terms::Source,
                    dg::DG, cache) where {Source}

    # Reset du
    Trixi.@trixi_timeit Trixi.timer() "reset ∂u/∂t" Trixi.reset_du!(du, dg, cache)

    # Calculate volume integral
    Trixi.@trixi_timeit Trixi.timer() "volume integral" begin
        Trixi.calc_volume_integral!(du, u, mesh,
                                    Trixi.have_nonconservative_terms(equations),
                                    equations,
                                    dg.volume_integral, dg, cache)
    end

    # Prolong solution to interfaces
    Trixi.@trixi_timeit Trixi.timer() "prolong2interfaces" begin
        Trixi.prolong2interfaces!(cache, u, mesh, equations, dg)
    end

    # Calculate interface fluxes
    Trixi.@trixi_timeit Trixi.timer() "interface flux" begin
        Trixi.calc_interface_flux!(cache.elements.surface_flux_values, mesh,
                                   Trixi.have_nonconservative_terms(equations),
                                   equations,
                                   dg.surface_integral, dg, cache)
    end

    # Prolong solution to boundaries
    Trixi.@trixi_timeit Trixi.timer() "prolong2boundaries" begin
        Trixi.prolong2boundaries!(cache, u, mesh, equations,
                                  dg.surface_integral, dg)
    end

    # Calculate boundary fluxes
    Trixi.@trixi_timeit Trixi.timer() "boundary flux" begin
        Trixi.calc_boundary_flux!(cache, t, boundary_conditions, mesh, equations,
                                  dg.surface_integral, dg)
    end

    # Prolong solution to mortars
    Trixi.@trixi_timeit Trixi.timer() "prolong2mortars" begin
        Trixi.prolong2mortars!(cache, u, mesh, equations, dg.mortar, dg)
    end

    # Calculate mortar fluxes
    Trixi.@trixi_timeit Trixi.timer() "mortar flux" begin
        Trixi.calc_mortar_flux!(cache.elements.surface_flux_values, mesh,
                                Trixi.have_nonconservative_terms(equations), equations,
                                dg.mortar, dg.surface_integral, dg, cache)
    end

    # Calculate surface integrals
    Trixi.@trixi_timeit Trixi.timer() "surface integral" begin
        Trixi.calc_surface_integral!(du, u, mesh, equations,
                                     dg.surface_integral, dg, cache)
    end

    # Apply Jacobian from mapping to reference element
    Trixi.@trixi_timeit Trixi.timer() "Jacobian" Trixi.apply_jacobian!(du, mesh,
                                                                       equations,
                                                                       dg, cache)

    # Calculate source terms
    Trixi.@trixi_timeit Trixi.timer() "source terms" begin
        calc_sources_2d_manifold_in_3d!(du, u, t, source_terms, equations, dg, cache)
    end

    return nothing
end

# Weak-form kernel for 3D equations solved in 2D manifolds
@inline function Trixi.weak_form_kernel!(du, u,
                                         element,
                                         mesh::Union{StructuredMesh{2},
                                                     UnstructuredMesh2D,
                                                     P4estMesh{2}, T8codeMesh{2}},
                                         nonconservative_terms::False,
                                         equations::AbstractEquations{3},
                                         dg::DGSEM, cache, alpha = true)
    # true * [some floating point value] == [exactly the same floating point value]
    # This can (hopefully) be optimized away due to constant propagation.
    @unpack derivative_dhat = dg.basis
    @unpack contravariant_vectors = cache.elements

    for j in eachnode(dg), i in eachnode(dg)
        u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)

        flux1 = flux(u_node, 1, equations)
        flux2 = flux(u_node, 2, equations)
        flux3 = flux(u_node, 3, equations)

        # Compute the contravariant flux by taking the scalar product of the
        # first contravariant vector Ja^1 and the flux vector
        Ja11, Ja12,
        Ja13 = Trixi.get_contravariant_vector(1, contravariant_vectors, i,
                                              j,
                                              element)
        contravariant_flux1 = Ja11 * flux1 + Ja12 * flux2 + Ja13 * flux3
        for ii in eachnode(dg)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_dhat[ii, i],
                                             contravariant_flux1, equations, dg, ii, j,
                                             element)
        end

        # Compute the contravariant flux by taking the scalar product of the
        # second contravariant vector Ja^2 and the flux vector
        Ja21, Ja22,
        Ja23 = Trixi.get_contravariant_vector(2, contravariant_vectors, i,
                                              j,
                                              element)
        contravariant_flux2 = Ja21 * flux1 + Ja22 * flux2 + Ja23 * flux3
        for jj in eachnode(dg)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_dhat[jj, j],
                                             contravariant_flux2, equations, dg, i, jj,
                                             element)
        end
    end

    return nothing
end

function calc_sources_2d_manifold_in_3d!(du, u, t, source_terms::Nothing,
                                         equations::AbstractEquations{3}, dg::DG, cache)
    return nothing
end

function calc_sources_2d_manifold_in_3d!(du, u, t, source_terms,
                                         equations::AbstractEquations{3}, dg::DG, cache)
    @unpack node_coordinates, contravariant_vectors, inverse_jacobian = cache.elements

    Trixi.@threaded for element in eachelement(dg, cache)
        for j in eachnode(dg), i in eachnode(dg)
            u_local = Trixi.get_node_vars(u, equations, dg, i, j, element)
            du_local = Trixi.get_node_vars(du, equations, dg, i, j, element)
            x_local = Trixi.get_node_coords(node_coordinates, equations, dg,
                                            i, j, element)
            contravariant_normal_vector = Trixi.get_contravariant_vector(3,
                                                                         contravariant_vectors,
                                                                         i, j,
                                                                         element) *
                                          inverse_jacobian[i, j, element]
            source = source_terms(u_local, du_local, x_local, t, equations,
                                  contravariant_normal_vector)
            Trixi.add_to_node_vars!(du, source, equations, dg, i, j, element)
        end
    end

    return nothing
end
end # @muladd
