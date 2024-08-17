# Weak-form kernel for 3D equations solved in 2D manifolds
@inline function Trixi.weak_form_kernel!(du, u,
                                         element,
                                         mesh::Union{StructuredMesh{2}, UnstructuredMesh2D,
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
        Ja11, Ja12, Ja13 = Trixi.get_contravariant_vector(1, contravariant_vectors, i, j,
                                                          element)
        contravariant_flux1 = Ja11 * flux1 + Ja12 * flux2 + Ja13 * flux3
        for ii in eachnode(dg)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_dhat[ii, i],
                                             contravariant_flux1, equations, dg, ii, j,
                                             element)
        end

        # Compute the contravariant flux by taking the scalar product of the
        # second contravariant vector Ja^2 and the flux vector
        Ja21, Ja22, Ja23 = Trixi.get_contravariant_vector(2, contravariant_vectors, i, j,
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
