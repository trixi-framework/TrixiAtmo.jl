# Specialization of max_dt function for 3D equations in 2D manifolds
function Trixi.max_dt(u, t,
                mesh::Union{StructuredMesh{2}, UnstructuredMesh2D, P4estMesh{2},
                            T8codeMesh{2}},
                constant_speed::False, equations::AbstractEquations{3}, dg::DG, cache)
    # to avoid a division by zero if the speed vanishes everywhere,
    # e.g. for steady-state linear advection
    max_scaled_speed = nextfloat(zero(t))

    @unpack contravariant_vectors, inverse_jacobian = cache.elements

    for element in eachelement(dg, cache)
        max_lambda1 = max_lambda2 = zero(max_scaled_speed)
        for j in eachnode(dg), i in eachnode(dg)
            u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
            lambda1, lambda2, lambda3 = max_abs_speeds(u_node, equations)

            # Local speeds transformed to the reference element
            Ja11, Ja12, Ja13 = Trixi.get_contravariant_vector(1, contravariant_vectors, i, j,
                                                  element)
            lambda1_transformed = abs(Ja11 * lambda1 + Ja12 * lambda2 + Ja13 * lambda3)
            Ja21, Ja22, Ja23 = Trixi.get_contravariant_vector(2, contravariant_vectors, i, j,
                                                  element)
            lambda2_transformed = abs(Ja21 * lambda1 + Ja22 * lambda2 + Ja23 * lambda3)

            inv_jacobian = abs(inverse_jacobian[i, j, element])

            max_lambda1 = max(max_lambda1, lambda1_transformed * inv_jacobian)
            max_lambda2 = max(max_lambda2, lambda2_transformed * inv_jacobian)
        end

        max_scaled_speed = max(max_scaled_speed, max_lambda1 + max_lambda2)
    end

    return 2 / (nnodes(dg) * max_scaled_speed)
end