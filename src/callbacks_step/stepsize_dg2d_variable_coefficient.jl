@muladd begin

function Trixi.max_dt(u, t,
    mesh::P4estMesh{2},
    constant_speed::False,
    equations::AbstractVariableCoefficientEquations{2},
    dg::DG,
    cache,
)
    # to avoid a division by zero if the speed vanishes everywhere,
    # e.g. for steady-state linear advection
    max_scaled_speed = nextfloat(zero(t))

    (; contravariant_vectors, inverse_jacobian) = cache.elements
    (; aux_node_vars) = cache.aux_vars

    Trixi.@batch reduction = (max, max_scaled_speed) for element in eachelement(dg, cache)
        max_lambda1 = max_lambda2 = zero(max_scaled_speed)
        for j in eachnode(dg), i in eachnode(dg)
            u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
            aux_node = get_node_aux_vars(aux_node_vars, equations, dg, i, j, element)
            lambda1, lambda2 = max_abs_speeds(u_node, aux_node, equations)
        
            # Local speeds transformed to the reference element
            Ja11, Ja12 = Trixi.get_contravariant_vector(1, contravariant_vectors, i, j, element)
            lambda1_transformed = abs(Ja11 * lambda1 + Ja12 * lambda2)
            Ja21, Ja22 = Trixi.get_contravariant_vector(2, contravariant_vectors, i, j, element)
            lambda2_transformed = abs(Ja21 * lambda1 + Ja22 * lambda2)

            inv_jacobian = abs(inverse_jacobian[i, j, element])

            max_lambda1 = max(max_lambda1, lambda1_transformed * inv_jacobian)
            max_lambda2 = max(max_lambda2, lambda2_transformed * inv_jacobian)
        end
        
        max_scaled_speed = max(max_scaled_speed, max_lambda1 + max_lambda2)
    end
    
    return 2 / (nnodes(dg) * max_scaled_speed)
end

function Trixi.analyze(::typeof(Trixi.entropy_timederivative), du, u, t,
                 mesh::Union{TreeMesh{2}, StructuredMesh{2}, StructuredMeshView{2},
                             UnstructuredMesh2D, P4estMesh{2}, T8codeMesh{2}},
                 equations::PerturbationEulerEquations2DAuxVars, dg::DG, cache)
    # Calculate ∫(∂S/∂u ⋅ ∂u/∂t)dΩ
    (; aux_node_vars) = cache.aux_vars
    Trixi.integrate_via_indices(u, mesh, equations, dg, cache,
                          du) do u, i, j, element, equations, dg, du
        u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
        aux_node = get_node_aux_vars(aux_node_vars, equations, dg, i, j, element)
        du_node = Trixi.get_node_vars(du, equations, dg, i, j, element)
        dot(cons2entropy(u_node, aux_node, equations), du_node)
    end
end

end #@muladd