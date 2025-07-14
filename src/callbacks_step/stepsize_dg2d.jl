@muladd begin
#! format: noindent

# Specialization of max_dt function for 3D equations in 2D manifolds
function Trixi.max_dt(u, t,
                      mesh::Union{StructuredMesh{2}, UnstructuredMesh2D, P4estMesh{2},
                                  T8codeMesh{2}},
                      constant_speed::False, equations::AbstractEquations{3}, dg::DG,
                      cache)
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
            Ja11, Ja12, Ja13 = Trixi.get_contravariant_vector(1, contravariant_vectors,
                                                              i,
                                                              j,
                                                              element)
            lambda1_transformed = abs(Ja11 * lambda1 + Ja12 * lambda2 + Ja13 * lambda3)
            Ja21, Ja22, Ja23 = Trixi.get_contravariant_vector(2, contravariant_vectors,
                                                              i,
                                                              j,
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

# Specialization of max_dt function for covariant formulation on 2D manifolds
function Trixi.max_dt(u, t, mesh::P4estMesh{2}, constant_speed::False,
                      equations::AbstractCovariantEquations{2},
                      dg::DG, cache)
    (; aux_node_vars) = cache.auxiliary_variables

    # to avoid a division by zero if the speed vanishes everywhere,
    # e.g. for steady-state linear advection
    max_scaled_speed = nextfloat(zero(t))

    # Because the covariant form computes max_abs_speeds using the contravariant 
    # velocity components already, there is no need to transform them here
    for element in eachelement(dg, cache)
        max_lambda1 = max_lambda2 = zero(max_scaled_speed)
        for j in eachnode(dg), i in eachnode(dg)
            u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
            aux_node = get_node_aux_vars(aux_node_vars, equations, dg, i, j, element)
            lambda1, lambda2 = Trixi.max_abs_speeds(u_node, aux_node, equations)
            max_lambda1 = max(max_lambda1, lambda1)
            max_lambda2 = max(max_lambda2, lambda2)
        end

        max_scaled_speed = max(max_scaled_speed, max_lambda1 + max_lambda2)
    end
    return 2 / (nnodes(dg) * max_scaled_speed)
end

function Trixi.max_dt(u, t, mesh::DGMultiMesh,
                constant_speed::False, equations::AbstractCovariantEquations{NDIMS},
                dg::DGMulti, cache) where {NDIMS}
    @unpack md = mesh
    rd = dg.basis
    (; aux_values) = cache

    dt_min = Inf
    for e in eachelement(mesh, dg, cache)
        max_speeds = ntuple(_ -> nextfloat(zero(t)), NDIMS)
        for i in 1:nnodes(dg)
            u_node = u[i, e]
            aux_node = aux_values[i, e]
            detg = area_element(aux_node, equations)
            lambda_i = max_abs_speeds(u_node, aux_node, equations)
            max_speeds = max.(max_speeds, lambda_i)
        end
        dt_min = min(dt_min, 1 / sum(max_speeds))
    end

    # This mimics `max_dt` for `TreeMesh`, except that `nnodes(dg)` is replaced by
    # `polydeg+1`. This is because `nnodes(dg)` returns the total number of
    # multi-dimensional nodes for DGMulti solver types, while `nnodes(dg)` returns
    # the number of 1D nodes for `DGSEM` solvers.
    return 2 * dt_min * Trixi.dt_polydeg_scaling(dg)
end
end # @muladd
