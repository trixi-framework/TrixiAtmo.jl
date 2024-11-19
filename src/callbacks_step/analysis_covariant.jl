@muladd begin
#! format: noindent

function Trixi.analyze(::typeof(Trixi.entropy_timederivative), du, u, t,
                       mesh::P4estMesh{2},
                       equations::AbstractCovariantEquations{2}, dg::DG, cache)
    # Calculate ∫(∂S/∂u ⋅ ∂u/∂t)dΩ
    Trixi.integrate_via_indices(u, mesh, equations, dg, cache,
                                du) do u, i, j, element, equations, dg, du
        aux_vars_node = get_node_aux_vars(cache.elements.auxiliary_variables, equations,
                                          dg, i, j, element)
        u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
        du_node = Trixi.get_node_vars(du, equations, dg, i, j, element)
        dot(cons2entropy(u_node, aux_vars_node, equations), du_node)
    end
end

function Trixi.calc_error_norms(func, u, t, analyzer, mesh::P4estMesh{2},
                                equations::AbstractCovariantEquations{2},
                                initial_condition, dg::DGSEM, cache, cache_analysis)
    (; weights) = dg.basis
    (; node_coordinates) = cache.elements

    # Set up data structures
    l2_error = zero(func(Trixi.get_node_vars(u, equations, dg, 1, 1, 1), equations))
    linf_error = copy(l2_error)
    total_volume = zero(real(mesh))

    # Iterate over all elements for error calculations
    for element in eachelement(dg, cache)

        # Calculate errors at each volume quadrature node
        for j in eachnode(dg), i in eachnode(dg)
            x = Trixi.get_node_coords(node_coordinates, equations, dg, i, j, element)

            aux_vars_node = get_node_aux_vars(cache.elements.auxiliary_variables, equations,
                                              dg, i, j, element)

            u_exact = initial_condition(x, t, aux_vars_node, equations)

            u_numerical = Trixi.get_node_vars(u, equations, dg, i, j, element)

            diff = func(u_exact, equations) - func(u_numerical, equations)

            sqrtG = volume_element(aux_vars_node, equations)

            l2_error += diff .^ 2 * (weights[i] * weights[j] * sqrtG)
            linf_error = @. max(linf_error, abs(diff))
            total_volume += weights[i] * weights[j] * sqrtG
        end
    end

    # For L2 error, divide by total volume
    l2_error = @. sqrt(l2_error / total_volume)

    return l2_error, linf_error
end
end # muladd
