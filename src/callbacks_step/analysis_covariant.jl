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

            u_exact = spherical2contravariant(initial_condition(x, t, equations),
                                              equations, cache.elements, i, j, element)

            u_numerical = Trixi.get_node_vars(u, equations, dg, i, j, element)

            diff = func(u_exact, equations) - func(u_numerical, equations)

            J = volume_element(cache.elements, i, j, element)

            l2_error += diff .^ 2 * (weights[i] * weights[j] * J)
            linf_error = @. max(linf_error, abs(diff))
            total_volume += weights[i] * weights[j] * J
        end
    end

    # For L2 error, divide by total volume
    l2_error = @. sqrt(l2_error / total_volume)

    return l2_error, linf_error
end

function Trixi.analyze(::typeof(Trixi.entropy_timederivative), du, u, t,
                       mesh::P4estMesh{2},
                       equations::AbstractCovariantEquations{2}, dg::DG, cache)
    # Calculate ∫(∂S/∂u ⋅ ∂u/∂t)dΩ
    Trixi.integrate_via_indices(u, mesh, equations, dg, cache,
                                du) do u, i, j, element, equations, dg, du
        u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
        du_node = Trixi.get_node_vars(du, equations, dg, i, j, element)
        dot(cons2entropy(u_node, equations, cache.elements, i, j, element), du_node)
    end
end

# Contribution to entropy time derivative from source term
function entropy_timederivative_source end
Trixi.pretty_form_utf(::typeof(entropy_timederivative_source)) = "∑∂S/∂U ⋅ src"
Trixi.pretty_form_ascii(::typeof(entropy_timederivative_source)) = "dsdu_dot_source"

function Trixi.analyze(::typeof(entropy_timederivative_source), du, u, t,
                       semi::Trixi.AbstractSemidiscretization)
    mesh, equations, dg, cache = Trixi.mesh_equations_solver_cache(semi)
    # Calculate ∫(∂S/∂u ⋅ source)dΩ
    Trixi.integrate_via_indices(u, mesh, equations, dg, cache,
                                du) do u, i, j, element, equations, dg, du
        u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
        x_node = Trixi.get_node_coords(cache.elements.node_coordinates, equations,
                                       dg,
                                       i, j, element)
        du_node = semi.source_terms(u_node, x_node, t, equations,
                                    cache.elements, i, j, element)
        dot(cons2entropy(u_node, equations, cache.elements, i, j, element), du_node)
    end
end
