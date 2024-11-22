@muladd begin
#! format: noindent

# Entropy time derivative which uses auxiliary variables
function Trixi.analyze(::typeof(Trixi.entropy_timederivative), du, u, t,
                       mesh::P4estMesh{2},
                       equations::AbstractCovariantEquations{2}, dg::DG, cache)
    (; aux_node_vars) = cache.auxiliary_variables

    # Calculate ∫(∂S/∂u ⋅ ∂u/∂t)dΩ
    Trixi.integrate_via_indices(u, mesh, equations, dg, cache,
                                du) do u, i, j, element, equations, dg, du_node
        # Get auxiliary variables, solution variables, and time derivative at given node
        a_node = get_node_aux_vars(aux_node_vars, equations, dg, i, j, element)
        u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
        du_node = Trixi.get_node_vars(du, equations, dg, i, j, element)

        # compute ∂S/∂u ⋅ ∂u/∂t, where the entropy variables ∂S/∂u depend on the solution 
        # and auxiliary variables
        dot(cons2entropy(u_node, a_node, equations), du_node)
    end
end

# L2 and Linf error calculation for the covariant form
function Trixi.calc_error_norms(func, u, t, analyzer, mesh::P4estMesh{2},
                                equations::AbstractCovariantEquations{2},
                                initial_condition, dg::DGSEM, cache, cache_analysis)
    (; weights) = dg.basis
    (; node_coordinates) = cache.elements
    (; aux_node_vars) = cache.auxiliary_variables

    # Set up data structures
    l2_error = zero(func(Trixi.get_node_vars(u, equations, dg, 1, 1, 1), equations))
    linf_error = copy(l2_error)
    total_volume = zero(real(mesh))

    # Iterate over all elements for error calculations
    for element in eachelement(dg, cache)

        # Calculate errors at each volume quadrature node
        for j in eachnode(dg), i in eachnode(dg)
            x_node = Trixi.get_node_coords(node_coordinates, equations, dg, i, j,
                                           element)

            # Convert exact solution into contravariant components using geometric
            # information stored in aux vars
            a_node = get_node_aux_vars(aux_node_vars, equations, dg, i, j, element)
            u_exact = initial_condition(x_node, t, a_node, equations)

            # Compute the difference as usual
            u_numerical = Trixi.get_node_vars(u, equations, dg, i, j, element)
            diff = func(u_exact, equations) - func(u_numerical, equations)

            # For the L2 error, integrate with respect to volume element stored in aux vars 
            sqrtG = volume_element(a_node, equations)
            l2_error += diff .^ 2 * (weights[i] * weights[j] * sqrtG)

            # Compute Linf error as usual
            linf_error = @. max(linf_error, abs(diff))

            # Increment total volume according to the volume element stored in aux vars
            total_volume += weights[i] * weights[j] * sqrtG
        end
    end

    # For L2 error, divide by total volume
    l2_error = @. sqrt(l2_error / total_volume)

    return l2_error, linf_error
end
end # muladd
