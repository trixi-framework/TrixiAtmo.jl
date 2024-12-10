@muladd begin
#! format: noindent

# For the covariant form, we want to integrate using the exact area element 
# √G = (det(AᵀA))^(1/2), which is stored in cache.auxiliary_variables, not the approximate 
# area element used in the Cartesian formulation, which stored in cache.elements.
function Trixi.integrate_via_indices(func::Func, u,
                                     mesh::Union{StructuredMesh{2},
                                                 StructuredMeshView{2},
                                                 UnstructuredMesh2D, P4estMesh{2},
                                                 T8codeMesh{2}},
                                     equations::AbstractCovariantEquations{2},
                                     dg::DGSEM, cache, args...;
                                     normalize = true) where {Func}
    (; weights) = dg.basis
    (; aux_node_vars) = cache.auxiliary_variables

    # Initialize integral with zeros of the right shape
    integral = zero(func(u, 1, 1, 1, equations, dg, args...))
    total_volume = zero(real(mesh))

    # Use quadrature to numerically integrate over entire domain
    for element in eachelement(dg, cache)
        for j in eachnode(dg), i in eachnode(dg)
            aux_node = get_node_aux_vars(aux_node_vars, equations, dg, i, j, element)
            sqrtG = area_element(aux_node, equations)
            integral += weights[i] * weights[j] * sqrtG *
                        func(u, i, j, element, equations, dg, args...)
            total_volume += weights[i] * weights[j] * sqrtG
        end
    end

    # Normalize with total volume
    if normalize
        integral = integral / total_volume
    end

    return integral
end

# Entropy time derivative for cons2entropy function which depends on auxiliary variables
function Trixi.analyze(::typeof(Trixi.entropy_timederivative), du, u, t,
                       mesh::P4estMesh{2},
                       equations::AbstractCovariantEquations{2}, dg::DG, cache)
    (; aux_node_vars) = cache.auxiliary_variables

    # Calculate ∫(∂S/∂u ⋅ ∂u/∂t)dΩ
    Trixi.integrate_via_indices(u, mesh, equations, dg, cache,
                                du) do u, i, j, element, equations, dg, du_node
        # Get auxiliary variables, solution variables, and time derivative at given node
        aux_node = get_node_aux_vars(aux_node_vars, equations, dg, i, j, element)
        u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
        du_node = Trixi.get_node_vars(du, equations, dg, i, j, element)

        # compute ∂S/∂u ⋅ ∂u/∂t, where the entropy variables ∂S/∂u depend on the solution 
        # and auxiliary variables
        dot(cons2entropy(u_node, aux_node, equations), du_node)
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
            aux_node = get_node_aux_vars(aux_node_vars, equations, dg, i, j, element)
            u_exact = initial_condition(x_node, t, aux_node, equations)

            # Compute the difference as usual
            u_numerical = Trixi.get_node_vars(u, equations, dg, i, j, element)
            diff = func(u_exact, equations) - func(u_numerical, equations)

            # For the L2 error, integrate with respect to area element stored in aux vars 
            sqrtG = area_element(aux_node, equations)
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
