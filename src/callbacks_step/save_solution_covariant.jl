# Convert to another set of variables using a solution_variables function that depends on 
# auxiliary variables
function convert_variables(u, solution_variables, mesh::P4estMesh{2},
                           equations, dg, cache)
    (; aux_node_vars) = cache.auxiliary_variables

    # Extract the number of solution variables to be output 
    # (may be different than the number of conservative variables) 
    n_vars = length(solution_variables(Trixi.get_node_vars(u, equations, dg, 1, 1, 1),
                                       get_node_aux_vars(aux_node_vars, equations, dg, 1, 1,
                                                         1), equations))
    # Allocate storage for output
    data = Array{eltype(u)}(undef, n_vars, nnodes(dg), nnodes(dg), nelements(dg, cache))

    # Loop over all nodes and convert variables, passing in auxiliary variables
    Trixi.@threaded for element in eachelement(dg, cache)
        for j in eachnode(dg), i in eachnode(dg)
            u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
            aux_node = get_node_aux_vars(aux_node_vars, equations, dg, i, j, element)
            u_node_converted = solution_variables(u_node, aux_node, equations)
            for v in 1:n_vars
                data[v, i, j, element] = u_node_converted[v]
            end
        end
    end
    return data
end

# Specialized save_solution_file method that supports a solution_variables function which 
# depends on auxiliary variables. The conversion must be defined as solution_variables(u, 
# aux_vars, equations), and an additional method must be defined as solution_variables(u,
# equations) = u, such that no conversion is done when auxiliary variables are not provided.
function Trixi.save_solution_file(u_ode, t, dt, iter,
    semi::SemidiscretizationHyperbolic{<:Trixi.AbstractMesh, <:AbstractCovariantEquations},
    solution_callback,
    element_variables = Dict{Symbol, Any}(),
    node_variables = Dict{Symbol, Any}();
    system = "")

    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    u = Trixi.wrap_array_native(u_ode, semi)

    # Perform the variable conversion at each node
    data = convert_variables(u, solution_callback.solution_variables, mesh, equations,
                             solver, cache)

    # Call the existing Trixi.save_solution_file, which will use solution_variables(u, 
    # equations). Since the variables are already converted above, we therefore require 
    # solution_variables(u, equations) = u.
    Trixi.save_solution_file(data, t, dt, iter, mesh, equations, solver, cache,
                       solution_callback, element_variables, 
                       node_variables, system = system)
end