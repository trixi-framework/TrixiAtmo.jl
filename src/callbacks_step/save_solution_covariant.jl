@muladd begin
#! format: noindent

# Convert to another set of variables using a solution_variables function that depends on 
# auxiliary variables
function convert_variables(u, solution_variables, mesh::P4estMesh{2},
                           equations::AbstractCovariantEquations{2}, dg, cache)
    (; aux_node_vars) = cache.auxiliary_variables

    # Extract the number of solution variables to be output 
    # (may be different than the number of conservative variables) 
    n_vars = length(solution_variables(Trixi.get_node_vars(u, equations, dg, 1, 1, 1),
                                       get_node_aux_vars(aux_node_vars, equations, dg,
                                                         1, 1,
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
end # muladd
