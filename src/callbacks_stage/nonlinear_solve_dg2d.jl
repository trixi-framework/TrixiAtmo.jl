using NonlinearSolve
using Trixi: get_node_vars

@muladd begin
function nonlinear_solve_dg_2d!(u, nonlin_fct, variables_index_vector,
                                tolerance, equations, dg::DGSEM, cache)
    for element in eachelement(dg, cache)
        for j in eachnode(dg), i in eachnode(dg)
            u_node = get_node_vars(u, equations, dg, i, j, element)
            initial_guess = u_node[variables_index_vector]
            prob = NonlinearProblem(nonlin_fct, initial_guess, (u_node, equations))
            # TODO: abs and rel?
            res = solve(prob, SimpleNewtonRaphson(); abstol = tolerance, reltol = tolerance)

            # similar to set_node_vars!
            for index in eachindex(variables_index_vector)
                u[variables_index_vector[index], i, j, element] = res[index]
            end
        end
    end
end
end  # muladd end
