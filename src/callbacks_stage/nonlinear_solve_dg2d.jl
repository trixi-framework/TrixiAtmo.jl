using NLsolve
using Trixi: get_node_vars, @batch



@muladd begin

function nonlinear_solve_dg_2d!(u, residual, tolerance, variables_index_vector, equations, dg::DGSEM, cache)

    @batch for element in eachelement(dg, cache)

        for j in eachnode(dg), i in eachnode(dg)

            u_node = get_node_vars(u, equations, dg, i, j, element)
            initial_guess = u_node[variables_index_vector]
            nl_sol = nlsolve(residual(u_node, equations), initial_guess, ftol = tolerance, method=:newton).zero
            
            # similar to set_node_vars!
            for index in eachindex(variables_index_vector)
                u[variables_index_vector[index], i, j, element] = nl_sol[index]
            end
        end
    end
end

end  # muladd end