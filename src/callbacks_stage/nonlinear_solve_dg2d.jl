using NLsolve
using Trixi: get_node_vars


# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin

function nonlinear_solve_dg_2d!(u, residual, tolerance, variables_index_vector, equations, dg::DGSEM, cache)

    #TODO @threaded not defined
    for element in eachelement(dg, cache)
        
        initial_guess = similar(variables_index_vector)

        for j in eachnode(dg), i in eachnode(dg)

            u_node = get_node_vars(u, equations, dg, i, j, element)
            initial_guess = u_node[variables_index_vector]
            nl_sol = nlsolve(residual(u_node, equations), initial_guess, ftol = tolerance)
            
            # similar to set_node_vars!
            for index in eachindex(variables_index_vector)
                u[variables_index_vector[index], i, j, element] = nl_sol.zero[index]
            end
        end
    end
end

end  # muladd end