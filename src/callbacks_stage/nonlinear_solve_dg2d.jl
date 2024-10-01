#using NLsolve
using Trixi: get_node_vars, @batch
using LinearAlgebra
using StaticArrays



### Warning: Hardcoded for performance reasons.
### will not work with differently structured problems

@muladd begin

function nonlinear_solve_dg_2d!(u, residual, jacobian, variables_index_vector, tolerance, equations, dg::DGSEM, cache)

    max_iterations = 50

    # iterate over every DGSEM element
    @batch for element in eachelement(dg, cache)
        # allocate static arrays for every thread
        res_vector = SVector{3, Real}
        jac_matrix = SMatrix{3, 3, Real}

        # iterate over every node
        for j in eachnode(dg), i in eachnode(dg)

            u_node = get_node_vars(u, equations, dg, i, j, element)
            initial_guess = SVector(u_node[7], u_node[8], u_node[9])

            # define residual function and jacobian function at the given node
            res = residual(u_node, equations)
            jac = jacobian(u_node, equations)

            # newton method
            for iteration in range(1, max_iterations)
                res_vector = res(res_vector, initial_guess)

                if (maximum(abs.(res_vector)) < tolerance)
                    break
                end

                jac_matrix = jac(jac_matrix, initial_guess)
                initial_guess += - jac_matrix \ res_vector

                #= warnings seem to have allocations...
                if iteration == max_iterations
                    @warn "newton method: tolerance not met"
                end
                =#
            end

            #= Old slow solver
            nl_sol = nlsolve(residual(u_node, equations),
                             jacobian(u_node, equations),
                             initial_guess, ftol = tolerance, method = :newton).zero
            =#
            
            # similar to set_node_vars!
            for index in eachindex(variables_index_vector)
                u[variables_index_vector[index], i, j, element] = initial_guess[index]
            end
        end
    end
end

end  # muladd end