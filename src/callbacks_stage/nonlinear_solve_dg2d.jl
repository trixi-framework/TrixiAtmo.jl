using Trixi: get_node_vars, @batch
using LinearAlgebra
using StaticArrays



@muladd begin

function nonlinear_solve_dg_2d!(u, residual, jacobian, variables_index_vector, tolerance,
                                equations::CompressibleRainyEulerEquations2D, dg::DGSEM, cache)
    max_iterations = 20

    # iterate over every DGSEM element
    @batch for element in eachelement(dg, cache)

        # iterate over every node
        for j in eachnode(dg), i in eachnode(dg)

            u_node = get_node_vars(u, equations, dg, i, j, element) 
            guess = SVector(u_node[7], u_node[8], u_node[9])

            # newton method
            for iteration in range(1, max_iterations)
                res_vector = residual(u_node, guess, equations)

                if (maximum(abs.(res_vector)) < tolerance)
                    break
                end

                jac_matrix = jacobian(u_node, guess, equations)
                guess += - jac_matrix \ res_vector

                #= warnings seem to have allocations...
                if iteration == max_iterations
                    @warn "newton method: tolerance not met"
                end
                =#
            end
            
            # similar to set_node_vars!
            for index in eachindex(variables_index_vector)
                u[variables_index_vector[index], i, j, element] = guess[index]
            end
        end
    end
end

end  # muladd end