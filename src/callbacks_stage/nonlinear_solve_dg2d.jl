using Trixi: get_node_vars, @batch, get_inverse_jacobian, set_node_vars!, get_node_coords, each_quad_node
using LinearAlgebra
using StaticArrays



@muladd begin

function nonlinear_solve_dg_2d!(u, residual, jacobian, variables_index_vector, tolerance,
                                equations::AbstractCompressibleRainyEulerEquations, dg::DGSEM, cache, mesh)
    max_iterations = 20
    rain_threshold = 1e-4
    #=threshold = 0.0
    @unpack weights = dg.basis
    @unpack inverse_jacobian = cache.elements=#
    
    # iterate over every DGSEM element
    @batch for element in eachelement(dg, cache)
        # iterate over every node
        for j in eachnode(dg), i in eachnode(dg)

            u_node = get_node_vars(u, equations, dg, i, j, element)
            guess = SVector(u_node[7], u_node[8], u_node[9])

            x_node = get_node_coords(cache.elements.node_coordinates, equations, dg, i, j, element)

            # keep rain positive
            if (u_node[3] < 0.0 || x_node[2] < 50.0)
                u[3, i, j, element] = 0.0
            end

            if (u_node[3] > rain_threshold)
                u[3, i, j, element] = rain_threshold * 0.1
            end

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


function nonlinear_solve_dg_2d!(u, residual, jacobian, variables_index_vector, tolerance,
                                equations::AbstractCompressibleRainyEulerEquations, dg::DGMulti, cache, mesh::DGMultiMesh)
    max_iterations = 20

    # iterate over every node
    @batch for element in eachelement(mesh, dg)
        for j in each_quad_node(mesh, dg)

            u_node = u[j, element]
            
            # keep rain positive
            if (u_node[3] < 0.0)
                u_node = SVector(u_node[1], u_node[2], 0.0, u_node[4], u_node[5], u_node[6],
                                        u_node[7], u_node[8], u_node[9])
            end

            guess  = SVector(u_node[7], u_node[8], u_node[9])

            # newton method
            for iteration in range(1, max_iterations)
                res_vector = residual(u_node, guess, equations)

                if (maximum(abs.(res_vector)) < tolerance)
                    break
                end

                jac_matrix = jacobian(u_node, guess, equations)
                guess += - jac_matrix \ res_vector
            end

            u[j, element] = SVector(u_node[1], u_node[2], u_node[3], u_node[4], u_node[5], u_node[6],
                             guess[1],  guess[2],  guess[3])
        end
    end
end

end  # muladd end