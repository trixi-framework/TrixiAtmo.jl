using Trixi: get_node_vars, @batch, get_inverse_jacobian, set_node_vars!, get_node_coords
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
        #= positivity preserving limiter zhang shu test determine minimum value
        value_min = typemax(eltype(u))
        for j in eachnode(dg), i in eachnode(dg)
            u_node = get_node_vars(u, equations, dg, i, j, element)
            value_min = min(value_min, u_node[3])
        end
        
        # detect if limiting is necessary
        if (value_min < threshold)
            # compute mean value
            u_mean = zero(get_node_vars(u, equations, dg, 1, 1, element))
            total_volume = zero(eltype(u))
            for j in eachnode(dg), i in eachnode(dg)
                volume_jacobian = abs(inv(get_inverse_jacobian(inverse_jacobian, mesh,
                                                           i, j, element)))
                u_node = get_node_vars(u, equations, dg, i, j, element)
                u_mean += u_node * weights[i] * weights[j] * volume_jacobian
                total_volume += weights[i] * weights[j] * volume_jacobian
            end
            
            # normalize with the total volume
            u_mean = u_mean / total_volume

            # We compute the value directly with the mean values, as we assume that
            # Jensen's inequality holds (e.g. pressure for compressible Euler equations).
            value_mean = u_mean[3]
            theta = (value_mean - threshold) / (value_mean - value_min)
            for j in eachnode(dg), i in eachnode(dg)
                u_node = get_node_vars(u, equations, dg, i, j, element)
                set_node_vars!(u, theta * u_node + (1 - theta) * u_mean,
                               equations, dg, i, j, element)
            end
        end
        
        # positivity preserving limiter zhang shu test determine minimum value
        value_min = typemax(eltype(u))
        for j in eachnode(dg), i in eachnode(dg)
            u_node = get_node_vars(u, equations, dg, i, j, element)
            value_min = min(value_min, u_node[2])
        end
        
        # detect if limiting is necessary
        if (value_min < 1e-5)
            # compute mean value
            u_mean = zero(get_node_vars(u, equations, dg, 1, 1, element))
            total_volume = zero(eltype(u))
            for j in eachnode(dg), i in eachnode(dg)
                volume_jacobian = abs(inv(get_inverse_jacobian(inverse_jacobian, mesh,
                                                           i, j, element)))
                u_node = get_node_vars(u, equations, dg, i, j, element)
                u_mean += u_node * weights[i] * weights[j] * volume_jacobian
                total_volume += weights[i] * weights[j] * volume_jacobian
            end
            
            # normalize with the total volume
            u_mean = u_mean / total_volume

            # We compute the value directly with the mean values, as we assume that
            # Jensen's inequality holds (e.g. pressure for compressible Euler equations).
            value_mean = u_mean[2]
            theta = (value_mean - 1e-5) / (value_mean - value_min)
            for j in eachnode(dg), i in eachnode(dg)
                u_node = get_node_vars(u, equations, dg, i, j, element)
                set_node_vars!(u, theta * u_node + (1 - theta) * u_mean,
                               equations, dg, i, j, element)
            end
        end

        # positivity preserving limiter zhang shu test determine minimum value
        value_min = typemax(eltype(u))
        for j in eachnode(dg), i in eachnode(dg)
            u_node = get_node_vars(u, equations, dg, i, j, element)
            value_min = min(value_min, u_node[1])
        end
        
        # detect if limiting is necessary
        if (value_min < 1e-5)
            # compute mean value
            u_mean = zero(get_node_vars(u, equations, dg, 1, 1, element))
            total_volume = zero(eltype(u))
            for j in eachnode(dg), i in eachnode(dg)
                volume_jacobian = abs(inv(get_inverse_jacobian(inverse_jacobian, mesh,
                                                           i, j, element)))
                u_node = get_node_vars(u, equations, dg, i, j, element)
                u_mean += u_node * weights[i] * weights[j] * volume_jacobian
                total_volume += weights[i] * weights[j] * volume_jacobian
            end
            
            # normalize with the total volume
            u_mean = u_mean / total_volume

            # We compute the value directly with the mean values, as we assume that
            # Jensen's inequality holds (e.g. pressure for compressible Euler equations).
            value_mean = u_mean[1]
            theta = (value_mean - 1e-5) / (value_mean - value_min)
            for j in eachnode(dg), i in eachnode(dg)
                u_node = get_node_vars(u, equations, dg, i, j, element)
                set_node_vars!(u, theta * u_node + (1 - theta) * u_mean,
                               equations, dg, i, j, element)
            end
        end =#

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

end  # muladd end