
# This method is called to determine whether the callback should be activated
function mean_temperature_callback(u, t, integrator)
    # TODO: Here we could do a check and only activate this periodically
    return true
end

# Actual callback function
function mean_temperature_callback(integrator)
    semi = integrator.p
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    u = Trixi.wrap_array(integrator.u, semi)

    update_mean_temperature!(u, mesh, equations, solver, cache)
    
    return nothing
end

# This function computes the mean temperature R*T in each element and stores it in the last variable
# of the solution vector
function update_mean_temperature!(u, mesh::Trixi.AbstractMesh{2}, equations, solver, cache)
    @unpack weights = solver.basis
    @unpack inverse_jacobian = cache.elements

    Trixi.@threaded for element in eachelement(solver, cache)
        # compute the element local mean temperature
        RT_mean = zero(eltype(u))
        total_volume = zero(eltype(u))
        for j in eachnode(solver), i in eachnode(solver)
            volume_jacobian = abs(inv(Trixi.get_inverse_jacobian(inverse_jacobian,
                                                                    mesh,
                                                                    i, j, element)))
            u_node = get_node_vars(u, equations, solver, i, j, element)
            p_node = pressure(u_node, equations)
            rho_node = u_node[1]
            RT_node = p_node / rho_node

            RT_mean += RT_node * weights[i] * weights[j] * volume_jacobian
             
            total_volume += weights[i] * weights[j] * volume_jacobian
        end
        # normalize with the total volume
        RT_mean = RT_mean / total_volume

        # Update the solution vector
        for j in eachnode(solver), i in eachnode(solver)
            u[end,i,j,element] = RT_mean
        end
    end

    return nothing
end

"""
    MeanTemperatureCallback()

    Callback to compute the mean temperature `R*T` in each element and store it in the last variable
    of the solution vector.
"""
function MeanTemperatureCallback()
    return Trixi.DiscreteCallback(mean_temperature_callback, mean_temperature_callback)
end