# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Variables of the vertical Jacobian at one node: the unknowns `(rho, rho_w, rho_theta)` of
# the linear systems of the columns, laid out `[column, v, node, layer]`, where `rho_w` is
# the component of the momentum in the unit contravariant direction `vertical_direction` of
# the vertical reference coordinate (see `compute_vertical_direction_and_norm`). Equations
# support the vertical Jacobian by implementing `cons2hevi`, `hevi2cons`, and
# `hevi_pressure_derivative`.

@inline function get_column_node_vars(u_vertical, column, node, layer)
    return SVector(u_vertical[column, 1, node, layer],
                   u_vertical[column, 2, node, layer],
                   u_vertical[column, 3, node, layer])
end

@inline function set_column_node_vars!(u_vertical, u_node, column, node, layer)
    u_vertical[column, 1, node, layer] = u_node[1]
    u_vertical[column, 2, node, layer] = u_node[2]
    u_vertical[column, 3, node, layer] = u_node[3]
    return nothing
end

# Unknowns of the vertical Jacobian of the conservative variables `u`. The map is linear,
# so that it also maps the right-hand side of the linear system.
@inline function cons2hevi(u, vertical_direction,
                           equations::CompressibleEulerPotentialTemperatureEquations3D)
    rho_v = SVector(u[2], u[3], u[4])
    return SVector(u[1], dot(vertical_direction, rho_v), u[5])
end

# Conservative variables of the solution of the linear system from the solution
# `x_vertical` of the vertical Jacobian. The horizontal components of the momentum are not
# unknowns of the vertical Jacobian and are solved with `J = 0`, i.e., their right-hand side
# `b` is multiplied by `gamma_dt`.
@inline function hevi2cons(x_vertical, b, vertical_direction, gamma_dt,
                           equations::CompressibleEulerPotentialTemperatureEquations3D)
    rho, rho_w, rho_theta = x_vertical
    rho_v = SVector(b[2], b[3], b[4])
    rho_v = gamma_dt * rho_v +
            (rho_w - gamma_dt * dot(vertical_direction, rho_v)) * vertical_direction
    return SVector(rho, rho_v[1], rho_v[2], rho_v[3], rho_theta)
end

# Potential temperature and derivative of the pressure `p = K rho_theta^gamma` with respect
# to `rho_theta` of the unknowns `u_vertical_node` of the vertical Jacobian
@inline function hevi_pressure_derivative(u_vertical_node,
                                          equations::CompressibleEulerPotentialTemperatureEquations3D)
    rho, _, rho_theta = u_vertical_node
    theta = rho_theta / rho
    dp_drho_theta = (equations.gamma * equations.K) * rho_theta^(equations.gamma - 1)
    return theta, dp_drho_theta
end

# Copy the unknowns of the vertical Jacobian of the node `(i, j, k, element)` of `u`,
# multiplied by `scaling`, to the columns. The `scaling` is one for the state and the norm
# `n` for the right-hand side of the linear system, i.e., the diagonal similarity
# transformation `N` of the vertical operator.
Base.@propagate_inbounds function copy_node_to_columns!(u_vertical, u, scaling,
                                                        vertical_direction_and_norm,
                                                        horizontal_index_of_element,
                                                        vertical_index_of_element,
                                                        nnodes_, equations, solver,
                                                        i, j, k, element)
    column = column_index(horizontal_index_of_element, nnodes_, i, j, element)
    layer = vertical_index_of_element[element]
    vertical_direction = get_vertical_direction(vertical_direction_and_norm,
                                                i, j, k, element)
    u_node = get_node_vars(u, equations, solver, i, j, k, element)
    set_column_node_vars!(u_vertical,
                          scaling * cons2hevi(u_node, vertical_direction, equations),
                          column, k, layer)
    return nothing
end

# Copy the solution of the linear system, divided by the norm `n`, to `b`
Base.@propagate_inbounds function copy_node_from_columns!(b, u_vertical,
                                                          vertical_direction_and_norm,
                                                          horizontal_index_of_element,
                                                          vertical_index_of_element,
                                                          nnodes_, gamma_dt, equations,
                                                          solver, i, j, k, element)
    column = column_index(horizontal_index_of_element, nnodes_, i, j, element)
    layer = vertical_index_of_element[element]
    vertical_direction = get_vertical_direction(vertical_direction_and_norm,
                                                i, j, k, element)
    inv_norm_Ja3 = inv(vertical_direction_and_norm[4, i, j, k, element])
    x_vertical = inv_norm_Ja3 * get_column_node_vars(u_vertical, column, k, layer)
    b_node = get_node_vars(b, equations, solver, i, j, k, element)
    set_node_vars!(b,
                   hevi2cons(x_vertical, b_node, vertical_direction, gamma_dt,
                             equations),
                   equations, solver, i, j, k, element)
    return nothing
end

# Unit contravariant direction `nu` of the vertical reference coordinate at the node
# `(i, j, k, element)`
@inline function get_vertical_direction(vertical_direction_and_norm, i, j, k, element)
    return SVector(vertical_direction_and_norm[1, i, j, k, element],
                   vertical_direction_and_norm[2, i, j, k, element],
                   vertical_direction_and_norm[3, i, j, k, element])
end
end # @muladd
