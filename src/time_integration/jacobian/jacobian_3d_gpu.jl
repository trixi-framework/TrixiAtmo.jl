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
# support the vertical Jacobian by implementing `cons2hevi`, `cons2hevi_b`, `hevi2cons`,
# and `hevi_pressure_derivative`. The Jacobian is the one of the potential temperature
# form; for other equations, the right-hand side and the solution of the linear system are
# transformed with the linearization `T = d(cons) / d(rho, rho_w, rho_theta, ...)` of the
# change of variables at the state `u` from which the Jacobian is assembled, i.e., the
# linear system `(I / gamma_dt - T J T^{-1}) x = b` is solved.

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

# Equations supported by the vertical Jacobian
const EquationsJacobianVertical = Union{CompressibleEulerPotentialTemperatureEquations3D,
                                        CompressibleEulerMoistPotentialTemperatureEquations3D,
                                        CompressibleEulerEnergyEquations3D,
                                        CompressibleEulerVectorInvariantEquations3D}

# Potential temperature and derivative of the pressure `p = K rho_theta^gamma` with respect
# to `rho_theta` of the unknowns `u_vertical_node` of the vertical Jacobian
@inline function hevi_pressure_derivative(u_vertical_node,
                                          equations::EquationsJacobianVertical)
    rho, _, rho_theta = u_vertical_node
    theta = rho_theta / rho
    dp_drho_theta = (equations.gamma * equations.K) * rho_theta^(equations.gamma - 1)
    return theta, dp_drho_theta
end

###############################################################################
# Equations with the potential temperature (density) as prognostic variable and the
# momentum, i.e., with the first variables `(rho, rho_v1, rho_v2, rho_v3, rho_theta)`, or
# `rho_theta_v` for the `CompressibleEulerMoistPotentialTemperatureEquations3D`, on which
# the pressure depends as on `rho_theta` for dry air. The change of variables is linear.
const PotentialTemperatureEquations3D = Union{CompressibleEulerPotentialTemperatureEquations3D,
                                              CompressibleEulerMoistPotentialTemperatureEquations3D}

# Unknowns of the vertical Jacobian of the state `u`
@inline function cons2hevi(u, vertical_direction,
                           equations::PotentialTemperatureEquations3D)
    rho_v = SVector(u[2], u[3], u[4])
    return SVector(u[1], dot(vertical_direction, rho_v), u[5])
end

# Right-hand side of the linear system of the vertical Jacobian for the right-hand side `b`
# of the linear system of the conservative variables, linearized at the state `u`
@inline function cons2hevi_b(b, u, vertical_direction,
                             equations::PotentialTemperatureEquations3D)
    return cons2hevi(b, vertical_direction, equations)
end

# Solution of the linear system of the conservative variables from the solution
# `x_vertical` of the vertical Jacobian, linearized at the state `u`. The variables that
# are not unknowns of the vertical Jacobian, i.e., the horizontal components of the
# momentum and the moisture variables, are solved with `J = 0`, i.e., their right-hand side
# `b` is multiplied by `gamma_dt`.
@inline function hevi2cons(x_vertical, b, u, vertical_direction, gamma_dt,
                           equations::PotentialTemperatureEquations3D)
    rho, rho_w, rho_theta = x_vertical
    rho_v = SVector(b[2], b[3], b[4])
    rho_v = gamma_dt * rho_v +
            (rho_w - gamma_dt * dot(vertical_direction, rho_v)) * vertical_direction
    return SVector(rho, rho_v[1], rho_v[2], rho_v[3], rho_theta,
                   ntuple(v -> gamma_dt * b[5 + v], Val(nvariables(equations) - 5))...)
end

###############################################################################
# `CompressibleEulerVectorInvariantEquations3D` with the variables
# `(rho, v1, v2, v3, rho_theta)` and the vertical momentum `rho_w = rho (nu . v)`

@inline function cons2hevi(u, vertical_direction,
                           equations::CompressibleEulerVectorInvariantEquations3D)
    rho = u[1]
    v = SVector(u[2], u[3], u[4])
    return SVector(rho, rho * dot(vertical_direction, v), u[5])
end

# `d(rho_w) = (nu . v) d(rho) + rho nu . dv`
@inline function cons2hevi_b(b, u, vertical_direction,
                             equations::CompressibleEulerVectorInvariantEquations3D)
    rho = u[1]
    v = SVector(u[2], u[3], u[4])
    b_v = SVector(b[2], b[3], b[4])
    return SVector(b[1],
                   dot(vertical_direction, v) * b[1] +
                   rho * dot(vertical_direction, b_v),
                   b[5])
end

# The vertical component of the velocity follows from `d(rho_w)`, while its horizontal
# components are solved with `J = 0`
@inline function hevi2cons(x_vertical, b, u, vertical_direction, gamma_dt,
                           equations::CompressibleEulerVectorInvariantEquations3D)
    x_rho, x_rho_w, x_rho_theta = x_vertical
    rho = u[1]
    v = SVector(u[2], u[3], u[4])
    b_v = SVector(b[2], b[3], b[4])
    x_w = (x_rho_w - dot(vertical_direction, v) * x_rho) / rho
    x_v = gamma_dt * b_v +
          (x_w - gamma_dt * dot(vertical_direction, b_v)) * vertical_direction
    return SVector(x_rho, x_v[1], x_v[2], x_v[3], x_rho_theta)
end

###############################################################################
# `CompressibleEulerEnergyEquations3D` with the variables `(rho, rho_v1, rho_v2, rho_v3, E)`
# of the energy `E = p / (gamma - 1) + |rho_v|^2 / (2 rho)` and the potential temperature
# density `rho_theta = (p / K)^(1 / gamma)` with `K = p_0 (R / p_0)^gamma`

# Velocity, pressure, and potential temperature density of the state `u`
@inline function hevi_velocity_pressure_rho_theta(u,
                                                  equations::CompressibleEulerEnergyEquations3D)
    rho = u[1]
    v = SVector(u[2], u[3], u[4]) / rho
    p = (equations.gamma - 1) * (u[5] - 0.5f0 * rho * dot(v, v))
    rho_theta = (p / equations.K)^inv(equations.gamma)
    return v, p, rho_theta
end

@inline function cons2hevi(u, vertical_direction,
                           equations::CompressibleEulerEnergyEquations3D)
    _, _, rho_theta = hevi_velocity_pressure_rho_theta(u, equations)
    rho_v = SVector(u[2], u[3], u[4])
    return SVector(u[1], dot(vertical_direction, rho_v), rho_theta)
end

# `d(rho_theta) = rho_theta / (gamma p) dp` with
# `dp = (gamma - 1) (|v|^2 / 2 d(rho) - v . d(rho_v) + dE)`
@inline function cons2hevi_b(b, u, vertical_direction,
                             equations::CompressibleEulerEnergyEquations3D)
    (; gamma) = equations
    v, p, rho_theta = hevi_velocity_pressure_rho_theta(u, equations)
    b_rho_v = SVector(b[2], b[3], b[4])
    dp = (gamma - 1) * (0.5f0 * dot(v, v) * b[1] - dot(v, b_rho_v) + b[5])
    return SVector(b[1], dot(vertical_direction, b_rho_v), rho_theta / (gamma * p) * dp)
end

# `dE = gamma p / ((gamma - 1) rho_theta) d(rho_theta) - |v|^2 / 2 d(rho) + v . d(rho_v)`,
# where the horizontal components of `d(rho_v)` are solved with `J = 0`
@inline function hevi2cons(x_vertical, b, u, vertical_direction, gamma_dt,
                           equations::CompressibleEulerEnergyEquations3D)
    (; gamma) = equations
    x_rho, x_rho_w, x_rho_theta = x_vertical
    v, p, rho_theta = hevi_velocity_pressure_rho_theta(u, equations)
    b_rho_v = SVector(b[2], b[3], b[4])
    x_rho_v = gamma_dt * b_rho_v +
              (x_rho_w - gamma_dt * dot(vertical_direction, b_rho_v)) *
              vertical_direction
    x_E = gamma * p / ((gamma - 1) * rho_theta) * x_rho_theta -
          0.5f0 * dot(v, v) * x_rho + dot(v, x_rho_v)
    return SVector(x_rho, x_rho_v[1], x_rho_v[2], x_rho_v[3], x_E)
end

###############################################################################
# Column and layer of the node `(i, j, k)` of `element` in the arrays of the columns
Base.@propagate_inbounds function get_column_and_layer(horizontal_index_of_element,
                                                       vertical_index_of_element,
                                                       nnodes_,
                                                       i, j, element)
    column = column_index(horizontal_index_of_element, nnodes_, i, j, element)
    layer = vertical_index_of_element[element]
    return column, layer
end

# Unit contravariant direction `nu` of the vertical reference coordinate at the node
# `(i, j, k, element)`
Base.@propagate_inbounds function get_vertical_direction(vertical_direction_and_norm,
                                                         i, j, k, element)
    return SVector(vertical_direction_and_norm[1, i, j, k, element],
                   vertical_direction_and_norm[2, i, j, k, element],
                   vertical_direction_and_norm[3, i, j, k, element])
end

# Norm `n = |Ja^3|` at the node `(i, j, k, element)`, i.e., the entry of the diagonal
# similarity transformation `N` of the vertical operator
Base.@propagate_inbounds function get_vertical_norm(vertical_direction_and_norm,
                                                    i, j, k, element)
    return vertical_direction_and_norm[4, i, j, k, element]
end
end # @muladd
