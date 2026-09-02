@doc raw"""
    initial_condition_adiabatic(x, t,
                                equations::AbstractCompressibleEulerEquations{1, 4})

An atmosphere at rest with constant potential temperature ``\theta = 300\ \mathrm{K}`` in
hydrostatic balance, used to assess the well-balanced property of a discretization.
"""
function initial_condition_adiabatic(x, t,
                                     equations::AbstractCompressibleEulerEquations{1, 4})
    g = equations.g
    c_p = equations.c_p
    c_v = equations.c_v
    p0 = 100_000
    R = c_p - c_v    # gas constant (dry air)
    potential_temperature = 300

    # Exner pressure, solves hydrostatic equation for x[1]
    exner = 1 - g / (c_p * potential_temperature) * x[1]

    # pressure
    p = p0 * exner^(c_p / R)

    # temperature
    T = potential_temperature * exner
    # density
    rho = p / (R * T)
    v1 = 0

    return prim2cons(SVector(rho, v1, p, g * x[1]), equations)
end
