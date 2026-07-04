@muladd begin
#! format: noindent

# Test cases for linearized analytical solution by
# -  Baldauf, Michael and Brdar, Slavko (2013)
#    An analytic solution for linear gravity waves in a channel as a test
#    for numerical models using the non-hydrostatic, compressible {E}uler equations
#    [DOI: 10.1002/qj.2105] (https://doi.org/10.1002/qj.2105)
function initial_condition_gravity_waves_generator(parameters::Parameters{RealType};
                                                   perturbation_center_x = RealType(100000),
                                                   temperature_ref = RealType(250),
                                                   pressure_ref = RealType(10^5),
                                                   velocity_x = RealType(20),
                                                   velocity_z = RealType(0)) where {RealType}
    a = 5_000
    H = 10_000
    DeltaT = 0.001

    g = parameters.earth_gravitational_acceleration
    c_p = parameters.c_dry_air_const_pressure
    c_v = parameters.c_dry_air_const_volume

    R = c_p - c_v

    function initial_condition(x, t, equations)
        delta = g / (R * temperature_ref)
        Tb = DeltaT * sinpi(x[2] / H) * exp(-(x[1] - perturbation_center_x)^2 / a^2)
        rhos = pressure_ref / (temperature_ref * R)
        rho_b = rhos * (-Tb / temperature_ref)
        p = pressure_ref * exp(-delta * x[2])
        rho = rhos * exp(-delta * x[2]) + rho_b * exp(-0.5 * delta * x[2])

        # by convention return primitve variables
        return SVector(rho, velocity_x, velocity_z, p)
    end
    return initial_condition
end
end # @muladd
