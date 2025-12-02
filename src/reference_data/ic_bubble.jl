@muladd begin
#! format: noindent

# Warm bubble test case from
# - Wicker, L. J., and Skamarock, W. C. (1998)
#   A time-splitting scheme for the elastic equations incorporating
#   second-order Runge–Kutta time differencing
#   [DOI: 10.1175/1520-0493(1998)126%3C1992:ATSSFT%3E2.0.CO;2](https://doi.org/10.1175/1520-0493(1998)126%3C1992:ATSSFT%3E2.0.CO;2)
# See also
# - Bryan and Fritsch (2002)
#   A Benchmark Simulation for Moist Nonhydrostatic Numerical Models
#   [DOI: 10.1175/1520-0493(2002)130<2917:ABSFMN>2.0.CO;2](https://doi.org/10.1175/1520-0493(2002)130<2917:ABSFMN>2.0.CO;2)
# - Carpenter, Droegemeier, Woodward, Hane (1990)
#   Application of the Piecewise Parabolic Method (PPM) to
#   Meteorological Modeling
#   [DOI: 10.1175/1520-0493(1990)118<0586:AOTPPM>2.0.CO;2](https://doi.org/10.1175/1520-0493(1990)118<0586:AOTPPM>2.0.CO;2)

function initial_condition_dry_air_warm_bubble_generator(
    parameters::Parameters{RealType};
    perturbation_center_x     = RealType(10000),
    perturbation_center_z     = RealType(2000),
    perturbation_radius       = RealType(2000),
    potential_temperature_ref = RealType(300),
    pressure_ref              = RealType(10^5),
    velocity_x                = RealType(20),
    velocity_z                = RealType(0)) where RealType

    g = parameters.earth_gravitational_acceleration
    c_p = parameters.c_dry_air_const_pressure
    c_v = parameters.c_dry_air_const_volume

    function initial_condtion(x, t, equations)
        # distance of current x to center of perturbation
        r = sqrt((x[1] - perturbation_center_x)^2 + (x[2] - perturbation_center_z)^2)

        # perturbation in potential temperature
        potential_temperature_perturbation = RealType(0)
        if r <= perturbation_radius
            potential_temperature_perturbation = 2 * cospi(0.5f0 * r / perturbation_radius)^2
        end
        potential_temperature = potential_temperature_ref +
                                potential_temperature_perturbation

        # Exner pressure, solves hydrostatic equation for x[2]
        exner = 1 - g / (c_p * potential_temperature) * x[2]

        # pressure
        R = c_p - c_v    # gas constant (dry air)
        p = pressure_ref * exner^(c_p / R)

        # temperature
        T = potential_temperature * exner

        # density
        rho = p / (R * T)

        return SVector(rho, velocity_x, velocity_z, p)
    end
end
end # @muladd
