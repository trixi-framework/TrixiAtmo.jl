@muladd begin
#! format: noindent

function initial_condition_isothermal_generator(parameters::Parameters{RealType};
                                                temperature_ref = RealType(285),
                                                pressure_ref = RealType(10^5)) where {RealType}
    a = parameters.earth_radius
    g = parameters.earth_gravitational_acceleration
    c_p = parameters.c_dry_air_const_pressure
    c_v = parameters.c_dry_air_const_volume

    R = c_p - c_v

    function initial_condition(x, t, equations)
        r = norm(x)

        # pressure, geopotential formulation
        p = pressure_ref *
            exp(g *
                (a^2 / r - a) /
                (R * temperature_ref))

        # density (via ideal gas law)
        rho = p / (R * temperature_ref)

        return SVector(rho, 0, 0, 0, p)
    end
    return initial_condition
end
end # @muladd
