# Kunz / Wassermann 2D
#
function initial_condition_rainy_mountain_generator(
        parameters::Parameters{RealType};
        mean_flow_velocity = RealType(10),
        relative_humidity  = RealType(0.95)
    ) where RealType

    # Test case parameters (currently fixed)
    N = RealType(1.1e-2)     # Brunt-Väisälä frequency, constant dry static stability
    zm = RealType(5_000)     # height of constant humidity
    T_s = RealType(283.15)   # not T_ref!

    # Physical parameters (taken from TrixiAtmo parameters)
    g    = parameters.earth_gravitational_acceleration
    c_pd = parameters.c_dry_air_const_pressure
    c_vd = parameters.c_dry_air_const_volume
    c_pv = parameters.c_vapour_const_pressure
    c_vv = parameters.c_vapour_const_volume
    p_0  = parameters.ref_pressure
    eps  = parameters.tol_eps
    
    R_d  = c_pd - c_vd
    R_v  = c_pv - c_vv


    function initial_condition(x, t, equations)
        # Exner pressure, solves hydrostatic equation for x[2]
        exner = 1 + g^2 / (c_pd * T_s * N^2) * (exp(-N^2 / g * x[2]) - 1)
        
        # pressure
        p_d = p_0 * exner^(c_pd/ R_d)

        # temperature
        potential_temperature = T_s * exp(N^2 / g * x[2])
        T = potential_temperature * exner

        # rho by ideal gas law
        rho_d = p_d / (R_d * T)

        RH = 0.1 #relative_humidity
        #if (x[2] > zm)
        #    RH = RH * (1.0 + 2 * inv(pi) * atan((zm - x[2]) / 500))
        #end

        e_s = saturation_vapour_pressure(T, equations, equations.microphysics)
        p_v = RH * e_s
        rho_v = p_v / (R_v * T)

        rho_total = rho_d + rho_v

        r_v = 0 #rho_v / rho_total
        r_l = 0  # no cloud water
        r_r = 0  # no rain

        v1 = mean_flow_velocity
        v2 = 0

        p = p_d + p_v  # Dalton's law

        r_v = 0

        return SVector(rho_d, v1, v2, p_d, r_v, r_l, r_r)
    end

    return initial_condition
end
