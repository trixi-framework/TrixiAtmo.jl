# By convention refers to conservative variables
# in the form (rho, rho u, rho v, [rho w, ] rho X) with X being the thermodynamic quantit

@muladd begin
#! format: noindent

@inline function st_hs_relaxation_td(u, k_T, temperature, T_equi,
                                     equations::CompressibleEulerAtmo{3},
                                     ::PotentialTemperature)
    return -k_T * var_td(u, equations) / temperature * (temperature - T_equi)
end

@inline function st_hs_relaxation_td(u, k_T, temperature, T_equi,
                                     equations::CompressibleEulerAtmo{3},
                                     ::EnergyTotal)
    rho_gas = vars_gas(u, equations)
    rho_condens = vars_condens(u, equations)
    rho_cv = rho_cv_total(rho_gas, rho_condens, equations.td_state)
    return -k_T * rho_cv * (temperature - T_equi)
end

function source_terms_hs_relaxation_generator(parameters::Parameters{RealType},
                                              equations::CompressibleEulerAtmo{3}) where {RealType}
    # equations (55)-(58) in the paper
    k_f = 1 / parameters.seconds_per_day         # Damping scale for momentum
    k_a = 1 / (40 * parameters.seconds_per_day)  # Polar relaxation scale
    k_s = 1 / (4 * parameters.seconds_per_day)   # Equatorial relaxation scale
    T_min = RealType(200)                        # Minimum equilibrium temperature
    T_equator = RealType(315)                    # Equatorial equilibrium temperature
    deltaT = RealType(60)                        # Latitudinal temperature difference
    deltaTheta = RealType(10)                    # Vertical temperature difference
    sigma_b = RealType(0.7)                      # Dimensionless damping height

    cp = parameters.c_dry_air_const_pressure
    cv = parameters.c_dry_air_const_volume
    p0 = parameters.ref_pressure

    R = cp - cv

    @inline function source_terms(u, aux, x, t, equations)
        lon, lat, r = cartesian_to_spherical_coordinates(x)
        rho = density_total(u, equations)
        p = pressure(u, equations)
        T = p / (rho * R)
        u_mom = vars_moment(u, equations)

        sigma = p / p0   # "p_0 instead of instantaneous surface pressure"
        delta_sigma = max(0, (sigma - sigma_b) / (1 - sigma_b))   # "height factor"
        k_v = k_f * delta_sigma
        k_T = k_a + (k_s - k_a) * delta_sigma * cos(lat)^4

        T_equi = max(T_min,
                     (T_equator - deltaT * sin(lat)^2 -
                      deltaTheta * log(sigma) * cos(lat)^2) *
                     sigma^(R / cp))

        # project onto r
        dotprod = dot(u_mom, x) / (r * r)

        du_mom = -k_v * (u_mom - dotprod * x)
        du_td = st_hs_relaxation_td(u, k_T, T, T_equi,
                                    equations, equations.td_equation)

        return SVector(du_mom..., du_td)
    end
    return source_terms
end
end # @muladd
