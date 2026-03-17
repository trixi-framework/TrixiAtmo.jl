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
    return initial_condition
end


# Moist bubble test case from paper:
# G.H. Bryan, J.M. Fritsch, A Benchmark Simulation for Moist Nonhydrostatic Numerical
# Models, MonthlyWeather Review Vol.130, 2917–2928, 2002,
# https://journals.ametsoc.org/view/journals/mwre/130/12/1520-0493_2002_130_2917_absfmn_2.0.co_2.xml.
function moist_state(y, dz, y0, r_t0, theta_e0, parameters::Parameters{RealType}) where {RealType}
    c_pd = parameters.c_dry_air_const_pressure
    c_vd = parameters.c_dry_air_const_volume
    c_pv = parameters.c_vapour_const_pressure
    c_vv = parameters.c_vapour_const_volume
    c_pl = parameters.c_liquid_water
    p_0  = parameters.ref_pressure
    g = parameters.earth_gravitational_acceleration
    L_00 = parameters.ref_latent_heat_vaporization
    R_d = c_pd - c_vd
    R_v = c_pv - c_vv

    (p, rho, T, r_t, r_v, rho_qv, theta_e) = y
    p0 = y0[1]

    F = zeros(7, 1)
    rho_d = rho / (1 + r_t)
    p_d = R_d * rho_d * T
    T_C = T - convert(RealType, 273.15)
    p_vs = convert(RealType, 611.2) *
           exp(convert(RealType, 17.62) * T_C / (convert(RealType, 243.12) + T_C))
    L = L_00 - (c_pl - c_pv) * T

    F[1] = (p - p0) / dz + g * rho
    F[2] = p - (R_d * rho_d + R_v * rho_qv) * T
    # H = 1 is assumed
    F[3] = (theta_e -
            T * (p_d / p_0)^(-R_d / (c_pd + c_pl * r_t)) *
            exp(L * r_v / ((c_pd + c_pl * r_t) * T)))
    F[4] = r_t - r_t0
    F[5] = rho_qv - rho_d * r_v
    F[6] = theta_e - theta_e0
    a = p_vs / (R_v * T) - rho_qv
    b = rho - rho_qv - rho_d
    # H=1 => phi=0
    F[7] = a + b - sqrt(a * a + b * b)

    return F
end

function generate_function_of_y(dz, y0, r_t0, theta_e0, parameters)
    function function_of_y(y)
        return moist_state(y, dz, y0, r_t0, theta_e0, parameters)
    end
end

struct AtmosphereLayerData{RealType <: Real}
    # structure:  1--> i-layer (z = total_height/precision *(i-1)),  2--> rho, rho_theta, rho_qv, rho_ql
    layer_data::Matrix{RealType} # Contains the layer data for each height
    total_height::RealType # Total height of the atmosphere
    preciseness::Int # Space between each layer data (dz)
end

function AtmosphereLayerData(parameters::Parameters{RealType};
                          moist = true,
                          total_height = 10010, preciseness = 10,
                          rho0 = RealType(1.4),
                          equivalent_potential_temperature = 320,
                          mixing_ratios = (0.2, 0.2)) where {RealType}

    c_pd = parameters.c_dry_air_const_pressure
    c_vd = parameters.c_dry_air_const_volume
    c_pv = parameters.c_vapour_const_pressure
    c_vv = parameters.c_vapour_const_volume
    c_pl = parameters.c_liquid_water
    p0  = parameters.ref_pressure
    R_d = c_pd - c_vd
    R_v = c_pv - c_vv

    if moist
        r_t0, r_v0 = convert.(RealType, mixing_ratios)
    else
        r_t0 = RealType(0)
        r_v0 = RealType(0)
    end
    theta_e0 = equivalent_potential_temperature

    rho_qv0 = rho0 * r_v0
    T0 = theta_e0
    y0 = [p0, rho0, T0, r_t0, r_v0, rho_qv0, theta_e0]

    n = convert(Int, total_height / preciseness)
    dz = convert(RealType, 0.01)
    layer_data = zeros(RealType, n + 1, 4)

    F = generate_function_of_y(dz, y0, r_t0, theta_e0, parameters)
    sol = nlsolve(F, y0)
    p, rho, T, r_t, r_v, rho_qv, theta_e = sol.zero

    rho_d = rho / (1 + r_t)
    rho_ql = rho - rho_d - rho_qv
    kappa_M = (R_d * rho_d + R_v * rho_qv) /
              (c_pd * rho_d + c_pv * rho_qv + c_pl * rho_ql)
    rho_theta = rho * (p0 / p)^kappa_M * T * (1 + (R_v / R_d) * r_v) / (1 + r_t)

    layer_data[1, :] = [rho, rho_theta, rho_qv, rho_ql]
    for i in (1:n)
        y0 = deepcopy(sol.zero)
        dz = preciseness
        F = generate_function_of_y(dz, y0, r_t0, theta_e0, parameters)
        sol = nlsolve(F, y0)
        p, rho, T, r_t, r_v, rho_qv, theta_e = sol.zero

        rho_d = rho / (1 + r_t)
        rho_ql = rho - rho_d - rho_qv
        kappa_M = (R_d * rho_d + R_v * rho_qv) /
                  (c_pd * rho_d + c_pv * rho_qv + c_pl * rho_ql)
        rho_theta = rho * (p0 / p)^kappa_M * T * (1 + (R_v / R_d) * r_v) / (1 + r_t)

        layer_data[i + 1, :] = [rho, rho_theta, rho_qv, rho_ql]
    end

    return AtmosphereLayerData{RealType}(layer_data, total_height, dz)
end

function initial_condition_bryan_fritsch_bubble_generator(
    parameters::Parameters{RealType};
    moist                     = true,
    perturbation_center_x     = RealType(10000),
    perturbation_center_z     = RealType(2000),
    perturbation_radius       = RealType(2000)) where RealType

    c_pd = parameters.c_dry_air_const_pressure
    c_vd = parameters.c_dry_air_const_volume
    c_pv = parameters.c_vapour_const_pressure
    c_vv = parameters.c_vapour_const_volume
    c_pl = parameters.c_liquid_water
    p_0  = parameters.ref_pressure
    L_00 = parameters.ref_latent_heat_vaporization  # RealType(2.5008e6)

    R_d = c_pd - c_vd
    R_v = c_pv - c_vv
    kappa = 1 - c_vd / c_pd
   
    # Create background atmosphere data set
    atmosphere_data = AtmosphereLayerData(parameters; moist = moist)

    function perturb_moist_profile!(x, rho, rho_theta, rho_qv, rho_ql)
        Δθ = 2

        r = sqrt((x[1] - perturbation_center_x)^2 + (x[2] - perturbation_center_z)^2)
        rho_d = rho - rho_qv - rho_ql
        kappa_M = (R_d * rho_d + R_v * rho_qv) /
                (c_pd * rho_d + c_pv * rho_qv + c_pl * rho_ql)
        p_loc = p_0 * (R_d * rho_theta / p_0)^(1 / (1 - kappa_M))
        T_loc = p_loc / (R_d * rho_d + R_v * rho_qv)
        rho_e = (c_vd * rho_d + c_vv * rho_qv + c_pl * rho_ql) * T_loc + L_00 * rho_qv

        # Assume pressure stays constant
        if (r < perturbation_radius && Δθ > 0)
            # Calculate background density potential temperature
            θ_dens = rho_theta / rho * (p_loc / p_0)^(kappa_M - kappa)
            # Calculate perturbed density potential temperature
            θ_dens_new = θ_dens * (1 + Δθ * cospi(0.5f0 * r / perturbation_radius)^2 / 300)
            rt = (rho_qv + rho_ql) / rho_d
            rv = rho_qv / rho_d
            # Calculate moist potential temperature
            θ_loc = θ_dens_new * (1 + rt) / (1 + (R_v / R_d) * rv)
            # Adjust varuables until the temperature is met
            if rt > 0
                while true
                    T_loc = θ_loc * (p_loc / p_0)^kappa
                    T_C = T_loc - convert(RealType, 273.15)
                    # SaturVapor
                    pvs = convert(RealType, 611.2) *
                        exp(convert(RealType, 17.62) * T_C / (convert(RealType, 243.12) + T_C))
                    rho_d_new = (p_loc - pvs) / (R_d * T_loc)
                    rvs = pvs / (R_v * rho_d_new * T_loc)
                    θ_new = θ_dens_new * (1 + rt) / (1 + (R_v / R_d) * rvs)
                    if abs(θ_new - θ_loc) <= θ_loc * 1.0e-12
                        break
                    else
                        θ_loc = θ_new
                    end
                end
            else
                rvs = 0
                T_loc = θ_loc * (p_loc / p_0)^kappa
                rho_d_new = p_loc / (R_d * T_loc)
                θ_new = θ_dens_new * (1 + rt) / (1 + (R_v / R_d) * rvs)
            end
            rho_qv = rvs * rho_d_new
            rho_ql = (rt - rvs) * rho_d_new
            rho = rho_d_new * (1 + rt)
            rho_d = rho - rho_qv - rho_ql
            kappa_M = (R_d * rho_d + R_v * rho_qv) /
                    (c_pd * rho_d + c_pv * rho_qv + c_pl * rho_ql)
            rho_theta = rho * θ_dens_new * (p_loc / p_0)^(kappa - kappa_M)
            rho_e = (c_vd * rho_d + c_vv * rho_qv + c_pl * rho_ql) * T_loc + L_00 * rho_qv
        end
        return SVector(rho, rho_e, rho_qv, rho_ql, T_loc, p_loc)
    end

    function initial_condition(x, t, equations)
        @unpack layer_data, preciseness, total_height = atmosphere_data
        dz = preciseness
        z = x[2]
        if (z > total_height && !(isapprox(z, total_height)))
            error("The atmosphere does not match the simulation domain")
        end
        n = convert(Int, floor((z + eps()) / dz)) + 1
        z_l = (n - 1) * dz
        (rho_l, rho_theta_l, rho_qv_l, rho_ql_l) = layer_data[n, :]
        z_r = n * dz
        if (z_l == total_height)
            z_r = z_l + dz
            n = n - 1
        end
        (rho_r, rho_theta_r, rho_qv_r, rho_ql_r) = layer_data[n + 1, :]
        rho = (rho_r * (z - z_l) + rho_l * (z_r - z)) / dz
        rho_theta = rho *
                    (rho_theta_r / rho_r * (z - z_l) + rho_theta_l / rho_l * (z_r - z)) /
                    dz
        rho_qv = rho * (rho_qv_r / rho_r * (z - z_l) + rho_qv_l / rho_l * (z_r - z)) / dz
        rho_ql = rho * (rho_ql_r / rho_r * (z - z_l) + rho_ql_l / rho_l * (z_r - z)) / dz

        rho, rho_e, rho_qv, rho_ql, T_loc, p_loc = perturb_moist_profile!(x, rho, rho_theta,
                                                                rho_qv, rho_ql)
        v1 = 0
        v2 = 0
        return SVector(rho, v1, v2, p_loc, rho_qv / rho, rho_ql / rho)
    end
    return initial_condition
end

end # @muladd
