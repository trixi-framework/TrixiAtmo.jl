# Implemented by Fabian Höck following
# Sabine Doppler, Philip L. Lederer, Joachim Schöberl, Henry von Wahl
# A discontinuous Galerkin approach for atmospheric flows with implicit condensation
# Journal of Computational Physics, Volume 499, 2024
# https://doi.org/10.1016/j.jcp.2023.112713

using Trixi
using NLsolve: nlsolve
import Trixi: varnames, cons2prim, cons2entropy,
              flux, max_abs_speeds, max_abs_speed_naive,
              boundary_condition_slip_wall, LaplaceDiffusion2D

@muladd begin
#! format: noindent

###  equation, parameters and constants  ###

@doc raw"""
    CompressibleRainyEulerEquations2D{RealT <: Real} <:
        AbstractCompressibleRainyEulerEquations{2, 9}

The compressible Euler equations in two dimensions with gravity and separate densities for
dry air ($\rho_d$), moist air ($\rho_m$), and rain ($\rho_r).
```math
\begin{alignat}{4}
    \partial_t \, &&(\rho_d) &+ \nabla \cdot (\rho_d\, v) &&= & 0\,,\\
    \partial_t \, &&(\rho_m) &+ \nabla \cdot (\rho_m\, v) &&= & -S_r \,,\\
    \partial_t \, &&(\rho_r) &+ \nabla \cdot (\rho_r\, (v-v_r\, e_z)) &&= & S_r\,,\\
    \partial_t\,&&(\rho v) &+ \nabla \cdot (\rho v \otimes v - \rho_r\, v_r \,v\otimes e_z + p\cdot\textrm{Id}) &&= & -\rho\, g\, e_z\,,\\
    \partial_t\, &&(E) &+ \nabla \cdot ((E + p)\,v - (c_l\, T + 1/2\, \langle v, v\rangle)\, \rho_r\, v_r\, e_z) &&= &\,\,-\rho\, g\, \langle e_z, v\rangle\,.
\end{alignat}
```
Here, moist air is the sum of vapor ($\rho_v$) and cloud water ($\rho_c$), which are given implicitly by the nonlinear system
```math
\begin{align}
    (c_{vd}\, \rho_d + c_{vv}\, \rho_v + c_l\, \rho_l)\,T + \rho_v\, L_{\textrm{ref}} + 1/2\,\rho\,\langle v, v\rangle - E &= 0\,,\\
    \min \left( \frac{e_s(T)}{R_v\,T}, \rho_m\right) - \rho_v &=0\,,\\
    \rho_m - \rho_v - \rho_c &=0\,,
\end{align}
```
where $e_s(T)$ is the modeled saturation vapor pressure, $v_r$ the modeled terminal rain,
fall velocity, and $S_r$ a modeled source term for rain conversion.

## Reference
S. Doppler et al. (2024). A discontinuous Galerkin approach for atmospheric flows with
implicit condensation. Journal of Computational Physics. 499:112713.
[DOI: 10.1016/j.jcp.2023.112713](https://doi.org/10.1016/j.jcp.2023.112713)
"""
struct CompressibleRainyEulerEquations2D{RealT <: Real} <:
       AbstractCompressibleRainyEulerEquations{2, 9}
    # Specific heat capacities:
    c_liquid_water           :: RealT
    c_dry_air_const_pressure :: RealT
    c_dry_air_const_volume   :: RealT
    c_vapour_const_pressure  :: RealT
    c_vapour_const_volume    :: RealT

    # Gas constants:
    R_dry_air :: RealT
    R_vapour  :: RealT
    eps       :: RealT

    # Reference values:
    ref_saturation_pressure  :: RealT
    ref_temperature          :: RealT
    ref_latent_heat_vap_temp :: RealT
    ref_pressure             :: RealT

    # Other:
    gravity          :: RealT
    rain_water_distr :: RealT
    v_mean_rain      :: RealT
end

function CompressibleRainyEulerEquations2D(; RealT = Float64)
    # Specific heat capacities:
    c_liquid_water = 4186
    c_dry_air_const_pressure = 1004
    c_dry_air_const_volume = 717
    c_vapour_const_pressure = 1885
    c_vapour_const_volume = 1424

    # Gas constants:
    R_dry_air = c_dry_air_const_pressure - c_dry_air_const_volume
    R_vapour = c_vapour_const_pressure - c_vapour_const_volume
    eps = R_dry_air / R_vapour

    # Reference values:
    ref_saturation_pressure = convert(RealT, 610.7)    # This needs to be adjusted if ref_temperature is changed!
    ref_temperature = convert(RealT, 273.15)
    ref_latent_heat_vap_temp = 2500000 #3147620.0
    ref_pressure = 100000

    # Other:
    gravity = EARTH_GRAVITATIONAL_ACCELERATION
    rain_water_distr = 8000000
    v_mean_rain = 130

    return CompressibleRainyEulerEquations2D{RealT}(c_liquid_water,
                                                    c_dry_air_const_pressure,
                                                    c_dry_air_const_volume,
                                                    c_vapour_const_pressure,
                                                    c_vapour_const_volume, R_dry_air,
                                                    R_vapour, eps,
                                                    ref_saturation_pressure,
                                                    ref_temperature,
                                                    ref_latent_heat_vap_temp,
                                                    ref_pressure, gravity,
                                                    rain_water_distr, v_mean_rain)
end

###  conversion  ###

@inline function cons2prim(u, equations::CompressibleRainyEulerEquations2D)
    # densities
    rho_dry, rho_moist, rho_rain, rho, rho_inv = densities(u, equations)

    # velocity
    v1, v2 = velocities(u, rho_inv, equations)

    # energy density
    energy = energy_density(u, equations)

    # nonlinear system
    rho_vapour, rho_cloud, temperature = cons2nonlinearsystemsol(u, equations)

    return SVector(rho_dry, rho_moist, rho_rain, v1, v2, energy, rho_vapour, rho_cloud,
                   temperature)
end

# converts consverved to entropy variables
@inline function cons2entropy(u, equations::CompressibleRainyEulerEquations2D)
    # constants
    c_vd = equations.c_dry_air_const_volume
    c_vv = equations.c_vapour_const_volume
    c_l = equations.c_liquid_water
    R_d = equations.R_dry_air
    R_v = equations.R_vapour
    L_ref = equations.ref_latent_heat_vap_temp

    # densities
    rho_dry, rho_moist, rho_rain, rho, rho_inv = densities(u, equations)

    # nonlinear system
    rho_vapour, rho_cloud, temperature = cons2nonlinearsystemsol(u, equations)
    ln_temperature = log(temperature)
    inv_temperature = inv(temperature)

    # velocity
    v1, v2 = velocities(u, rho_inv, equations)
    v_squared_temp = 0.5f0 * (v1^2 + v2^2) * inv_temperature

    # check for zero density
    rho_vapour_log = 0

    if (rho_vapour > 0.0)
        rho_vapour_log = log(rho_vapour)
    end

    omega_dry = c_vd * ln_temperature - R_d * log(rho_dry) + v_squared_temp - c_vd - R_d
    omega_vapour = c_vv * ln_temperature - R_v * rho_vapour_log + v_squared_temp -
                   c_vv - R_v - L_ref * inv_temperature
    omega_liquid = c_l * ln_temperature + v_squared_temp - c_l
    omega_momentum_1 = -v1 * inv_temperature
    omega_momentum_2 = -v2 * inv_temperature
    omega_energy = inv_temperature

    return SVector(omega_dry, omega_vapour + omega_liquid, omega_liquid,
                   omega_momentum_1, omega_momentum_2, omega_energy, 0, 0, 0)
end

# adapted from compressible_moist_euler_2d.jl
@inline function cons2eq_pot_temp(u, equations::CompressibleRainyEulerEquations2D{RealT}) where {RealT}
    # constants
    c_l = equations.c_liquid_water
    c_pd = equations.c_dry_air_const_pressure
    c_pv = equations.c_vapour_const_pressure
    R_d = equations.R_dry_air
    R_v = equations.R_vapour
    ref_p = equations.ref_pressure
    ref_temp = equations.ref_temperature
    ref_L = equations.ref_latent_heat_vap_temp

    # densities
    rho_dry, rho_moist, rho_rain, rho, rho_inv = densities(u, equations)

    # velocity
    v1, v2 = velocities(u, rho_inv, equations)

    # energy density
    energy = energy_density(u, equations)

    # nonlinear system
    rho_vapour, rho_cloud, temperature = cons2nonlinearsystemsol(u, equations)

    # pressure
    p = pressure(u, equations)

    p_v = rho_vapour * R_v * temperature
    p_d = p - p_v
    T_C = temperature - ref_temp
    p_vs = convert(RealT, 611.2) * exp(convert(RealT, 17.62) * T_C / (convert(RealT, 243.12) + T_C))
    H = p_v / p_vs
    r_v = rho_vapour / rho_dry
    r_c = rho_cloud / rho_dry
    r_r = rho_rain / rho_dry
    L_v = ref_L + (c_pv - c_l) * temperature
    c_p = c_pd + (r_v + r_c + r_r) * c_l

    # equivalent potential temperature
    eq_pot = (temperature * (ref_p / p_d)^(R_d / c_p) * H^(-r_v * R_v / c_p) *
              exp(L_v * r_v * inv(c_p * temperature)))

    return SVector(rho_dry, rho_vapour, rho_cloud, rho_rain, v1, v2, eq_pot, p)
end

@inline function cons2nonlinearsystemsol(u,
                                         equations::CompressibleRainyEulerEquations2D)
    # constants
    c_vd = equations.c_dry_air_const_volume
    #ref_temp = equations.ref_temperature

    # densities
    rho_dry, rho_moist, rho_rain, rho, rho_inv = densities(u, equations)

    # velocity
    v1, v2 = velocities(u, rho_inv, equations)

    # energy density
    energy = energy_density(u, equations)

    # recover temperature explicitly from energy when other variables are zero
    if (rho_moist == 0.0 && rho_rain == 0.0)
        # energy density definition without ref_temp for dry case
        energy_kinetic = 0.5f0 * (v1^2 + v2^2) * rho
        temperature = (energy - energy_kinetic) / (c_vd * rho) #+ ref_temp

        if (temperature < 0.0)
            error("temp negative")
        end

        return SVector(0, 0, temperature)
    else
        # experimental and overly simple positivity check
        rho_vapour = u[7]
        rho_cloud = u[8]
        temperature = u[9]

        if (rho_vapour < 0.0 && isapprox(rho_vapour, 0.0, atol = 1e-15))
            rho_vapour = 0
        end

        if (rho_cloud < 0.0 && isapprox(rho_cloud, 0.0, atol = 1e-15))
            rho_cloud = 0
        end

        return SVector(rho_vapour, rho_cloud, temperature)
    end
end

# for convenience TODO rename
@inline function cons2speeds(u, equations::CompressibleRainyEulerEquations2D)
    # densities
    rho_dry, rho_moist, rho_rain, rho, rho_inv = densities(u, equations)

    # velocity
    v1, v2 = velocities(u, rho_inv, equations)

    # get speed of sound 
    v_sound = speed_of_sound(u, equations)[1]

    # get terminal velocity rain
    v_r = terminal_velocity_rain(rho_moist, rho_rain, equations)

    return SVector(v1, v2, v_sound, v_r)
end

###  varnames  ###

varnames(::typeof(cons2cons), ::CompressibleRainyEulerEquations2D) = ("rho_dry",
                                                                      "rho_moist",
                                                                      "rho_rain",
                                                                      "rho_v1",
                                                                      "rho_v2",
                                                                      "energy_density",
                                                                      "rho_vapour",
                                                                      "rho_cloud",
                                                                      "temperature")

varnames(::typeof(cons2prim), ::CompressibleRainyEulerEquations2D) = ("rho_dry",
                                                                      "rho_moist",
                                                                      "rho_rain",
                                                                      "v1", "v2",
                                                                      "energy_density",
                                                                      "rho_vapour",
                                                                      "rho_cloud",
                                                                      "temperature")

varnames(::typeof(cons2eq_pot_temp), ::CompressibleRainyEulerEquations2D) = ("rho_dry",
                                                                             "rho_vapour",
                                                                             "rho_cloud",
                                                                             "rho_rain",
                                                                             "v1", "v2",
                                                                             "eq_pot_temp",
                                                                             "pressure")

###  physics variables  ###

@inline function densities(u, equations::CompressibleRainyEulerEquations2D)
    rho_dry = u[1]
    rho_moist = u[2]
    rho_rain = u[3]
    rho = rho_dry + rho_moist + rho_rain
    rho_inv = inv(rho)

    return SVector(rho_dry, rho_moist, rho_rain, rho, rho_inv)
end

@inline function velocities(u, rho_inv, equations::CompressibleRainyEulerEquations2D)
    return SVector(u[4] * rho_inv, u[5] * rho_inv)
end

@inline function energy_density(u, equations::CompressibleRainyEulerEquations2D)
    return u[6]
end

@inline function pressure(u, equations::CompressibleRainyEulerEquations2D)
    # constants
    R_d = equations.R_dry_air
    R_v = equations.R_vapour

    # densities
    rho_dry, rho_moist, rho_rain, rho, rho_inv = densities(u, equations)

    rho_vapour, _, temperature = cons2nonlinearsystemsol(u, equations)

    p = (R_d * rho_dry + R_v * rho_vapour) * temperature

    return p
end

@inline function speed_of_sound(u, equations::CompressibleRainyEulerEquations2D)
    # constants
    c_l = equations.c_liquid_water
    c_vd = equations.c_dry_air_const_volume
    c_vv = equations.c_vapour_const_volume
    R_d = equations.R_dry_air
    R_v = equations.R_vapour

    # densities
    rho_dry, _, rho_rain, _, rho_inv = densities(u, equations)

    # recover rho_vapour, rho_cloud, temperature from nonlinear system
    rho_vapour, rho_cloud, temperature = cons2nonlinearsystemsol(u, equations)
    if (rho_vapour < 0.0)
        error("rho vapour less than zero")
    end
    if (rho_cloud < 0.0)
        error("rho cloud less than zero")
    end

    # formula
    p = pressure(u, equations)
    q_v = rho_vapour / rho_dry
    q_l = (rho_cloud + rho_rain) / rho_dry
    gamma_m = 1 + (R_d + R_v * q_v) / (c_vd + c_vv * q_v + c_l * q_l)

    if (rho_inv < 0.0)
        error("rho less than zero")
    elseif (p < 0.0)
        error("pressure less than zero")
    end

    v_sound = sqrt(gamma_m * p * rho_inv)

    return SVector(v_sound, gamma_m)
end

@inline function terminal_velocity_rain(rho_moist, rho_rain,
                                        equations::CompressibleRainyEulerEquations2D{RealT}) where {RealT}
    # constants
    N_0 = equations.rain_water_distr
    v_0 = equations.v_mean_rain

    # formula ( \Gamma(4.5) / 6 ~= 1.9386213994279082 )
    if (rho_rain > 0.0)
        v_terminal_rain = v_0 * convert(RealT, 1.9386213994279082) *
                          (rho_rain / (pi * (rho_moist + rho_rain) * N_0))^(0.125f0)
    else
        v_terminal_rain = 0
    end

    return v_terminal_rain
end

@inline function saturation_vapour_pressure(temperature,
                                            equations::CompressibleRainyEulerEquations2D)
    # constants
    c_l = equations.c_liquid_water
    c_pv = equations.c_vapour_const_pressure
    R_v = equations.R_vapour
    ref_s_p = equations.ref_saturation_pressure
    ref_temp = equations.ref_temperature
    ref_L = equations.ref_latent_heat_vap_temp

    # testing 
    if (temperature < 0.0)
        display(temperature)
        error("temp less than zero")
    end

    # Clausius Clapeyron formula
    p_vapour_saturation = ref_s_p * (temperature / ref_temp)^((c_pv - c_l) / R_v)
    p_vapour_saturation *= exp(((ref_L - (c_pv - c_l) * ref_temp) / R_v) *
                               (1 / ref_temp - 1 / temperature))

    return p_vapour_saturation
end

@inline function saturation_vapour_pressure_derivative(temperature,
                                                       equations::CompressibleRainyEulerEquations2D)
    # constants
    c_l = equations.c_liquid_water
    c_pv = equations.c_vapour_const_pressure
    R_v = equations.R_vapour
    ref_s_p = equations.ref_saturation_pressure
    ref_temp = equations.ref_temperature
    ref_L = equations.ref_latent_heat_vap_temp

    # testing 
    if (temperature < 0.0)
        display(temperature)
        error("temp less than zero")
    end

    const_1 = (c_pv - c_l) / R_v
    const_2 = (ref_L - (c_pv - c_l) * ref_temp) / R_v

    p_vapour_saturation_derivative = ref_s_p / (ref_temp^const_1)
    p_vapour_saturation_derivative *= (const_1 * temperature^(const_1 - 1) +
                                       const_2 * temperature^(const_1 - 2))
    p_vapour_saturation_derivative *= exp(const_2 * (1 / ref_temp - 1 / temperature))

    return p_vapour_saturation_derivative
end

###  pde discretization  ###

@inline function flux(u, orientation::Integer,
                      equations::CompressibleRainyEulerEquations2D)
    # constants
    c_l = equations.c_liquid_water

    #densities
    rho_dry, rho_moist, rho_rain, rho, rho_inv = densities(u, equations)

    # recover rho_vapour, rho_cloud, temperature from nonlinear system
    _, _, temperature = cons2nonlinearsystemsol(u, equations)

    # velocity
    v1, v2 = velocities(u, rho_inv, equations)
    v_r = terminal_velocity_rain(rho_moist, rho_rain, equations)

    # pressure
    p = pressure(u, equations)

    # energy density
    energy = energy_density(u, equations)

    # flux for orientation cases 
    if (orientation == 1)
        # mass
        f1 = rho_dry * v1
        f2 = rho_moist * v1
        f3 = rho_rain * v1

        # momentum
        f4 = rho * v1 * v1 + p
        f5 = rho * v1 * v2

        # energy
        f6 = (energy + p) * v1

    else
        # mass
        f1 = rho_dry * v2
        f2 = rho_moist * v2
        f3 = rho_rain * (v2 - v_r)

        # momentum
        f4 = rho * v1 * v2 - rho_rain * v_r * v1
        f5 = rho * v2 * v2 + p - rho_rain * v_r * v2

        # energy
        f6 = (energy + p) * v2 -
             (c_l * temperature + 0.5f0 * (v1^2 + v2^2)) * rho_rain * v_r
    end

    return SVector(f1, f2, f3, f4, f5, f6,
                   0, 0, 0)
end

@inline function flux(u, normal_direction::AbstractVector,
                      equations::CompressibleRainyEulerEquations2D)
    # constants
    c_l = equations.c_liquid_water

    # densities
    rho_dry, rho_moist, rho_rain, rho, rho_inv = densities(u, equations)

    # recover rho_vapour, rho_cloud, temperature from nonlinear system
    _, _, temperature = cons2nonlinearsystemsol(u, equations)

    # velocity
    v1, v2 = velocities(u, rho_inv, equations)
    v_r = terminal_velocity_rain(rho_moist, rho_rain, equations)

    # normal velocities
    v_normal = v1 * normal_direction[1] + v2 * normal_direction[2]
    v_r_normal = v_r * normal_direction[2]

    # pressure
    p = pressure(u, equations)

    # energy density
    energy = energy_density(u, equations)

    # flux
    # mass
    f1 = rho_dry * v_normal
    f2 = rho_moist * v_normal
    f3 = rho_rain * (v_normal - v_r_normal)

    # momentum
    f4 = rho * v_normal * v1 + p * normal_direction[1] - rho_rain * v_r_normal * v1
    f5 = rho * v_normal * v2 + p * normal_direction[2] - rho_rain * v_r_normal * v2

    # energy
    f6 = (energy + p) * v_normal -
         (c_l * temperature + 0.5f0 * (v1^2 + v2^2)) * rho_rain * v_r_normal

    return SVector(f1, f2, f3, f4, f5, f6,
                   0, 0, 0)
end

# no Coriolis term
@inline function source_terms_rainy(u, x, t,
                                    equations::CompressibleRainyEulerEquations2D{RealT}) where {RealT}
    # constants
    R_v = equations.R_vapour
    ref_temp = equations.ref_temperature
    g = equations.gravity

    # name needed variables
    rho_v2 = u[5]

    # densities
    rho_dry, rho_moist, rho_rain, rho, rho_inv = densities(u, equations)

    # recover rho_vapour, rho_cloud, temperature from nonlinear system
    rho_vapour, rho_cloud, temperature = cons2nonlinearsystemsol(u, equations)

    rho_vs = saturation_vapour_pressure(temperature, equations) / (R_v * temperature)

    # source terms phase change
    S_evaporation = (convert(RealT, 3.86e-3) - convert(RealT, 9.41e-5) * (temperature - ref_temp)) *
                    (1 + convert(RealT, 9.1) * rho_rain^(0.1875f0))
    S_evaporation *= (rho_vs - rho_vapour) * rho_rain^(0.5f0)
    S_auto_conversion = convert(RealT, 0.001) * rho_cloud
    S_accretion = convert(RealT, 1.72) * rho_cloud * rho_rain^(0.875f0)
    S_rain = S_auto_conversion + S_accretion - S_evaporation

    return SVector(0, -S_rain, S_rain, 0,
                   -rho * g, -rho_v2 * g, 0, 0, 0)
end

# no phase changes and no Coriolis term
@inline function source_terms_no_phase_change(u, x, t,
                                              equations::CompressibleRainyEulerEquations2D)
    # constants
    g = equations.gravity

    # name needed variables
    rho_v2 = u[5]

    # densities
    rho_dry, rho_moist, rho_rain, rho, rho_inv = densities(u, equations)

    return SVector(0, 0, 0, 0,
                   -g * rho, -g * rho_v2, 0, 0, 0)
end

@inline function max_abs_speeds(u, equations::CompressibleRainyEulerEquations2D)
    # name needed variables
    v1, v2, v_sound, v_r = cons2speeds(u, equations)

    return SVector((abs(v1) + v_sound), (abs(v2) + v_sound + abs(v_r)))
end

@inline function max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
                                     equations::CompressibleRainyEulerEquations2D)
    # name needed variables
    v1_ll, v2_ll, v_sound_ll, v_r_ll = cons2speeds(u_ll, equations)
    v1_rr, v2_rr, v_sound_rr, v_r_rr = cons2speeds(u_rr, equations)

    # calculate upper bounds for left and right speed
    v_ll_max = abs(v1_ll * normal_direction[1] + v2_ll * normal_direction[2])
    v_ll_max += abs(v_r_ll * normal_direction[2])

    v_rr_max = abs(v1_rr * normal_direction[1] + v2_rr * normal_direction[2])
    v_rr_max += abs(v_r_rr * normal_direction[2])

    return max(v_ll_max, v_rr_max) +
           max(v_sound_ll, v_sound_rr) * norm(normal_direction)
end

@inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
                                     equations::CompressibleRainyEulerEquations2D)
    # name needed variables
    v1_ll, v2_ll, v_sound_ll, v_r_ll = cons2speeds(u_ll, equations)
    v1_rr, v2_rr, v_sound_rr, v_r_rr = cons2speeds(u_rr, equations)

    if (orientation == 1)
        v_ll = abs(v1_ll)
        v_rr = abs(v1_rr)
    else
        v_ll = abs(v2_ll)
        v_rr = abs(v2_rr)
    end
    # experimental
    return max(v_ll, v_rr) + max(v_sound_ll, v_sound_rr) + max(abs(v_r_ll), abs(v_r_rr))
end

###  boundary conditions  ###

# adapted from compressible_moist_euler_2d.jl which was probably adapted from compressible_euler_2d.jl
@inline function boundary_condition_slip_wall(u_inner, normal_direction::AbstractVector,
                                              x, t,
                                              surface_flux_function,
                                              equations::CompressibleRainyEulerEquations2D)
    norm_ = norm(normal_direction)
    # Normalize the vector without using `normalize` since we need to multiply by the `norm_` later
    normal = normal_direction / norm_

    # rotate the internal solution state
    u_local = rotate_to_x(u_inner, normal, equations)

    # name needed variables
    rho_v1 = u_local[4]

    # densities
    rho_dry_local, rho_moist_local, rho_rain_local, rho_local, rho_inv_local = densities(u_local,
                                                                                         equations)

    # velocities
    v_normal = rho_v1 * rho_inv_local
    v_sound, gamma = speed_of_sound(u_local, equations)

    # pressure
    p_local = pressure(u_local, equations)

    # Get the solution of the pressure Riemann problem
    # See Section 6.3.3 of
    # Eleuterio F. Toro (2009)
    # Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction
    # [DOI: 10.1007/b79761](https://doi.org/10.1007/b79761)
    if v_normal <= 0.0
        p_star = p_local *
                 (1 + 0.5f0 * (gamma - 1) * v_normal / v_sound)^(2 * gamma *
                                                                 inv(gamma - 1))
    else # v_normal > 0.0
        A = 2 / ((gamma + 1) * rho_local)
        B = p_local * (gamma - 1) / (gamma + 1)
        p_star = p_local +
                 0.5f0 * v_normal / A *
                 (v_normal + sqrt(v_normal^2 + 4 * A * (p_local + B)))
    end

    # For the slip wall we directly set the flux as the normal velocity is zero
    return SVector(0, 0, 0,
                   p_star * normal[1] * norm_,
                   p_star * normal[2] * norm_,
                   0, 0, 0, 0)
end

# same as in compressible_euler_2d.jl
@inline function boundary_condition_slip_wall(u_inner, normal_direction::AbstractVector,
                                              direction, x, t,
                                              surface_flux_function,
                                              equations::CompressibleRainyEulerEquations2D)
    # flip sign of normal to make it outward pointing, then flip the sign of the normal flux back
    # to be inward pointing on the -x and -y sides due to the orientation convention used by StructuredMesh
    if isodd(direction)
        boundary_flux = -boundary_condition_slip_wall(u_inner, -normal_direction,
                                                      x, t, surface_flux_function,
                                                      equations)
    else
        boundary_flux = boundary_condition_slip_wall(u_inner, normal_direction,
                                                     x, t, surface_flux_function,
                                                     equations)
    end

    return boundary_flux
end

# same as in compressible_euler_2d.jl
@inline function rotate_to_x(u, normal_vector,
                             equations::CompressibleRainyEulerEquations2D)
    # cos and sin of the angle between the x-axis and the normalized normal_vector are
    # the normalized vector's x and y coordinates respectively (see unit circle).
    c = normal_vector[1]
    s = normal_vector[2]

    return SVector(u[1], u[2], u[3],
                   c * u[4] + s * u[5],
                   -s * u[4] + c * u[5],
                   u[6], u[7], u[8], u[9])
end

# for parabolic terms (LaplaceDiffusion2D)
@inline function boundary_condition_laplace(flux_inner, u_inner, normal::AbstractVector,
                                            x, t, operator_type::Trixi.Gradient,
                                            equations_parabolic::LaplaceDiffusion2D)
    return u_inner
end

@inline function boundary_condition_laplace(flux_inner, u_inner, normal::AbstractVector,
                                            x, t, operator_type::Trixi.Divergence,
                                            equations_parabolic::LaplaceDiffusion2D)
    return flux_inner
end

###  Nonlinear System Residual  ###

# in preparation for a callback to solve the nonlinear system
@inline function saturation_residual(u, guess,
                                     equations::CompressibleRainyEulerEquations2D)
    # constants
    c_l = equations.c_liquid_water
    c_vd = equations.c_dry_air_const_volume
    c_vv = equations.c_vapour_const_volume
    R_v = equations.R_vapour
    L_ref = equations.ref_latent_heat_vap_temp

    # densities
    rho_dry, rho_moist, rho_rain, rho, rho_inv = densities(u, equations)

    # velocity
    v1, v2 = velocities(u, rho_inv, equations)

    # energy density
    energy = energy_density(u, equations)

    # define residual
    residual1 = (c_vd * rho_dry + c_vv * guess[1] + c_l * (guess[2] + rho_rain)) *
                guess[3]
    residual1 += guess[1] * L_ref
    residual1 -= (energy - rho * 0.5f0 * (v1^2 + v2^2))

    residual2 = min(saturation_vapour_pressure(guess[3], equations) / (R_v * guess[3]),
                    rho_moist)
    residual2 -= guess[1]
    residual2 *= 1e7

    residual3 = rho_moist
    residual3 -= guess[1] + guess[2]
    residual3 *= 1e7

    return SVector(residual1, residual2, residual3)
end

@inline function saturation_residual_jacobian(u, guess,
                                              equations::CompressibleRainyEulerEquations2D)
    # constants
    c_l = equations.c_liquid_water
    c_vd = equations.c_dry_air_const_volume
    c_vv = equations.c_vapour_const_volume
    R_v = equations.R_vapour
    L_ref = equations.ref_latent_heat_vap_temp

    # densities
    rho_dry, rho_moist, rho_rain, rho, rho_inv = densities(u, equations)

    # saturation
    svp = saturation_vapour_pressure(guess[3], equations)
    svp_t = saturation_vapour_pressure_derivative(guess[3], equations)

    # define jacobian
    J_11 = c_vv * guess[3] + L_ref
    J_12 = c_l * guess[3]
    J_13 = c_vd * rho_dry + c_vv * guess[1] + c_l * (guess[2] + rho_rain)

    J_21 = -10000000
    J_22 = 0

    if (svp / (R_v * guess[3]) < rho_moist)
        J_23 = (svp_t * guess[3] - svp) / (R_v * guess[3]^2) * 1e7
    else
        J_23 = 0
    end

    J_31 = -10000000
    J_32 = -10000000
    J_33 = 0

    return SMatrix{3, 3}(J_11, J_21, J_31, J_12, J_22, J_32, J_13, J_23, J_33)
end

# Low Mach number approximate Riemann solver (LMARS) from
# X. Chen, N. Andronova, B. Van Leer, J. E. Penner, J. P. Boyd, C. Jablonowski, S.
# Lin, A Control-Volume Model of the Compressible Euler Equations with a Vertical Lagrangian
# Coordinate Monthly Weather Review Vol. 141.7, pages 2526–2544, 2013,
# https://journals.ametsoc.org/view/journals/mwre/141/7/mwr-d-12-00129.1.xml.
# adapted from compressible_moist_euler_2d.jl, does NOT work with rain!
@inline function flux_LMARS(u_ll, u_rr, normal_direction::AbstractVector,
                            equations::CompressibleRainyEulerEquations2D)
    # constants
    a = 360

    # densities
    rho_dry_ll, rho_moist_ll, rho_rain_ll, rho_ll, rho_inv_ll = densities(u_ll,
                                                                          equations)
    rho_dry_rr, rho_moist_rr, rho_rain_rr, rho_rr, rho_inv_rr = densities(u_rr,
                                                                          equations)

    # pressure
    p_ll = pressure(u_ll, equations)
    p_rr = pressure(u_rr, equations)

    # velocities
    v1_ll, v2_ll = velocities(u_ll, rho_inv_ll, equations)
    v1_rr, v2_rr = velocities(u_rr, rho_inv_rr, equations)

    v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]
    norm_ = norm(normal_direction)

    # diffusion parameter 0.0 < beta <= 1.0
    beta = 1

    # interface flux components
    rho = 0.5f0 * (rho_ll + rho_rr)
    p_interface = 0.5f0 * (p_ll + p_rr) - beta * 0.5f0 * a * rho * (v_rr - v_ll) / norm_
    v_interface = 0.5f0 * (v_ll + v_rr) - beta * 1 / (2 * a * rho) * (p_rr - p_ll) * norm_

    if (v_interface > 0)
        f1, f2, _, f4, f5, f6, _, _, _ = u_ll * v_interface
        f6 += p_ll * v_interface
    else
        f1, f2, _, f4, f5, f6, _, _, _ = u_rr * v_interface
        f6 += p_rr * v_interface
    end

    return SVector(f1, f2, 0,
                   f4 + p_interface * normal_direction[1],
                   f5 + p_interface * normal_direction[2],
                   f6, 0, 0, 0)
end

"""
    flux_ec_rain(u_ll, u_rr, orientation_or_normal_direction,
                 equations::CompressibleRainyEulerEquations2D)

Entropy-conserving flux including rain, derived in:
Fabian Höck
A Discontinuous Galerkin Method for Moist Atmospheric Dynamics with Rain
Master's thesis, University of Cologne, 2025
"""
@inline function flux_ec_rain(u_ll, u_rr, normal_direction::AbstractVector,
                              equations::CompressibleRainyEulerEquations2D)
    # constants
    c_l = equations.c_liquid_water
    c_vd = equations.c_dry_air_const_volume
    c_vv = equations.c_vapour_const_volume
    R_d = equations.R_dry_air
    R_v = equations.R_vapour
    L_ref = equations.ref_latent_heat_vap_temp

    # densities and temperatures
    rho_dry_ll, rho_moist_ll, rho_rain_ll, rho_ll, rho_inv_ll = densities(u_ll,
                                                                          equations)
    rho_dry_rr, rho_moist_rr, rho_rain_rr, rho_rr, rho_inv_rr = densities(u_rr,
                                                                          equations)
    rho_vapour_ll, rho_cloud_ll, temperature_ll = cons2nonlinearsystemsol(u_ll,
                                                                          equations)
    rho_vapour_rr, rho_cloud_rr, temperature_rr = cons2nonlinearsystemsol(u_rr,
                                                                          equations)
    inv_temperature_ll = inv(temperature_ll)
    inv_temperature_rr = inv(temperature_rr)

    # velocities
    v1_ll, v2_ll = velocities(u_ll, rho_inv_ll, equations)
    v1_rr, v2_rr = velocities(u_rr, rho_inv_rr, equations)
    vr_ll = terminal_velocity_rain(rho_vapour_ll + rho_cloud_ll, rho_rain_ll, equations)
    vr_rr = terminal_velocity_rain(rho_vapour_rr + rho_cloud_rr, rho_rain_rr, equations)

    # velocity averages
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v1_square_avg = 0.5f0 * (v1_ll^2 + v1_rr^2)
    v2_square_avg = 0.5f0 * (v2_ll^2 + v2_rr^2)
    K_avg = 0.5f0 * (v1_square_avg + v2_square_avg)
    v_dot_n_avg = normal_direction[1] * v1_avg + normal_direction[2] * v2_avg
    vr_normal_avg = 0.5f0 * normal_direction[2] * (vr_ll + vr_rr)

    # density averages
    rho_dry_avg = 0.5f0 * (rho_dry_ll + rho_dry_rr)
    rho_vapour_avg = 0.5f0 * (rho_vapour_ll + rho_vapour_rr)
    rho_cloud_avg = 0.5f0 * (rho_cloud_ll + rho_cloud_rr)
    rho_rain_avg = 0.5f0 * (rho_rain_ll + rho_rain_rr)

    # density log means
    rho_dry_log = 0
    rho_vapour_log = 0
    rho_cloud_log = 0
    rho_rain_log = 0

    if (!(rho_dry_ll == 0.0) && !(rho_dry_rr == 0.0))
        rho_dry_log = ln_mean(rho_dry_ll, rho_dry_rr)
    end

    if (!(rho_vapour_ll == 0.0) && !(rho_vapour_rr == 0.0))
        rho_vapour_log = ln_mean(rho_vapour_ll, rho_vapour_rr)
    end

    if (!(rho_cloud_ll == 0.0) && !(rho_cloud_rr == 0.0))
        rho_cloud_log = ln_mean(rho_cloud_ll, rho_cloud_rr)
    end

    if (!(rho_rain_ll == 0.0) && !(rho_rain_rr == 0.0))
        rho_rain_log = ln_mean(rho_rain_ll, rho_rain_rr)
    end

    # other averages
    inv_temperature_avg = 0.5f0 * (inv_temperature_ll + inv_temperature_rr)
    inv_temperature_log = inv_ln_mean(inv_temperature_ll, inv_temperature_rr)
    p_int = inv(inv_temperature_avg) * (R_d * rho_dry_avg + R_v * rho_vapour_avg)

    # density flux
    f_vapour = rho_vapour_log * v_dot_n_avg
    f_cloud = rho_cloud_avg * v_dot_n_avg
    f_dry = rho_dry_log * v_dot_n_avg
    f_rain = rho_rain_avg * (v_dot_n_avg - vr_normal_avg)
    f_moist = f_vapour + f_cloud
    f_rho = f_dry + f_moist + f_rain

    # momentum flux
    f_rhov1 = f_rho * v1_avg + p_int * normal_direction[1]
    f_rhov2 = f_rho * v2_avg + p_int * normal_direction[2]

    # energy flux
    f_energy = (c_vd * inv_temperature_log - K_avg) * f_dry +
               (c_vv * inv_temperature_log - K_avg + L_ref) * f_vapour +
               (c_l * inv_temperature_log - K_avg) * (f_cloud + f_rain) +
               (v1_avg * f_rhov1 + v2_avg * f_rhov2)

    return SVector(f_dry, f_moist, f_rain, f_rhov1, f_rhov2, f_energy, 0, 0, 0)
end

@inline function flux_ec_rain(u_ll, u_rr, orientation::Int,
                              equations::CompressibleRainyEulerEquations2D)
    if (orientation == 1)
        return flux_ec_rain(u_ll, u_rr, SVector(1, 0), equations)
    else
        return flux_ec_rain(u_ll, u_rr, SVector(0, 1), equations)
    end
end

# adapted from ShallowWaterEquations2D (Recommended with rain!)
@inline function boundary_condition_simple_slip_wall(u_inner, orientation, direction, x,
                                                     t, surface_flux_function,
                                                     equations::CompressibleRainyEulerEquations2D)
    ## get the appropriate normal vector from the orientation
    if orientation == 1
        u_boundary = SVector(u_inner[1], u_inner[2], u_inner[3], -u_inner[4],
                             u_inner[5], u_inner[6], u_inner[7], u_inner[8], u_inner[9])
    else # orientation == 2
        u_boundary = SVector(u_inner[1], u_inner[2], u_inner[3], u_inner[4],
                             -u_inner[5], u_inner[6], u_inner[7], u_inner[8],
                             u_inner[9])
    end

    # Calculate boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
    end

    return flux
end

# hydrostatic base state residual
function generate_hydrostatic_residual(pressure_lower, humidity_rel0, z, dz,
                                       equations::CompressibleRainyEulerEquations2D)
    # equations constants
    c_pd = equations.c_dry_air_const_pressure
    R_d = equations.R_dry_air
    R_v = equations.R_vapour
    eps = equations.eps
    ref_pressure = equations.ref_pressure
    g = equations.gravity

    function hydrostatic_residual!(residual, guess)
        # variables
        pressure, rho_dry, rho_vapour, temperature = guess

        rho_vs = saturation_vapour_pressure(temperature, equations) / (R_v * temperature)

        # pressure derivative residual approximation
        residual[1] = (pressure - pressure_lower) / dz + (rho_dry + rho_vapour) * g

        # pressure residual
        residual[2] = pressure - temperature * (rho_dry * R_d + rho_vapour * R_v)

        # hydrostatic dry potential temperature residual
        residual[3] = theta_d(z, equations) -
                      temperature * (ref_pressure / pressure)^(R_d / c_pd)

        # humidity residual
        residual[4] = rho_vs * (rho_dry + rho_vapour / eps) * humidity_rel0
        residual[4] -= rho_vapour * (rho_dry + rho_vs / eps)
        residual[4] *= 1000
    end

    return hydrostatic_residual!
end

function generate_perturbation_residual(pressure_hydrostatic, H_init, z,
                                        equations::CompressibleRainyEulerEquations2D)
    # equations constants
    c_pd = equations.c_dry_air_const_pressure
    R_d = equations.R_dry_air
    R_v = equations.R_vapour
    eps = equations.eps
    ref_pressure = equations.ref_pressure

    function perturbation_residual!(residual, guess)
        # variables
        rho_dry, rho_vapour, temperature = guess

        rho_vs = saturation_vapour_pressure(temperature, equations) / (R_v * temperature)
        pressure = (rho_dry * R_d + rho_vapour * R_v) * temperature

        # humidity residual
        residual[1] = rho_vs * (rho_dry + rho_vapour / eps) * H_init
        residual[1] -= rho_vapour * (rho_dry + rho_vs / eps)
        residual[1] *= 30

        # hydrostatic dry potential temperature residual
        residual[2] = theta_d(z, equations) -
                      temperature * (ref_pressure / pressure_hydrostatic)^(R_d / c_pd)

        # pressure residual
        residual[3] = pressure_hydrostatic - pressure
    end

    return perturbation_residual!
end

# hydrostatic dry potential temperature
function theta_d(z, equations::CompressibleRainyEulerEquations2D{RealT}) where {RealT}
    # constants
    c_pd = equations.c_dry_air_const_pressure
    R_d = equations.R_dry_air
    ref_pressure = equations.ref_pressure

    # problem specific constants
    surface_temperature = 283
    surface_pressure = 85000
    stratification = convert(RealT, 1.3e-5)

    # dry potential temperature at surface
    Theta0 = surface_temperature * (ref_pressure / surface_pressure)^(R_d / c_pd)
    # at height z
    theta_d = Theta0 * exp(stratification * z)

    return theta_d
end

# for approximating the dz pressure gradient
struct AtmosphereLayersRainyBubble{RealT <: Real}
    layer_data   :: Matrix{RealT}
    total_height :: RealT
    precision    :: RealT
end

function AtmosphereLayersRainyBubble(equations::CompressibleRainyEulerEquations2D;
                                     total_height, precision = 1, RealT = Float64)
    # constants
    humidity_rel0 = convert(RealT, 0.2)      # hydrostatic relative humidity
    surface_pressure = 85000

    # surface layer with initial guesses for rho_dry, rho_vapour and temperature
    surface_layer = [surface_pressure, convert(RealT, 1.4), convert(RealT, 0.04), 300]

    # allocate layer_data
    n = convert(Int, total_height / precision)
    layer_data = zeros(RealT, n + 1, 4)

    # solve (slightly above) surface layer
    dz = convert(RealT, 0.01)
    z = convert(RealT, 0.01)
    residual_function! = generate_hydrostatic_residual(surface_pressure, humidity_rel0, z,
                                                       dz, equations)
    layer_data[1, :] .= nlsolve(residual_function!, surface_layer).zero

    # adjust to chosen precision
    dz = precision

    # iterate up the atmosphere
    for i in (1:n)
        z += dz
        residual_function! = generate_hydrostatic_residual(layer_data[i, 1], humidity_rel0,
                                                           z, dz, equations)
        guess = deepcopy(layer_data[i, :])
        layer_data[i + 1, :] .= nlsolve(residual_function!, guess, ftol = 1e-10,
                                        iterations = 20).zero
    end

    return AtmosphereLayersRainyBubble{RealT}(layer_data, total_height, precision)
end

function initial_condition_bubble_rainy(x, t, equations::CompressibleRainyEulerEquations2D;
                                        atmosphere_layers = layers)
    # equations constants
    c_vd = equations.c_dry_air_const_volume
    c_vv = equations.c_vapour_const_volume
    ref_L = equations.ref_latent_heat_vap_temp

    # problem specific constants
    humidity_rel_bar = convert(RealT, 0.2)                # background relative humidity field
    humidity_max = 1

    # bubble parameters
    radius_outer, radius_inner = 300, 200      # radii of humidity bubble
    x_center, z_center = 1200, 800      # center of humidity bubble

    # radius relative to bubble center
    r = sqrt((x[1] - x_center)^2 + (x[2] - z_center)^2)

    # humidity definition
    if (r > radius_outer)
        # outside the bubble
        humidity = humidity_rel_bar
    elseif (r > radius_inner)
        # outer layers of the bubble
        humidity = humidity_rel_bar +
                   (humidity_max - humidity_rel_bar) *
                   cos(pi * (r - radius_inner) / (2 * (radius_outer - radius_inner)))^2
    else
        # inner layer
        humidity = humidity_max
    end

    # get atmosphere layer and height information
    @unpack layer_data, total_height, precision = atmosphere_layers
    dz = precision
    z = x[2]
    n = convert(Int, floor((z + eps()) / dz)) + 1
    z_lower = (n - 1) * dz
    z_upper = n * dz

    if (z_lower == total_height)
        z_upper = z_lower + dz
        n = n - 1
    end

    if (n == 0)
        n = 1
    end

    # check height consistency
    if (z > total_height && !(isapprox(z, total_height)))
        error("The atmosphere does not match the simulation domain")
    end

    # get hydrostatic pressures and approximate between lower and upper data point
    pressure_hydrostatic_lower = layer_data[n, 1]
    pressure_hydrostatic_upper = layer_data[n + 1, 1]
    pressure_hydrostatic = (pressure_hydrostatic_upper * (z - z_lower) +
                            pressure_hydrostatic_lower * (z_upper - z)) / dz

    # solve perturbation
    residual_function! = generate_perturbation_residual(pressure_hydrostatic, humidity, z,
                                                        equations)
    rho_dry, rho_vapour, temperature = nlsolve(residual_function!, layer_data[n, 2:4],
                                               ftol = 1e-9, iterations = 20).zero

    energy_density = (c_vd * rho_dry + c_vv * rho_vapour) * temperature + rho_vapour * ref_L

    return SVector(rho_dry, rho_vapour, 0, 0, 0, energy_density, rho_vapour, 0,
                   temperature)
end
end  # muladd end
