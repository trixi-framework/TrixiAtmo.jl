using OrdinaryDiffEq
using Trixi
using TrixiAtmo
using TrixiAtmo: saturation_vapour_pressure, saturation_vapour_pressure_derivative,
                 terminal_velocity_rain, RainLimiterDG



function initial_condition_convergence_test_rainy(x, t, equations::CompressibleRainyEulerExplicitEquations2D)
    # constants
    c_l   = equations.c_liquid_water
    c_vd  = equations.c_dry_air_const_volume
    c_vv  = equations.c_vapour_const_volume
    R_v   = equations.R_vapour
    ref_L = equations.ref_latent_heat_vap_temp

    # define rho like in dry convergence test
    c = 2.0
    A = 0.1
    L = 2.0
    f = 1.0 / L
    ω = 2 * pi * f
    rho = c + A * sin(ω * (x[1] + x[2] - t))

    # define variables of rho
    temperature = rho + 250.0
    rho_vapour  = saturation_vapour_pressure(temperature, equations) / (R_v * temperature)
    rho_cloud   = rho / c_l * 3000
    rho_rain    = rho / c_l * 1000
    rho_moist   = rho_vapour + rho_cloud
    rho_dry     = rho - rho_moist - rho_rain

    # define matching energydensity with v1 := 1 and v2 := 1 , initially
    energy  = (c_vd * rho_dry + c_vv * rho_vapour + c_l * (rho_cloud + rho_rain)) * temperature
    energy += rho_vapour * ref_L + rho

    return SVector(rho_dry, rho_vapour, rho_cloud, rho_rain, rho, rho, energy)
end


function source_terms_convergence_test_rainy(u, x, t, equations::CompressibleRainyEulerExplicitEquations2D)
    # constants
    c_l   = equations.c_liquid_water
    c_vd  = equations.c_dry_air_const_volume
    c_vv  = equations.c_vapour_const_volume
    R_d   = equations.R_dry_air
    R_v   = equations.R_vapour
    ref_L = equations.ref_latent_heat_vap_temp
    N_0   = equations.rain_water_distr
    v_0   = equations.v_mean_rain

    # help constant for terminal rain velocity derivative ( \Gamma(4.5) / 6 ~= 1.9386213994279082 )
    c_help   = v_0 * 1.9386213994279082 * (pi * N_0)^(-0.125)

    # define rho like initial condition
    c = 2.0
    A = 0.1
    L = 2.0
    f = 1.0 / L
    ω = 2 * pi * f
    si, co = sincos(ω * (x[1] + x[2] - t))
    rho    = c + A * si
    rho_x  = ω * A * co
    rho_t  = -rho_x

    # define variables of rho
    temperature = rho + 250.0
    sat_vap_p   = saturation_vapour_pressure(temperature, equations)
    rho_vapour  = sat_vap_p / (R_v * temperature)
    rho_cloud   = rho / c_l * 3000
    rho_rain    = rho / c_l * 1000
    rho_moist   = rho_vapour + rho_cloud
    rho_dry     = rho - rho_moist - rho_rain
    vr          = terminal_velocity_rain(rho_moist, rho_rain, equations)

    # define needed derivatives
    sat_vap_p_t  = rho_t * saturation_vapour_pressure_derivative(temperature, equations)
    sat_vap_p_x  = rho_x * saturation_vapour_pressure_derivative(temperature, equations)

    rho_vapour_t = (sat_vap_p_t * temperature - rho_t * sat_vap_p) / (R_v * temperature^2)
    rho_vapour_x = (sat_vap_p_x * temperature - rho_x * sat_vap_p) / (R_v * temperature^2)

    rho_cloud_t  = rho_t / c_l * 3000
    rho_cloud_x  = rho_x / c_l * 3000

    rho_rain_t   = rho_t / c_l * 1000
    rho_rain_x   = rho_x / c_l * 1000

    rho_moist_t  = rho_vapour_t + rho_cloud_t
    rho_moist_x  = rho_vapour_x + rho_cloud_x
    
    rho_dry_t    = rho_t - rho_moist_t - rho_rain_t
    rho_dry_x    = rho_x - rho_moist_x - rho_rain_x

    energy_t     = (c_vd * rho_dry_t + c_vv * rho_vapour_t + c_l * (rho_cloud_t + rho_rain_t)) * temperature
    energy_t    += (c_vd * rho_dry   + c_vv * rho_vapour   + c_l * (rho_cloud   + rho_rain))   * rho_t
    energy_t    += rho_vapour_t * ref_L + rho_t

    energy_x     = (c_vd * rho_dry_x + c_vv * rho_vapour_x + c_l * (rho_cloud_x + rho_rain_x)) * temperature
    energy_x    += (c_vd * rho_dry   + c_vv * rho_vapour   + c_l * (rho_cloud   + rho_rain))   * rho_x
    energy_x    += rho_vapour_x * ref_L + rho_x

    pressure_x   = (rho_dry_x * R_d + rho_vapour_x * R_v) * temperature
    pressure_x  += (rho_dry   * R_d + rho_vapour   * R_v) * rho_x         # temperature_x = rho_x

    vr_x         = c_help * 0.125 * ((rho_rain_x * rho_moist - rho_rain * rho_moist_x) / (rho_moist + rho_rain)^2)
    vr_x        *= (rho_rain / (rho_moist + rho_rain))^(-0.875)

    rhor_vr__x   = rho_rain_x * vr + rho_rain * vr_x

    # calculate source terms for manufactured solution
    # density
    S_rho_dry    = rho_dry_t    + 2.0 * rho_dry_x
    S_rho_vapour = rho_vapour_t + 2.0 * rho_vapour_x
    S_rho_cloud  = rho_cloud_t  + 2.0 * rho_cloud_x
    S_rho_rain   = rho_rain_t   + 2.0 * rho_rain_x  - rhor_vr__x
    
    # "momentum"
    S_rho_v1    = rho_x       + pressure_x - rhor_vr__x
    S_rho_v2    = rho_x       + pressure_x - rhor_vr__x

    # "energy"
    S_energy    = energy_t    + 2.0 * (energy_x + pressure_x) - (c_l * rho_x * rho_rain * vr)
    S_energy   -= (c_l * temperature + 1) * rhor_vr__x

    return SVector(S_rho_dry, S_rho_vapour, S_rho_cloud, S_rho_rain, S_rho_v1, S_rho_v2, S_energy)
end


###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleRainyEulerExplicitEquations2D()

initial_condition = initial_condition_convergence_test_rainy

solver = DGSEM(polydeg = 3, surface_flux = flux_lax_friedrichs)

coordinates_min = (0.0, 0.0)
coordinates_max = (2.0, 2.0)

cells_per_dimension = (6, 6)

mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_convergence_test_rainy)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 2.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 1000,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        save_solution,
                        stepsize_callback)

stage_limiter! = RainLimiterDG()

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false, stage_limiter!),
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks);
summary_callback() # print the timer summary

# For copy-paste convenience:
#convergence_test("TrixiAtmo.jl/examples/convergence_test_elixirs/convergence_test_rainy_explicit.jl", 5)

#=
####################################################################################################
l2
rho_dry             rho_vapour          rho_cloud           rho_rain            rho_v1              rho_v2              energy_density
error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     EOC
2.39e-06  -         1.19e-07  -         3.86e-05  -         1.51e-05  -         7.49e-03  -         2.06e-02  -         5.17e+01  -
1.43e-07  4.06      3.43e-09  5.12      2.32e-06  4.06      9.69e-07  3.96      2.89e-04  4.70      1.38e-03  3.90      3.13e+00  4.04      
8.44e-09  4.09      8.83e-11  5.28      1.37e-07  4.08      1.20e-07  3.02      1.13e-05  4.68      3.09e-05  5.48      1.78e-01  4.14
5.28e-10  4.00      3.28e-12  4.75      8.55e-09  4.00      5.28e-09  4.51      6.34e-07  4.16      9.54e-07  5.02      1.39e-02  3.68
3.32e-11  3.99      7.50e-14  5.45      5.35e-10  4.00      2.02e-10  4.71      3.85e-08  4.04      3.80e-08  4.65      7.80e-04  4.16

mean      4.03      mean      5.15      mean      4.03      mean      4.05      mean      4.39      mean      4.76      mean      4.00
----------------------------------------------------------------------------------------------------
linf
rho_dry             rho_vapour          rho_cloud           rho_rain            rho_v1              rho_v2              energy_density
error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     EOC
1.68e-05  -         5.92e-07  -         2.69e-04  -         8.52e-05  -         2.33e-02  -         6.19e-02  -         3.71e+02  -
1.06e-06  3.99      1.25e-08  5.57      1.72e-05  3.97      4.33e-06  4.30      1.12e-03  4.37      3.08e-03  4.33      2.27e+01  4.03
5.98e-08  4.15      2.97e-10  5.39      9.68e-07  4.15      3.45e-07  3.65      4.43e-05  4.67      6.53e-05  5.56      1.08e+00  4.39
3.68e-09  4.02      1.28e-11  4.53      5.96e-08  4.02      2.91e-08  3.57      2.47e-06  4.16      2.37e-06  4.79      9.43e-02  3.52
2.30e-10  4.00      3.03e-13  5.40      3.71e-09  4.00      1.47e-09  4.31      1.50e-07  4.04      1.16e-07  4.36      5.52e-03  4.09

mean      4.04      mean      5.22      mean      4.04      mean      3.96      mean      4.31      mean      4.76      mean      4.01
----------------------------------------------------------------------------------------------------
Dict{Symbol, Any} with 3 entries:
  :variables => ("rho_dry", "rho_vapour", "rho_cloud", "rho_rain", "rho_v1", "rho_v2", "energy_density")
  :l2        => [4.03455, 5.14974, 4.0349, 4.04611, 4.39262, 4.76286, 4.0041]
  :linf      => [4.03983, 5.22416, 4.03676, 3.95607, 4.31093, 4.75739, 4.00916]
=#