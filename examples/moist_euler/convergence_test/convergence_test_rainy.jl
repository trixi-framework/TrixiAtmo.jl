using OrdinaryDiffEq
using Trixi
using TrixiAtmo
using TrixiAtmo: saturation_vapour_pressure, saturation_vapour_pressure_derivative,
                 saturation_residual, saturation_residual_jacobian,
                 terminal_velocity_rain

function initial_condition_convergence_test_rainy(x, t,
                                                  equations::CompressibleRainyEulerEquations2D)
    # constants
    c_l = equations.c_liquid_water
    c_vd = equations.c_dry_air_const_volume
    c_vv = equations.c_vapour_const_volume
    R_v = equations.R_vapour
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
    rho_vapour = saturation_vapour_pressure(temperature, equations) / (R_v * temperature)
    rho_cloud = rho / c_l * 3000
    rho_rain = rho / c_l * 1000
    rho_moist = rho_vapour + rho_cloud
    rho_dry = rho - rho_moist - rho_rain

    # define matching energydensity with v1 := 1 and v2 := 1 , initially
    energy = (c_vd * rho_dry + c_vv * rho_vapour + c_l * (rho_cloud + rho_rain)) *
             temperature
    energy += rho_vapour * ref_L + rho

    return SVector(rho_dry, rho_moist, rho_rain, rho, rho, energy, rho_vapour, rho_cloud,
                   temperature)
end

function source_terms_convergence_test_rainy(u, x, t,
                                             equations::CompressibleRainyEulerEquations2D)
    # constants
    c_l = equations.c_liquid_water
    c_vd = equations.c_dry_air_const_volume
    c_vv = equations.c_vapour_const_volume
    R_d = equations.R_dry_air
    R_v = equations.R_vapour
    ref_L = equations.ref_latent_heat_vap_temp
    N_0 = equations.rain_water_distr
    v_0 = equations.v_mean_rain

    # help constant for terminal rain velocity derivative ( \Gamma(4.5) / 6 ~= 1.9386213994279082 )
    c_help = v_0 * 1.9386213994279082 * (pi * N_0)^(-0.125)

    # define rho like initial condition
    c = 2.0
    A = 0.1
    L = 2.0
    f = 1.0 / L
    ω = 2 * pi * f
    si, co = sincos(ω * (x[1] + x[2] - t))
    rho = c + A * si
    rho_x = ω * A * co
    rho_t = -rho_x

    # define variables of rho
    temperature = rho + 250.0
    sat_vap_p = saturation_vapour_pressure(temperature, equations)
    rho_vapour = sat_vap_p / (R_v * temperature)
    rho_cloud = rho / c_l * 3000
    rho_rain = rho / c_l * 1000
    rho_moist = rho_vapour + rho_cloud
    rho_dry = rho - rho_moist - rho_rain
    vr = terminal_velocity_rain(rho_moist, rho_rain, equations)

    # define needed derivatives
    sat_vap_p_t = rho_t * saturation_vapour_pressure_derivative(temperature, equations)
    sat_vap_p_x = rho_x * saturation_vapour_pressure_derivative(temperature, equations)

    rho_vapour_t = (sat_vap_p_t * temperature - rho_t * sat_vap_p) / (R_v * temperature^2)
    rho_vapour_x = (sat_vap_p_x * temperature - rho_x * sat_vap_p) / (R_v * temperature^2)

    rho_cloud_t = rho_t / c_l * 3000
    rho_cloud_x = rho_x / c_l * 3000

    rho_rain_t = rho_t / c_l * 1000
    rho_rain_x = rho_x / c_l * 1000

    rho_moist_t = rho_vapour_t + rho_cloud_t
    rho_moist_x = rho_vapour_x + rho_cloud_x

    rho_dry_t = rho_t - rho_moist_t - rho_rain_t
    rho_dry_x = rho_x - rho_moist_x - rho_rain_x

    energy_t = (c_vd * rho_dry_t + c_vv * rho_vapour_t + c_l * (rho_cloud_t + rho_rain_t)) *
               temperature
    energy_t += (c_vd * rho_dry + c_vv * rho_vapour + c_l * (rho_cloud + rho_rain)) * rho_t
    energy_t += rho_vapour_t * ref_L + rho_t

    energy_x = (c_vd * rho_dry_x + c_vv * rho_vapour_x + c_l * (rho_cloud_x + rho_rain_x)) *
               temperature
    energy_x += (c_vd * rho_dry + c_vv * rho_vapour + c_l * (rho_cloud + rho_rain)) * rho_x
    energy_x += rho_vapour_x * ref_L + rho_x

    pressure_x = (rho_dry_x * R_d + rho_vapour_x * R_v) * temperature
    pressure_x += (rho_dry * R_d + rho_vapour * R_v) * rho_x         # temperature_x = rho_x

    vr_x = c_help * 0.125 *
           ((rho_rain_x * rho_moist - rho_rain * rho_moist_x) / (rho_moist + rho_rain)^2)
    vr_x *= (rho_rain / (rho_moist + rho_rain))^(-0.875)

    rhor_vr__x = rho_rain_x * vr + rho_rain * vr_x

    # calculate source terms for manufactured solution
    # density
    S_rho_dry = rho_dry_t + 2.0 * rho_dry_x
    S_rho_moist = rho_moist_t + 2.0 * rho_moist_x
    S_rho_rain = rho_rain_t + 2.0 * rho_rain_x - rhor_vr__x

    # "momentum"
    S_rho_v1 = rho_x + pressure_x - rhor_vr__x
    S_rho_v2 = rho_x + pressure_x - rhor_vr__x

    # "energy"
    S_energy = energy_t + 2.0 * (energy_x + pressure_x) - (c_l * rho_x * rho_rain * vr)
    S_energy -= (c_l * temperature + 1) * rhor_vr__x

    return SVector(S_rho_dry, S_rho_moist, S_rho_rain, S_rho_v1, S_rho_v2, S_energy, 0.0,
                   0.0, 0.0)
end

###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleRainyEulerEquations2D()

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

stage_limiter! = NonlinearSolveDG(saturation_residual, saturation_residual_jacobian,
                                  SVector(7, 8, 9), 1e-9)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false, stage_limiter!),
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks);

# For copy-paste convenience:
#convergence_test("TrixiAtmo.jl/examples/convergence_test_elixirs/convergence_test_rainy.jl", 5)

#=
####################################################################################################
l2
rho_dry             rho_moist           rho_rain            rho_v1              rho_v2              energy_density      rho_vapour          rho_cloud           temperature 

error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     EOC
2.40e-06  -         3.89e-05  -         1.52e-05  -         7.53e-03  -         2.09e-02  -         5.20e+01  -         1.75e-08  -         3.89e-05  -         2.16e-04  - 

1.43e-07  4.07      2.32e-06  4.07      9.72e-07  3.97      2.90e-04  4.70      1.39e-03  3.91      3.14e+00  4.05      1.15e-09  3.93      2.32e-06  4.07      1.42e-05  3.93
8.48e-09  4.08      1.37e-07  4.08      1.20e-07  3.02      1.13e-05  4.68      3.10e-05  5.49      1.78e-01  4.14      2.69e-11  5.42      1.37e-07  4.08      3.32e-07  5.41
5.29e-10  4.00      8.57e-09  4.00      5.27e-09  4.51      6.34e-07  4.16      9.54e-07  5.02      1.39e-02  3.68      1.24e-12  4.43      8.57e-09  4.00      1.54e-08  4.43
3.31e-11  4.00      5.34e-10  4.00      2.02e-10  4.70      3.85e-08  4.04      3.80e-08  4.65      7.79e-04  4.16      7.09e-14  4.13      5.34e-10  4.00      8.78e-10  4.13

mean      4.04      mean      4.04      mean      4.05      mean      4.39      mean      4.77      mean      4.01      mean      4.48      mean      4.04      mean      4.48
----------------------------------------------------------------------------------------------------
linf
rho_dry             rho_moist           rho_rain            rho_v1              rho_v2              energy_density      rho_vapour          rho_cloud           temperature 

error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     EOC
1.68e-05  -         2.74e-04  -         8.69e-05  -         2.34e-02  -         6.25e-02  -         3.85e+02  -         6.00e-08  -         2.74e-04  -         7.43e-04  - 

1.07e-06  3.97      1.73e-05  3.98      4.36e-06  4.32      1.13e-03  4.38      3.11e-03  4.33      2.29e+01  4.07      5.72e-09  3.39      1.73e-05  3.98      7.09e-05  3.39
6.01e-08  4.15      9.73e-07  4.16      3.46e-07  3.66      4.43e-05  4.67      6.53e-05  5.57      1.08e+00  4.40      1.63e-10  5.13      9.73e-07  4.16      2.00e-06  5.15
3.69e-09  4.02      5.98e-08  4.03      2.91e-08  3.57      2.47e-06  4.16      2.37e-06  4.79      9.44e-02  3.52      9.40e-12  4.12      5.97e-08  4.03      1.15e-07  4.12       
2.30e-10  4.01      3.71e-09  4.01      1.47e-09  4.31      1.50e-07  4.04      1.16e-07  4.35      5.51e-03  4.10      5.70e-13  4.04      3.71e-09  4.01      6.97e-09  4.04       

mean      4.04      mean      4.04      mean      3.96      mean      4.31      mean      4.76      mean      4.02      mean      4.17      mean      4.04      mean      4.18       
----------------------------------------------------------------------------------------------------
Dict{Symbol, Any} with 3 entries:
  :variables => ("rho_dry", "rho_moist", "rho_rain", "rho_v1", "rho_v2", "energy_density", "rho_vapour", "rho_cloud", "temperature")
  :l2        => [4.03641, 4.03834, 4.04901, 4.39454, 4.7665, 4.00677, 4.47752, 4.03837, 4.47639]
  :linf      => [4.03824, 4.04306, 3.96283, 4.31327, 4.76078, 4.023, 4.17127, 4.04306, 4.17577]
=#
