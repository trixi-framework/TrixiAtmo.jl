using OrdinaryDiffEq
using Trixi
using TrixiAtmo
using TrixiAtmo: saturation_vapour_pressure, saturation_vapour_pressure_derivative,
                 saturation_residual, saturation_residual_jacobian



function initial_condition_convergence_test_rainy_no_rain(x, t, equations::CompressibleRainyEulerExplicitEquations2D)
    # needed constants
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
    rho_cloud   = rho / c_l * 4000
    rho_moist   = rho_vapour + rho_cloud
    rho_dry     = rho - rho_moist

    # define matching energydensity with v1 := 1 and v2 := 1 , initially
    energy  = (c_vd * rho_dry + c_vv * rho_vapour + c_l * rho_cloud) * temperature
    energy += rho_vapour * ref_L + rho

    return SVector(rho_dry, rho_vapour, rho_cloud, 0.0, rho, rho, energy)
end


function source_terms_convergence_test_rainy_no_rain(u, x, t, equations::CompressibleRainyEulerExplicitEquations2D)
    # needed constants
    c_l   = equations.c_liquid_water
    c_vd  = equations.c_dry_air_const_volume
    c_vv  = equations.c_vapour_const_volume
    R_d   = equations.R_dry_air
    R_v   = equations.R_vapour
    ref_L = equations.ref_latent_heat_vap_temp

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
    rho_cloud   = rho / c_l * 4000
    rho_moist   = rho_vapour + rho_cloud
    rho_dry     = rho - rho_moist

    # define needed derivatives
    sat_vap_p_t  = rho_t * saturation_vapour_pressure_derivative(temperature, equations)
    sat_vap_p_x  = rho_x * saturation_vapour_pressure_derivative(temperature, equations)

    rho_vapour_t = (sat_vap_p_t * temperature - rho_t * sat_vap_p) / (R_v * temperature^2)
    rho_vapour_x = (sat_vap_p_x * temperature - rho_x * sat_vap_p) / (R_v * temperature^2)

    rho_cloud_t  = rho_t / c_l * 4000
    rho_cloud_x  = rho_x / c_l * 4000

    rho_moist_t  = rho_vapour_t + rho_cloud_t
    rho_moist_x  = rho_vapour_x + rho_cloud_x
    
    rho_dry_t    = rho_t - rho_moist_t
    rho_dry_x    = rho_x - rho_moist_x

    energy_t     = (c_vd * rho_dry_t + c_vv * rho_vapour_t + c_l * rho_cloud_t) * temperature
    energy_t    += (c_vd * rho_dry   + c_vv * rho_vapour   + c_l * rho_cloud  ) * rho_t
    energy_t    += rho_vapour_t * ref_L + rho_t

    energy_x     = (c_vd * rho_dry_x + c_vv * rho_vapour_x + c_l * rho_cloud_x) * temperature
    energy_x    += (c_vd * rho_dry   + c_vv * rho_vapour   + c_l * rho_cloud  ) * rho_x
    energy_x    += rho_vapour_x * ref_L + rho_x

    pressure_x   = (rho_dry_x * R_d + rho_vapour_x * R_v) * temperature
    pressure_x  += (rho_dry   * R_d + rho_vapour   * R_v) * rho_x         # temperature_x = rho_x

    # calculate source terms for manufactured solution
    # density
    S_rho_dry    = rho_dry_t    + 2.0 * rho_dry_x
    S_rho_vapour = rho_vapour_t + 2.0 * rho_vapour_x
    S_rho_cloud  = rho_cloud_t  + 2.0 * rho_cloud_x
    
    # "momentum"
    S_rho_v1    = rho_x       + pressure_x
    S_rho_v2    = rho_x       + pressure_x

    # "energy"
    S_energy    = energy_t    + 2.0 * (energy_x + pressure_x)

    return SVector(S_rho_dry, S_rho_vapour, S_rho_cloud, 0.0, S_rho_v1, S_rho_v2, S_energy)
end


###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleRainyEulerExplicitEquations2D()

initial_condition = initial_condition_convergence_test_rainy_no_rain

solver = DGSEM(polydeg = 3, surface_flux = flux_lax_friedrichs)

coordinates_min = (0.0, 0.0)
coordinates_max = (2.0, 2.0)

cells_per_dimension = (8, 8)

mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_convergence_test_rainy_no_rain)

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

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks);
summary_callback() # print the timer summary

# For copy-paste convenience:
#convergence_test("TrixiAtmo.jl/examples/convergence_test_elixirs/convergence_test_rainy_no_rain_explicit.jl", 4)
#=
####################################################################################################
l2
rho_dry             rho_vapour          rho_cloud           rho_rain            rho_v1              rho_v2              energy_density
error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     EOC
8.21e-07  -         1.47e-07  -         1.76e-05  -         0.00e+00  -         7.48e-03  -         7.48e-03  -         1.88e+01  -
4.37e-08  4.23      4.86e-09  4.92      9.38e-07  4.23      0.00e+00  NaN       2.13e-04  5.13      2.13e-04  5.13      1.01e+00  4.23      
2.69e-09  4.02      1.91e-10  4.67      5.78e-08  4.02      0.00e+00  NaN       7.04e-06  4.92      7.04e-06  4.92      6.20e-02  4.02
1.68e-10  4.00      5.06e-12  5.24      3.61e-09  4.00      0.00e+00  NaN       2.78e-07  4.66      2.78e-07  4.66      3.87e-03  4.00

mean      4.09      mean      4.94      mean      4.08      mean      NaN       mean      4.90      mean      4.90      mean      4.08
----------------------------------------------------------------------------------------------------
linf
rho_dry             rho_vapour          rho_cloud           rho_rain            rho_v1              rho_v2              energy_density
error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     EOC
6.66e-06  -         6.69e-07  -         1.42e-04  -         0.00e+00  -         1.55e-02  -         1.55e-02  -         1.54e+02  -
3.44e-07  4.28      1.96e-08  5.09      7.25e-06  4.30      0.00e+00  NaN       4.91e-04  4.98      4.91e-04  4.98      7.81e+00  4.30
1.93e-08  4.16      7.70e-10  4.67      4.13e-07  4.13      0.00e+00  NaN       1.81e-05  4.76      1.81e-05  4.76      4.43e-01  4.14
1.18e-09  4.03      2.47e-11  4.96      2.53e-08  4.03      0.00e+00  NaN       8.48e-07  4.42      8.48e-07  4.42      2.72e-02  4.03

mean      4.16      mean      4.91      mean      4.15      mean      NaN       mean      4.72      mean      4.72      mean      4.16
----------------------------------------------------------------------------------------------------
Dict{Symbol, Any} with 3 entries:
  :variables => ("rho_dry", "rho_vapour", "rho_cloud", "rho_rain", "rho_v1", "rho_v2", "energy_density")
  :l2        => [4.08558, 4.94267, 4.08233, NaN, 4.90466, 4.90466, 4.08271]
  :linf      => [4.1557, 4.90799, 4.15355, NaN, 4.72056, 4.72056, 4.15543]=#