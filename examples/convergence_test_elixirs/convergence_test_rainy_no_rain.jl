using OrdinaryDiffEq
using Trixi
using TrixiAtmo
using TrixiAtmo: saturation_vapour_pressure, saturation_vapour_pressure_derivative,
                 saturation_residual_custom, saturation_residual_jacobian_custom



function initial_condition_convergence_test_rainy_no_rain(x, t, equations::CompressibleRainyEulerEquations2D)
    # needed constants
    c_l      = equations.c_liquid_water
    c_vd     = equations.c_dry_air_const_volume
    c_vv     = equations.c_vapour_const_volume
    R_v      = equations.R_vapour
    ref_L    = equations.ref_latent_heat_vap_temp

    # define rho like in dry convergence test
    c = 2.0
    A = 0.1
    L = 2.0
    f = 1.0 / L
    ω = 2 * pi * f
    rho = c + A * sin(ω * (x[1] + x[2] - t))

    # define variables of rho
    temperature = rho
    rho_vapour  = saturation_vapour_pressure(temperature, equations) / (R_v * temperature)
    rho_cloud   = rho / c_l
    rho_moist   = rho_vapour + rho_cloud
    rho_dry     = rho - rho_moist

    # define matching energydensity with v1 := 1 and v2 := 1 , initially
    energy  = (c_vd * rho_dry + c_vv * rho_vapour + c_l * rho_cloud) * temperature
    energy += rho_vapour * ref_L + rho

    return SVector(rho_dry, rho_moist, 0.0, rho, rho, energy, rho_vapour, rho_cloud, temperature)
end


function source_terms_convergence_test_rainy_no_rain(u, x, t, equations::CompressibleRainyEulerEquations2D)
    # needed constants
    c_l      = equations.c_liquid_water
    c_vd     = equations.c_dry_air_const_volume
    c_vv     = equations.c_vapour_const_volume
    R_d      = equations.R_dry_air
    R_v      = equations.R_vapour
    ref_L    = equations.ref_latent_heat_vap_temp

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
    temperature = rho
    sat_vap_p   = saturation_vapour_pressure(temperature, equations)
    rho_vapour  = sat_vap_p / (R_v * temperature)
    rho_cloud   = rho / c_l
    rho_moist   = rho_vapour + rho_cloud
    rho_dry     = rho - rho_moist

    # define needed derivatives
    sat_vap_p_t  = rho_t * saturation_vapour_pressure_derivative(temperature, equations)
    sat_vap_p_x  = rho_x * saturation_vapour_pressure_derivative(temperature, equations)

    rho_vapour_t = (sat_vap_p_t * temperature - rho_t * sat_vap_p) / (R_v * temperature^2)
    rho_vapour_x = (sat_vap_p_x * temperature - rho_x * sat_vap_p) / (R_v * temperature^2)

    rho_cloud_t  = rho_t / c_l
    rho_cloud_x  = rho_x / c_l

    rho_moist_t  = rho_vapour_t + rho_cloud_t
    rho_moist_x  = rho_vapour_x + rho_cloud_x
    
    rho_dry_t    = rho_t - rho_moist_t
    rho_dry_x    = rho_x - rho_moist_x

    energy_t     = (c_vd * rho_dry_t + c_vv * rho_vapour_t + rho_t) * rho
    energy_t    += (c_vd * rho_dry   + c_vv * rho_vapour   + rho  ) * rho_t
    energy_t    += rho_vapour_t * ref_L + rho_t

    energy_x     = (c_vd * rho_dry_x + c_vv * rho_vapour_x + rho_x) * rho
    energy_x    += (c_vd * rho_dry   + c_vv * rho_vapour   + rho  ) * rho_x
    energy_x    += rho_vapour_x * ref_L + rho_x

    pressure_x   = (rho_dry_x * R_d + rho_vapour_x * R_v) * temperature
    pressure_x  += (rho_dry   * R_d + rho_vapour   * R_v) * rho_x         # temperature_x = rho_x

    # calculate source terms for manufactured solution
    # density
    S_rho_dry   = rho_dry_t   + 2.0 * rho_dry_x
    S_rho_moist = rho_moist_t + 2.0 * rho_moist_x
    
    # "momentum"
    S_rho_v1    = rho_x       + pressure_x
    S_rho_v2    = rho_x       + pressure_x

    # "energy"
    S_energy    = energy_t    + 2.0 * (energy_x + pressure_x)

    return SVector(S_rho_dry, S_rho_moist, 0.0, S_rho_v1, S_rho_v2, S_energy, 0.0, 0.0, 0.0)
end


###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleRainyEulerEquations2D()

initial_condition = initial_condition_convergence_test_rainy_no_rain

solver = DGSEM(polydeg = 3, surface_flux = flux_lax_friedrichs)

coordinates_min = (0.0, 0.0)
coordinates_max = (2.0, 2.0)

cells_per_dimension = (16, 16)

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

save_solution = SaveSolutionCallback(interval = analysis_interval,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        save_solution,
                        stepsize_callback)

stage_limiter! = NonlinearSolveDG(saturation_residual_custom, saturation_residual_jacobian_custom, SVector(7, 8, 9), 1e-9)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false, stage_limiter!),
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks);
summary_callback() # print the timer summary

# For copy-paste convenience:
#convergence_test("TrixiAtmo.jl/examples/convergence_test_elixirs/convergence_test_rainy_no_rain.jl", 3)