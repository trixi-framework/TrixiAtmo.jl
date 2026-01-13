using OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiAtmo
using TrixiAtmo: source_terms_no_phase_change, saturation_residual,
                 saturation_residual_jacobian, NonlinearSolveDG,
                 cons2eq_pot_temp, flux_chandrashekar,
                 boundary_condition_simple_slip_wall, initial_condition_moist_bubble
using NLsolve: nlsolve

equations_moist = CompressibleMoistEulerEquations2D(c_pd = 1004, c_vd = 717, c_pv = 1885,
                                                    c_vv = 1424,
                                                    gravity = EARTH_GRAVITATIONAL_ACCELERATION)

# Create background atmosphere data set
atmosphere_data = AtmosphereLayers(equations_moist)

# Create the initial condition with the initial data set
function initial_condition_moist(x, t, equations::CompressibleRainyEulerEquations2D)
    rho, rho_v1, rho_v2, rho_E, rho_qv, rho_ql, T_loc = initial_condition_moist_bubble(x, t,
                                                                                       equations_moist,
                                                                                       atmosphere_data)
    return SVector(rho - rho_qv - rho_ql, rho_qv + rho_ql, 0, rho_v1, rho_v2, rho_E,
                   rho_qv, rho_ql, T_loc)
end

###############################################################################
# semidiscretization of the compressible rainy Euler equations

equations = CompressibleRainyEulerEquations2D()

boundary_conditions = (x_neg = boundary_condition_periodic,
                       x_pos = boundary_condition_periodic,
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

polydeg = 3

surface_flux = flux_LMARS
volume_flux = flux_ec_rain

volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(polydeg, surface_flux, volume_integral)

coordinates_min = (0.0, 0.0)
coordinates_max = (20_000.0, 10_000.0)

cells_per_dimension = (128, 64)
mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                      periodicity = (true, false))

semi = SemidiscretizationHyperbolic(mesh, equations,
                                    initial_condition_moist, solver,
                                    source_terms = source_terms_no_phase_change,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 1000.0)

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000

analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:entropy_conservation_error,))

alive_callback = AliveCallback(analysis_interval = 1000)

save_solution = SaveSolutionCallback(interval = 1000,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = "out",
                                     solution_variables = cons2eq_pot_temp)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        stepsize_callback)

stage_limiter! = NonlinearSolveDG(saturation_residual, saturation_residual_jacobian,
                                  SVector(7, 8, 9))

###############################################################################
# run the simulation
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false, stage_limiter!),
            maxiters = 1.0e7,
            dt = 1.0,
            save_everystep = false, callback = callbacks);
