using OrdinaryDiffEqLowStorageRK
using Trixi, TrixiAtmo
using TrixiAtmo: source_terms_convergence_test_moist,
                 initial_condition_convergence_test_moist

###############################################################################
# semidiscretization of the compressible Euler equations

c_pd = 1004 # specific heat at constant pressure for dry air
c_vd = 717  # specific heat at constant volume for dry air
c_pv = 1885 # specific heat at constant pressure for moist air
c_vv = 1424 # specific heat at constant volume for moist air
equations = CompressibleMoistEulerEquations2D(c_pd = c_pd, c_vd = c_vd, c_pv = c_pv,
                                              c_vv = c_vv,
                                              gravity = EARTH_GRAVITATIONAL_ACCELERATION)

initial_condition = initial_condition_convergence_test_moist

polydeg = 4

surface_flux = flux_chandrashekar
volume_flux = flux_chandrashekar

volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(polydeg, surface_flux, volume_integral)

coordinates_min = (0.0, 0.0)
coordinates_max = (2.0, 2.0)

cells_per_dimension = (4, 4)

mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max, periodicity = true)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_convergence_test_moist,
                                    boundary_conditions = boundary_condition_periodic)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 2.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 100,
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
            maxiters = 1.0e7,
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks)
