# An idealized baroclinic instability test case
# For optimal results consider increasing the resolution to 16x16x8 trees per cube face.
#
# This elixir takes about 8 hours, using 16 threads of an AMD Ryzen 7 7800X3D.
#
# References:
# - Paul A. Ullrich, Thomas Melvin, Christiane Jablonowski, Andrew Staniforth (2013)
#   A proposed baroclinic wave test case for deep- and shallow-atmosphere dynamical cores
#   https://doi.org/10.1002/qj.2241

using OrdinaryDiffEqSSPRK
using Trixi, TrixiAtmo

equations = CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = 1004,
                                                                        c_v = 717,
                                                                        gravity = EARTH_GRAVITATIONAL_ACCELERATION)

initial_condition = initial_condition_baroclinic_instability

boundary_conditions = (; inside = boundary_condition_slip_wall,
                       outside = boundary_condition_slip_wall)

polydeg = 5
surface_flux = (FluxLMARS(340), flux_zero)
volume_flux = (flux_tec, flux_nonconservative_souza_etal)

solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))
trees_per_cube_face = (8, 4)
mesh = P4estMeshCubedSphereTopography(trees_per_cube_face..., EARTH_RADIUS, 30000,
                                      polydeg = polydeg, initial_refinement_level = 0)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_coriolis,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.
T = 10 # 10 days
tspan = (0.0, T * SECONDS_PER_DAY) # time in seconds for 10 days

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 100, save_initial_solution = true,
                                     save_final_solution = true)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution)

###############################################################################
# Use a Runge-Kutta method with automatic (error based) time step size control
# Enable threading of the RK method for better performance on multiple threads
tol = 1e-6
sol = solve(ode,
            SSPRK43(thread = Trixi.Threaded());
            abstol = tol, reltol = tol, ode_default_options()...,
            callback = callbacks)
