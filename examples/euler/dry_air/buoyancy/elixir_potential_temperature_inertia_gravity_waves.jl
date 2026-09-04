using OrdinaryDiffEqSSPRK
using Trixi, TrixiAtmo

equations = CompressibleEulerPotentialTemperatureEquationsWithGravity2D(c_p = 1004,
                                                                        c_v = 717,
                                                                        gravity = EARTH_GRAVITATIONAL_ACCELERATION)

# We have an isothermal background state with T0 = 250 K.
# The reference speed of sound can be computed as:
# cs = sqrt(gamma * R * T0)
cs = sqrt(equations.gamma * equations.R * 250)
surface_flux = (FluxLMARS(cs), flux_zero)
volume_flux = (flux_tec, flux_nonconservative_waruszewski_etal)
polydeg = 3
solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

boundary_conditions = (; x_neg = boundary_condition_periodic,
                       x_pos = boundary_condition_periodic,
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

coordinates_min = (0.0, 0.0)
coordinates_max = (300_000.0, 10_000.0)
cells_per_dimension = (60, 8)
mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                      periodicity = (true, false))
source_terms = nothing
initial_condition = initial_condition_gravity_waves
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms,
                                    boundary_conditions = boundary_conditions)
tspan = (0.0, 1800.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_integrals = (entropy,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        stepsize_callback)

sol = solve(ode,
            SSPRK43(thread = Trixi.Threaded());
            maxiters = 1.0e7,
            dt = 1e-1, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks, adaptive = false)
