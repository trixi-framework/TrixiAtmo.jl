# The inertia-gravity waves in a periodic channel of Baldauf and Brdar (2013), integrated
# with a Rosenbrock-W method and a vertically implicit Jacobian (horizontally explicit,
# vertically implicit time integration)
#
# References:
# - Michael Baldauf, Slavko Brdar (2013)
#   An analytic solution for linear gravity waves in a channel as a test for numerical
#   models using the non-hydrostatic, compressible Euler equations
#   https://doi.org/10.1002/qj.2105

using Trixi, TrixiAtmo

equations = CompressibleEulerPotentialTemperatureEquationsWithGravity2D(c_p = 1004,
                                                                        c_v = 717,
                                                                        gravity = EARTH_GRAVITATIONAL_ACCELERATION)

# We have an isothermal background state with T0 = 250 K, whose speed of sound
# cs = sqrt(gamma * R * T0) is used for the dissipation of the vertical Jacobian below
cs = sqrt(equations.gamma * equations.R * 250)
surface_flux = flux_es
volume_flux = (flux_ec, flux_nonconservative_waruszewski_etal)
polydeg = 3
solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

boundary_conditions = (; y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

coordinates_min = (0.0, 0.0)
coordinates_max = (300_000.0, 10_000.0)
cells_per_dimension = (60, 8)
mesh = P4estMesh(cells_per_dimension, polydeg = polydeg,
                 coordinates_min = coordinates_min, coordinates_max = coordinates_max,
                 periodicity = (true, false), initial_refinement_level = 0)
initial_condition = initial_condition_gravity_waves
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_conditions)
tspan = (0.0, 1800.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_integrals = (entropy,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback)

###############################################################################
# Use a Rosenbrock-W method with the vertically implicit Jacobian, such that the time step
# is restricted by the horizontal dynamics only
jacobian = VerticalJacobianTwoSided(c = cs)

dt = 2.0 # seconds
sol = Trixi.solve(ode, SSPRW32(); dt = dt, callback = callbacks, jacobian = jacobian)
