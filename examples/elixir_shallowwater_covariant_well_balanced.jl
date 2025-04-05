###############################################################################
# Entropy-stable DGSEM for the shallow water equations in covariant form on the
# cubed sphere: Well-balancedness test for an isolated mountain (checking that
# geopotential height remains constant after one day for zero initial velocity)
###############################################################################

using OrdinaryDiffEq, Trixi, TrixiAtmo

###############################################################################
# Parameters

# Initial condition is atmosphere at rest with constant total geopotential height
initial_condition = (x, t, equations) -> SVector(5960.0, 0.0, 0.0, 0.0)

polydeg = 3
cells_per_dimension = (5, 5)
n_saves = 10
tspan = (0.0, 1.0 * SECONDS_PER_DAY)

###############################################################################
# Spatial discretization

mesh = P4estMeshCubedSphere2D(cells_per_dimension[1], EARTH_RADIUS, polydeg = polydeg,
                              element_local_mapping = true)

equations = SplitCovariantShallowWaterEquations2D(EARTH_GRAVITATIONAL_ACCELERATION,
                                                  EARTH_ROTATION_RATE,
                                                  global_coordinate_system = GlobalCartesianCoordinates())

# Use entropy-conservative two-point flux for volume terms, dissipative surface flux
volume_flux = (flux_ec, flux_nonconservative_ec)
surface_flux = (FluxPlusDissipation(flux_ec, DissipationLocalLaxFriedrichs()),
                flux_nonconservative_ec)

# Create DG solver with polynomial degree = polydeg
solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

# Transform the initial condition to the proper set of conservative variables
initial_condition_transformed = transform_initial_condition(initial_condition, equations)

# A semidiscretization collects data structures and functions for the spatial
# discretization. Here, we pass in the additional keyword argument "auxiliary_field" to 
# specify the bottom topography, which is the same as for the standard isolated mountain 
# case.
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_transformed, solver,
                                    source_terms = source_terms_geometric_coriolis,
                                    auxiliary_field = bottom_topography_isolated_mountain)

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0 to T
ode = semidiscretize(semi, tspan)

# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation 
# setup and resets the timers
summary_callback = SummaryCallback()

# The AnalysisCallback allows to analyse the solution in regular intervals and prints the
# results. Note that entropy should be conserved at the semi-discrete level.
analysis_callback = AnalysisCallback(semi, interval = 200,
                                     save_analysis = true,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (entropy,))

# The SaveSolutionCallback allows to save the solution to a file in regular intervals
save_solution = SaveSolutionCallback(dt = (tspan[2] - tspan[1]) / n_saves,
                                     solution_variables = cons2prim_and_vorticity)

# The StepsizeCallback handles the re-calculation of the maximum Δt after each time step.
# Here we take a small CFL value to make the time integration error negligible when
# verifying well balancedness
stepsize_callback = StepsizeCallback(cfl = 0.01)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE 
# solver
callbacks = CallbackSet(summary_callback, analysis_callback, save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

# OrdinaryDiffEq's `solve` method evolves the solution in time and executes the passed 
# callbacks
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 300.0, save_everystep = false, callback = callbacks)
