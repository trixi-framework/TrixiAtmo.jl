###############################################################################
# DGSEM for the shallow water equations in covariant form on the cubed sphere
###############################################################################

using OrdinaryDiffEq, Trixi, TrixiAtmo

###############################################################################
# Parameters

initial_condition = initial_condition_geostrophic_balance
polydeg = 3
cells_per_dimension = 5
n_saves = 10
tspan = (0.0, 5.0 * SECONDS_PER_DAY)

###############################################################################
# Spatial discretization

mesh = P4estMeshCubedSphere2D(cells_per_dimension, EARTH_RADIUS, polydeg = polydeg,
                              initial_refinement_level = 0,
                              element_local_mapping = true)

equations = CovariantShallowWaterEquations2D(EARTH_GRAVITATIONAL_ACCELERATION,
                                             EARTH_ROTATION_RATE,
                                             global_coordinate_system = GlobalSphericalCoordinates())

# The covariant shallow water equations are treated as a nonconservative system in order to 
# handle flux-differencing formulations of the covariant derivative. With 
# VolumeIntegralWeakForm, there are actually no nonconservative terms, but we must still
# pass a no-op function "flux_nonconservative_weak_form" as the nonconservative surface flux
surface_flux = (flux_lax_friedrichs, flux_nonconservative_weak_form)

# Create DG solver with polynomial degree = p
solver = DGSEM(polydeg = polydeg, volume_integral = VolumeIntegralWeakForm(),
               surface_flux = surface_flux)

# Transform the initial condition to the proper set of conservative variables
initial_condition_transformed = transform_initial_condition(initial_condition, equations)

# A semidiscretization collects data structures and functions for the spatial discretization
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_transformed, solver,
                                    source_terms = source_terms_weak_form)

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0 to T
ode = semidiscretize(semi, tspan)

# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation 
# setup and resets the timers
summary_callback = SummaryCallback()

# The AnalysisCallback allows to analyse the solution in regular intervals and prints the 
# results
analysis_callback = AnalysisCallback(semi, interval = 200,
                                     save_analysis = true,
                                     extra_analysis_errors = (:conservation_error,))

# The SaveSolutionCallback allows to save the solution to a file in regular intervals
save_solution = SaveSolutionCallback(dt = (tspan[2] - tspan[1]) / n_saves,
                                     solution_variables = cons2cons)

# The StepsizeCallback handles the re-calculation of the maximum Δt after each time step
stepsize_callback = StepsizeCallback(cfl = 0.4)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE 
# solver
callbacks = CallbackSet(summary_callback, analysis_callback, save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

# OrdinaryDiffEq's `solve` method evolves the solution in time and executes the passed 
# callbacks
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 100.0, save_everystep = false, callback = callbacks);

# Print the timer summary
summary_callback()
