###############################################################################
# DGSEM for the linear advection equation on the cubed sphere
###############################################################################

using OrdinaryDiffEq, Trixi, TrixiAtmo

###############################################################################
# Spatial discretization

cells_per_dimension = 5
initial_condition = initial_condition_gaussian

equations = CovariantLinearAdvectionEquation2D()

# Create DG solver with polynomial degree = p and a local Lax-Friedrichs flux
solver = DGSEM(polydeg = 3, surface_flux = flux_lax_friedrichs,
               volume_integral = VolumeIntegralWeakForm())

# Create a 2D cubed sphere mesh the size of the Earth
mesh = P4estMeshCubedSphere2D(cells_per_dimension, EARTH_RADIUS, polydeg = 3,
                              initial_refinement_level = 0,
                              element_local_mapping = true)

# Convert initial condition given in terms of zonal and meridional velocity components to 
# one given in terms of contravariant velocity components
initial_condition_transformed = spherical2contravariant(initial_condition, equations)

# A semidiscretization collects data structures and functions for the spatial discretization
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_transformed, solver)

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0 to T
ode = semidiscretize(semi, (0.0, 12 * SECONDS_PER_DAY))

# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation setup
# and resets the timers
summary_callback = SummaryCallback()

# The AnalysisCallback allows to analyse the solution in regular intervals and prints the results
analysis_callback = AnalysisCallback(semi, interval = 10,
                                     save_analysis = true,
                                     extra_analysis_errors = (:conservation_error,))

# The SaveSolutionCallback allows to save the solution to a file in regular intervals
save_solution = SaveSolutionCallback(interval = 10,
                                     solution_variables = cons2cons)

# The StepsizeCallback handles the re-calculation of the maximum Δt after each time step
stepsize_callback = StepsizeCallback(cfl = 0.7)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE solver
callbacks = CallbackSet(summary_callback, analysis_callback, save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

# OrdinaryDiffEq's `solve` method evolves the solution in time and executes the passed callbacks
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 1.0, save_everystep = false, callback = callbacks);

# Print the timer summary
summary_callback()
