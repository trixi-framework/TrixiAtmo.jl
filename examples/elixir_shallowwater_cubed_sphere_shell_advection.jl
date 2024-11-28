
using OrdinaryDiffEq
using Trixi
using TrixiAtmo

###############################################################################
# semidiscretization of the linear advection equation

cells_per_dimension = 5
initial_condition = initial_condition_gaussian_cartesian

# We use the ShallowWaterEquations3D equations structure but modify the rhs! function to
# convert it to a variable-coefficient advection equation
equations = ShallowWaterEquations3D(gravity_constant = 0.0)

# Create DG solver with polynomial degree = 3 and (local) Lax-Friedrichs/Rusanov flux as surface flux
solver = DGSEM(polydeg = 3, surface_flux = flux_lax_friedrichs)

# Source term function to transform the Euler equations into a linear advection equation with variable advection velocity
function source_terms_convert_to_linear_advection(u, du, x, t,
                                                  equations::ShallowWaterEquations3D,
                                                  normal_direction)
    v1 = u[2] / u[1]
    v2 = u[3] / u[1]
    v3 = u[4] / u[1]

    s2 = du[1] * v1 - du[2]
    s3 = du[1] * v2 - du[3]
    s4 = du[1] * v3 - du[4]

    return SVector(0.0, s2, s3, s4, 0.0)
end

# Create a 2D cubed sphere mesh the size of the Earth
mesh = P4estMeshCubedSphere2D(cells_per_dimension, EARTH_RADIUS, polydeg = 3,
                              initial_refinement_level = 0,
                              element_local_mapping = false)

initial_condition_transformed = transform_initial_condition(initial_condition, equations)

# A semidiscretization collects data structures and functions for the spatial discretization
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_transformed, solver,
                                    source_terms = source_terms_convert_to_linear_advection)

###############################################################################
# ODE solvers, callbacks etc.

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
                                     solution_variables = cons2prim)

# The StepsizeCallback handles the re-calculation of the maximum Δt after each time step
stepsize_callback = StepsizeCallback(cfl = 0.7)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE solver
callbacks = CallbackSet(summary_callback, analysis_callback, save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

# OrdinaryDiffEq's `solve` method evolves the solution in time and executes the passed callbacks
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks);

# Print the timer summary
summary_callback()
