###############################################################################
# DGSEM for the linear advection equation on a prismed icosahedral grid
###############################################################################

using OrdinaryDiffEqLowStorageRK, Trixi, TrixiAtmo

###############################################################################
# Spatial discretization

initial_condition = initial_condition_barotropic_instability

equations = CovariantShallowWaterEquations2D(EARTH_GRAVITATIONAL_ACCELERATION,
                                             EARTH_ROTATION_RATE,
                                             global_coordinate_system = GlobalCartesianCoordinates())

###############################################################################
# Build DG solver.

polydeg = 4

dg = DGMulti(element_type = Tri(),
             approximation_type = Polynomial(),
             surface_flux = flux_lax_friedrichs,
             polydeg = polydeg)

###############################################################################
# Build mesh.

initial_refinement_level = 4

mesh = DGMultiMeshTriIcosahedron2D(dg, EARTH_RADIUS;
                                   initial_refinement_level = initial_refinement_level)

# Transform the initial condition to the proper set of conservative variables
initial_condition_transformed = transform_initial_condition(initial_condition, equations)

# A semidiscretization collects data structures and functions for the spatial discretization
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_transformed, dg,
                                    source_terms = source_terms_geometric_coriolis)

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0 to T
tspan = (0.0, 12.0 * SECONDS_PER_DAY)
ode = semidiscretize(semi, tspan)

# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation 
# setup and resets the timers
summary_callback = SummaryCallback()

# The AnalysisCallback allows to analyse the solution in regular intervals and prints the 
# results
analysis_callback = AnalysisCallback(semi, interval = 100,
                                     save_analysis = true,
                                     extra_analysis_errors = (:conservation_error,),
                                     uEltype = real(dg))

# The SaveSolutionCallback allows to save the solution to a file in regular intervals
save_solution = SaveSolutionCallback(interval = 300,
                                     solution_variables = cons2prim_and_vorticity)

# The StepsizeCallback handles the re-calculation of the maximum Δt after each time step
stepsize_callback = StepsizeCallback(cfl = 0.7)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE 
# solver
callbacks = CallbackSet(summary_callback, analysis_callback, save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

# OrdinaryDiffEq's `solve` method evolves the solution in time and executes the passed 
# callbacks
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 1.0, save_everystep = false, callback = callbacks)
