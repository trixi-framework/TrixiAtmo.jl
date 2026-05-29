###############################################################################
# DGSEM for the linear advection equation on a prismed icosahedral grid
###############################################################################

using OrdinaryDiffEqLowStorageRK, Trixi, TrixiAtmo

###############################################################################
# Spatial discretization

initial_condition = initial_condition_unsteady_solid_body_rotation

equations = SplitCovariantShallowWaterEquations2D(EARTH_GRAVITATIONAL_ACCELERATION,
                                             EARTH_ROTATION_RATE,
                                             global_coordinate_system = GlobalCartesianCoordinates())

###############################################################################
# Build DG solver.

polydeg = 3

volume_flux = (flux_ec, flux_nonconservative_ec)
surface_flux = (FluxPlusDissipation(flux_ec, DissipationLocalLaxFriedrichs()),
                flux_nonconservative_surface_simplified)

dg = DGMulti(polydeg = polydeg, element_type = Tri(), approximation_type = SBP(),
             surface_integral = SurfaceIntegralWeakForm(surface_flux),
             volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Build mesh.

initial_refinement_level = 3

mesh = DGMultiMeshTriIcosahedron2D(dg, EARTH_RADIUS;
                                   initial_refinement_level = initial_refinement_level)

# Transform the initial condition to the proper set of conservative variables
initial_condition_transformed = transform_initial_condition(initial_condition, equations)

# A semidiscretization collects data structures and functions for the spatial discretization
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_transformed, dg,
                                    metric_terms = MetricTermsCovariantSphere(
                                        christoffel_symbols = ChristoffelSymbolsCollocationDerivative()
                                    ),
                                    source_terms = source_terms_geometric_coriolis,
                                    auxiliary_field = bottom_topography_unsteady_solid_body_rotation)

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0 to T
tspan = (0.0, 15.0 * SECONDS_PER_DAY)
ode = semidiscretize(semi, tspan)

# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation 
# setup and resets the timers
summary_callback = SummaryCallback()

# The AnalysisCallback allows to analyse the solution in regular intervals and prints the 
# results
analysis_callback = AnalysisCallback(semi, interval = 100,
                                     save_analysis = true,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (entropy,),
                                     uEltype = real(dg))

# The SaveSolutionCallback allows to save the solution to a file in regular intervals
save_solution = SaveSolutionCallback(interval = 300,
                                     solution_variables = cons2prim_and_vorticity)

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
            dt = 1.0, save_everystep = false, callback = callbacks)
