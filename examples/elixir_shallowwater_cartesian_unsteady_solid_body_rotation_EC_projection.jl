using OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiAtmo

###############################################################################
# Entropy conservation for the spherical shallow water equations in Cartesian
# form obtained through the projection of the momentum onto the divergence-free
# tangential contravariant vectors

###############################################################################
# Parameters

initial_condition = initial_condition_unsteady_solid_body_rotation
polydeg = 3
cells_per_dimension = (16, 16)
n_saves = 10
tspan = (0.0, 15.0 * SECONDS_PER_DAY)

###############################################################################
# Spatial discretization
equations = ShallowWaterEquations3D(gravity = EARTH_GRAVITATIONAL_ACCELERATION,
                                    rotation_rate = EARTH_ROTATION_RATE)

# Create DG solver with polynomial degree = 3 and Wintemeyer et al.'s flux as surface flux
volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
surface_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
# For provably entropy-stable surface fluxes, use
# surface_flux = (FluxPlusDissipation(flux_wintermeyer_etal, DissipationLaxFriedrichsEntropyVariables()), 
#                 flux_nonconservative_wintermeyer_etal)

solver = DGSEM(polydeg = polydeg,
               surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

# Transform the initial condition to the proper set of conservative variables
initial_condition_transformed = transform_initial_condition(initial_condition, equations)

mesh = P4estMeshCubedSphere2D(cells_per_dimension[1], EARTH_RADIUS, polydeg = polydeg,
                              initial_refinement_level = 0)

# A semidiscretization collects data structures and functions for the spatial discretization
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_transformed, solver,
                                    metric_terms = MetricTermsInvariantCurl(),
                                    source_terms = source_terms_coriolis_lagrange_multiplier)

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0.0 to π
ode = semidiscretize(semi, tspan)

# Clean the initial condition
for element in eachelement(solver, semi.cache)
    for j in eachnode(solver), i in eachnode(solver)
        u0 = Trixi.wrap_array(ode.u0, semi)

        contravariant_normal_vector = Trixi.get_contravariant_vector(3,
                                                                     semi.cache.elements.contravariant_vectors,
                                                                     i, j, element)
        clean_solution_lagrange_multiplier!(u0[:, i, j, element], equations,
                                            contravariant_normal_vector)
    end
end

# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation setup
# and resets the timers
summary_callback = SummaryCallback()

# The AnalysisCallback allows to analyse the solution in regular intervals and prints the results
analysis_callback = AnalysisCallback(semi, interval = 10,
                                     save_analysis = true,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight, energy_total),
                                     analysis_polydeg = polydeg)

# The SaveSolutionCallback allows to save the solution to a file in regular intervals
save_solution = SaveSolutionCallback(dt = (tspan[2] - tspan[1]) / n_saves,
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
            save_everystep = false, callback = callbacks)
