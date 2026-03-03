###############################################################################
# Entropy-stable DGSEM for the 3D shallow water equations in Cartesian form on
# the cubed sphere: Steady geostrophic balance (Case 2, Williamson et al., 1992)
###############################################################################

using OrdinaryDiffEqLowStorageRK, Trixi, TrixiAtmo

###############################################################################
# Parameters

initial_condition = initial_condition_geostrophic_balance
polydeg = 3
cells_per_dimension = (5, 5)
n_saves = 10
tspan = (0.0, 5.0 * SECONDS_PER_DAY)

###############################################################################
# Spatial discretization

mesh = P4estMeshCubedSphere2D(cells_per_dimension[1], EARTH_RADIUS, polydeg = polydeg,
                              element_local_mapping = true)

equations = ShallowWaterEquations3D(gravity = EARTH_GRAVITATIONAL_ACCELERATION,
                                    rotation_rate = EARTH_ROTATION_RATE)

# Create DG solver with polynomial degree = polydeg
volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
surface_flux = (FluxPlusDissipation(flux_wintermeyer_etal, DissipationLocalLaxFriedrichs()),
                flux_nonconservative_wintermeyer_etal)

solver = DGSEM(polydeg = polydeg,
               surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

# Transform the initial condition to the proper set of conservative variables
initial_condition_transformed = transform_initial_condition(initial_condition, equations)

# A semidiscretization collects data structures and functions for the spatial discretization
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_transformed, solver,
                                    source_terms = source_terms_coriolis_lagrange_multiplier,
                                    boundary_conditions = boundary_condition_periodic)

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0 to T
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

# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation
# setup and resets the timers
summary_callback = SummaryCallback()

# The AnalysisCallback allows to analyse the solution in regular intervals and prints the
# results
analysis_callback = AnalysisCallback(semi, interval = 200,
                                     save_analysis = true,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (energy_internal,
                                                                 energy_kinetic, pressure))

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
            dt = 100.0, save_everystep = false, callback = callbacks)
