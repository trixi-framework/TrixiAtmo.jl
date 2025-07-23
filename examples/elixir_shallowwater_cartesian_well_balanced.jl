###############################################################################
# Entropy-stable DGSEM for the shallow water equations in Cartesian form on the
# cubed sphere: Well-balancedness test for an isolated mountain (checking that
# geopotential height remains constant after one day for zero initial velocity)
###############################################################################

using OrdinaryDiffEqLowStorageRK, Trixi, TrixiAtmo

equations = ShallowWaterEquations3D(gravity = EARTH_GRAVITATIONAL_ACCELERATION,
                                    rotation_rate = EARTH_ROTATION_RATE)

# Create DG solver with polynomial degree = 3 and Wintemeyer et al.'s flux as 
# surface and volume fluxes
polydeg = 3
volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
surface_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)

solver = DGSEM(polydeg = polydeg,
               surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

# Initial condition is atmosphere at rest with constant total geopotential height
function initial_condition_well_balancedness(x, t, equations::ShallowWaterEquations3D)
    # Constant water height
    H = 5960.0
    v1 = v2 = v3 = 0.0

    # Non-constant topography
    b = bottom_topography_isolated_mountain(x)

    return prim2cons(SVector(H, v1, v2, v3, b), equations)
end

initial_condition = initial_condition_well_balancedness

cells_per_dimension = (5, 5)
mesh = P4estMeshCubedSphere2D(cells_per_dimension[1], EARTH_RADIUS, polydeg = polydeg)

# A semidiscretization collects data structures and functions for the spatial discretization
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    metric_terms = MetricTermsInvariantCurl(),
                                    source_terms = source_terms_coriolis_lagrange_multiplier)

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0.0 to π
tspan = (0.0, pi)
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
