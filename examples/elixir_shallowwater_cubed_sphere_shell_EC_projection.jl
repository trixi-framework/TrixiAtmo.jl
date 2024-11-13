
using OrdinaryDiffEq
using Trixi
using TrixiAtmo

###############################################################################
# Entropy conservation for the spherical shallow water equations in Cartesian
# form obtained through the projection of the momentum onto the divergence-free
# tangential contravariant vectors

equations = ShallowWaterEquations3D(gravity_constant = 9.81)

# Create DG solver with polynomial degree = 3 and Wintemeyer et al.'s flux as surface flux
polydeg = 3
volume_flux = flux_wintermeyer_etal # flux_fjordholm_etal 
surface_flux = flux_wintermeyer_etal # flux_fjordholm_etal #flux_lax_friedrichs
solver = DGSEM(polydeg = polydeg,
               surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

# Initial condition for a Gaussian density profile with constant pressure
# and the velocity of a rotating solid body
function initial_condition_advection_sphere(x, t, equations::ShallowWaterEquations3D)
    # Gaussian density
    rho = 1.0 + exp(-20 * (x[1]^2 + x[3]^2))

    # Spherical coordinates for the point x
    if sign(x[2]) == 0.0
        signy = 1.0
    else
        signy = sign(x[2])
    end
    # Co-latitude
    colat = acos(x[3] / sqrt(x[1]^2 + x[2]^2 + x[3]^2))
    # Latitude (auxiliary variable)
    lat = -colat + 0.5 * pi
    # Longitude
    r_xy = sqrt(x[1]^2 + x[2]^2)
    if r_xy == 0.0
        phi = pi / 2
    else
        phi = signy * acos(x[1] / r_xy)
    end

    # Compute the velocity of a rotating solid body
    # (alpha is the angle between the rotation axis and the polar axis of the spherical coordinate system)
    v0 = 1.0 # Velocity at the "equator"
    alpha = 0.0 #pi / 4
    v_long = v0 * (cos(lat) * cos(alpha) + sin(lat) * cos(phi) * sin(alpha))
    vlat = -v0 * sin(phi) * sin(alpha)

    # Transform to Cartesian coordinate system
    v1 = -cos(colat) * cos(phi) * vlat - sin(phi) * v_long
    v2 = -cos(colat) * sin(phi) * vlat + cos(phi) * v_long
    v3 = sin(colat) * vlat

    return prim2cons(SVector(rho, v1, v2, v3, 0), equations)
end

initial_condition = initial_condition_advection_sphere

mesh = P4estMeshCubedSphere2D(5, 2.0, polydeg = polydeg,
                              initial_refinement_level = 0)

# A semidiscretization collects data structures and functions for the spatial discretization
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    metric_terms = MetricTermsInvariantCurl(),
                                    source_terms = source_terms_lagrange_multiplier)

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
                                     extra_analysis_integrals = (waterheight, energy_total))

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
