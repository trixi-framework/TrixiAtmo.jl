###############################################################################
# DGSEM for the shallow water equations on the cubed sphere
###############################################################################

using OrdinaryDiffEq, Trixi, TrixiAtmo

###############################################################################
# Problem definition

function initial_condition_geostrophic_balance(x, t,
                                               equations::CovariantShallowWaterEquations2D)
    (; gravitational_acceleration, rotation_rate) = equations

    radius = sqrt(x[1]^2 + x[2]^2 + x[3]^2)
    lat = asin(x[3] / radius)

    # compute zonal and meridional components of the velocity
    V = convert(eltype(x), 2π) * radius / (12 * SECONDS_PER_DAY)
    v_lon, v_lat = V * cos(lat), zero(eltype(x))

    # compute geopotential height 
    h = 1 / gravitational_acceleration *
        (2.94f4 - (radius * rotation_rate * V + 0.5f0 * V^2) * (sin(lat))^2)

    # convert to conservative variables
    return SVector(h, h * v_lon, h * v_lat)
end

###############################################################################
# Spatial discretization

polydeg = 3
cells_per_dimension = 4

element_local_mapping = true

mesh = P4estMeshCubedSphere2D(5, EARTH_RADIUS, polydeg = polydeg,
                              initial_refinement_level = 0,
                              element_local_mapping = element_local_mapping)

equations = CovariantShallowWaterEquations2D(EARTH_GRAVITATIONAL_ACCELERATION,
                                             EARTH_ROTATION_RATE)

initial_condition = initial_condition_geostrophic_balance
source_terms = source_terms_convergence_test

# Standard weak-form volume integral
volume_integral = VolumeIntegralWeakForm()

# Create DG solver with polynomial degree = p and a local Lax-Friedrichs flux
solver = DGSEM(polydeg = polydeg, surface_flux = flux_lax_friedrichs,
               volume_integral = volume_integral)

# A semidiscretization collects data structures and functions for the spatial discretization
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

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
stepsize_callback = StepsizeCallback(cfl = 0.4)

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
