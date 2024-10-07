###############################################################################
# DGSEM for the shallow water equations on the cubed sphere
###############################################################################

using OrdinaryDiffEq, Trixi, TrixiAtmo

###############################################################################
# Problem definition

function initial_condition_rossby_haurwitz(x, t,
                                           equations::CovariantShallowWaterEquations2D)
    (; gravitational_acceleration, rotation_rate) = equations

    radius = sqrt(x[1]^2 + x[2]^2 + x[3]^2)
    lon = atan(x[2], x[1])
    lat = asin(x[3] / radius)

    h_0 = 8.0f3
    K = 7.848f-6
    R = 4.0f0

    A = 0.5f0 * K * (2 * rotation_rate + K) * (cos(lat))^2 +
        0.25f0 * K^2 * (cos(lat))^(2 * R) *
        ((R + 1) * (cos(lat))^2 +
         (2 * R^2 - R - 2) - 2 * R^2 / ((cos(lat))^2))
    B = 2 * (rotation_rate + K) * K / ((R + 1) * (R + 2)) * (cos(lat))^R *
        ((R^2 + 2R + 2) - (R + 1)^2 * (cos(lat))^2)
    C = 0.25f0 * K^2 * (cos(lat))^(2 * R) * ((R + 1) * (cos(lat))^2 - (R + 2))

    h = h_0 +
        (1 / gravitational_acceleration) *
        (radius^2 * A + radius^2 * B * cos(R * lon) + radius^2 * C * cos(2 * R * lon))

    v_lon = radius * K * cos(lat) +
            radius * K * (cos(lat))^(R - 1) * (R * (sin(lat))^2 - (cos(lat))^2) *
            cos(R * lon)
    v_lat = -radius * K * R * (cos(lat))^(R - 1) * sin(lat) * sin(R * lon)

    # convert to conservative variables
    return SVector(h, h * v_lon, h * v_lat)
end

###############################################################################
# Spatial discretization

polydeg = 5
cells_per_dimension = 5
element_local_mapping = true

mesh = P4estMeshCubedSphere2D(cells_per_dimension, EARTH_RADIUS, polydeg = polydeg,
                              initial_refinement_level = 0,
                              element_local_mapping = element_local_mapping)

equations = CovariantShallowWaterEquations2D(EARTH_GRAVITATIONAL_ACCELERATION,
                                             EARTH_ROTATION_RATE)

initial_condition = initial_condition_rossby_haurwitz
source_terms = source_terms_convergence_test
tspan = (0.0, 7 * SECONDS_PER_DAY)

# Standard weak-form volume integral
volume_integral = VolumeIntegralWeakForm()

# Create DG solver with polynomial degree = p and a local Lax-Friedrichs flux
solver = DGSEM(polydeg = polydeg, surface_flux = flux_lax_friedrichs,
               volume_integral = volume_integral)

# A semidiscretization collects data structures and functions for the spatial discretization
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms)

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0 to T
ode = semidiscretize(semi, tspan)

# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation setup
# and resets the timers
summary_callback = SummaryCallback()

# The AnalysisCallback allows to analyse the solution in regular intervals and prints the results
analysis_callback = AnalysisCallback(semi, interval = 100,
                                     save_analysis = true,
                                     extra_analysis_errors = (:conservation_error,))

# The SaveSolutionCallback allows to save the solution to a file in regular intervals
save_solution = SaveSolutionCallback(dt = (tspan[2] - tspan[1]) / 50,
                                     solution_variables = cons2cons)

# The StepsizeCallback handles the re-calculation of the maximum Δt after each time step
stepsize_callback = StepsizeCallback(cfl = 0.4)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE solver
callbacks = CallbackSet(summary_callback, analysis_callback, save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

# OrdinaryDiffEq's `solve` method evolves the solution in time and executes the passed callbacks
tsave = LinRange(tspan..., 100)
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 100.0, save_everystep = false, saveat = tsave, callback = callbacks);

# Print the timer summary
summary_callback()
