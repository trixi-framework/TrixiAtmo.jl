###############################################################################
# DGSEM for the linear advection equation on the cubed sphere
###############################################################################

using OrdinaryDiffEq, Trixi, TrixiAtmo

###############################################################################
# Problem definition

function initial_condition_advection_earth(x, t, ::CovariantLinearAdvectionEquation2D)
    RealT = eltype(x)

    # set parameters
    a = EARTH_RADIUS  # radius of the sphere in metres
    V = convert(RealT, 2π) * a / (12 * SECONDS_PER_DAY)  # speed of rotation
    alpha = convert(RealT, π / 4)  # angle of rotation
    h_0 = 1000.0f0  # bump height in metres
    b_0 = 5.0f0 / (a^2)  # bump width
    lon_0, lat_0 = convert(RealT, 3π / 2), 0.0f0  # initial bump location

    # convert initial position to Cartesian coordinates
    x_0 = SVector(a * cos(lat_0) * cos(lon_0), a * cos(lat_0) * sin(lon_0), a * sin(lat_0))

    # compute Gaussian bump profile
    h = h_0 * exp(-b_0 * ((x[1] - x_0[1])^2 + (x[2] - x_0[2])^2 + (x[3] - x_0[3])^2))

    # get zonal and meridional components of the velocity
    lon, lat = atan(x[2], x[1]), asin(x[3] / a)
    v_lon = V * (cos(lat) * cos(alpha) + sin(lat) * cos(lon) * sin(alpha))
    v_lat = -V * sin(lon) * sin(alpha)

    return SVector(h, v_lon, v_lat)
end

###############################################################################
# Spatial discretization

polydeg = 3
cells_per_dimension = 2

element_local_mapping = true

mesh = P4estMeshCubedSphere2D(5, EARTH_RADIUS, polydeg = polydeg,
                              initial_refinement_level = 0,
                              element_local_mapping = element_local_mapping)

equations = CovariantLinearAdvectionEquation2D()

initial_condition = initial_condition_advection_earth

# Local Lax-Friedrichs surface flux
surface_flux = flux_lax_friedrichs

# Standard weak-form volume integral
volume_integral = VolumeIntegralWeakForm()

# Create DG solver with polynomial degree = p and a local Lax-Friedrichs flux
solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
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
