###############################################################################
# DGSEM for the shallow water equations on the cubed sphere
###############################################################################

using OrdinaryDiffEq, Trixi, TrixiAtmo, QuadGK

###############################################################################
# Problem definition
@inline function galewsky_velocity(θ, u_0, θ_0, θ_1)
    if (θ_0 < θ) && (θ < θ_1)
        u = u_0 / exp(-4 / (θ_1 - θ_0)^2) * exp(1 / (θ - θ_0) * 1 / (θ - θ_1))
    else
        u = zero(θ)
    end
    return u
end

@inline function galewsky_integrand(θ, u_0, θ_0, θ_1, a,
                                    equations::CovariantShallowWaterEquations2D)
    (; rotation_rate) = equations
    u = galewsky_velocity(θ, u_0, θ_0, θ_1)
    return u * (2 * rotation_rate * sin(θ) + u * tan(θ) / a)
end

function initial_condition_barotropic_instability(x, t,
                                                  equations::CovariantShallowWaterEquations2D)
    (; gravitational_acceleration, rotation_rate) = equations
    realT = eltype(x)
    radius = sqrt(x[1]^2 + x[2]^2 + x[3]^2)
    lon = atan(x[2], x[1])
    lat = asin(x[3] / radius)

    # compute zonal and meridional velocity components
    u_0 = 80.0f0
    lat_0 = convert(realT, π / 7)
    lat_1 = convert(realT, π / 2) - lat_0
    v_lon = galewsky_velocity(lat, u_0, lat_0, lat_1)
    v_lat = zero(eltype(x))

    # numerically integrate (here we use the QuadGK package) to get height
    galewsky_integral, _ = quadgk(latp -> galewsky_integrand(latp, u_0, lat_0, lat_1,
                                                             radius,
                                                             equations), -π / 2, lat)
    h = 10158.0f0 - radius / gravitational_acceleration * galewsky_integral

    # add perturbation to initiate instability
    α = convert(realT, 1 / 3)
    β = convert(realT, 1 / 15)
    lat_2 = convert(realT, π / 4)
    if (-π < lon) && (lon < π)
        h = h + 120.0f0 * cos(lat) * exp(-((lon / α)^2)) * exp(-((lat_2 - lat) / β)^2)
    end
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

initial_condition = initial_condition_barotropic_instability
source_terms = source_terms_convergence_test
tspan = (0.0, 6 * SECONDS_PER_DAY)

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
