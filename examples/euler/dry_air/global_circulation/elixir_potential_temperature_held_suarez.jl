# Held-Suarez test case
#
# References:
# - Souza et al. (2023):
#   The Flux-Differencing Discontinuous Galerkin Method Applied to an Idealized Fully
#   Compressible Nonhydrostatic Dry Atmosphere
#   https://doi.org/10.1029/2022MS003527

using OrdinaryDiffEqLowStorageRK
using Trixi, TrixiAtmo
using LinearAlgebra: norm

function initial_condition_isothermal(x, t,
                                      equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    # equation (60) in the paper
    T = 285

    @unpack p_0, R = equations

    r = norm(x)
    # Make sure that r is not smaller than radius_earth
    z = max(r - EARTH_RADIUS, 0.0)
    r = z + EARTH_RADIUS

    # pressure, geopotential formulation
    p = p_0 *
        exp(EARTH_GRAVITATIONAL_ACCELERATION *
            (EARTH_RADIUS^2 / r - EARTH_RADIUS) /
            (R * T))

    # density (via ideal gas law)
    rho = p / (R * T)

    # geopotential
    phi = EARTH_GRAVITATIONAL_ACCELERATION * (EARTH_RADIUS - EARTH_RADIUS^2 / r)

    return prim2cons(SVector(rho, 0, 0, 0, p, phi), equations)
end

@inline function source_terms_coriolis(u, x, t,
                                       equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    du0 = zero(eltype(u))

    # Coriolis term, -2Ω × ρv = -2 * angular_velocity * (0, 0, 1) × u[2:4]
    du2 = 2 * EARTH_ROTATION_RATE * u[3]
    du3 = -2 * EARTH_ROTATION_RATE * u[2]

    return SVector(du0, du2, du3, du0, du0, du0)
end

function cartesian_to_sphere(x)
    r = norm(x)
    lambda = atan(x[2], x[1])
    if lambda < 0
        lambda += 2 * pi
    end
    phi = asin(x[3] / r)

    return lambda, phi, r
end

@inline function source_terms_hs_relaxation(u, x, t,
                                            equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    # equations (55)-(58) in the paper
    k_f = 1 / SECONDS_PER_DAY         # Damping scale for momentum
    k_a = 1 / (40 * SECONDS_PER_DAY)  # Polar relaxation scale
    k_s = 1 / (4 * SECONDS_PER_DAY)   # Equatorial relaxation scale
    T_min = 200                       # Minimum equilibrium temperature
    T_equator = 315                   # Equatorial equilibrium temperature
    deltaT = 60                       # Latitudinal temperature difference
    deltaTheta = 10                   # Vertical temperature difference
    sigma_b = 0.7                     # Dimensionless damping height

    @unpack p_0, c_p, c_v, R = equations

    p = pressure(u, equations)
    lon, lat, r = cartesian_to_sphere(x)
    temperature = p / (u[1] * R)

    sigma = p / p_0   # "p_0 instead of instantaneous surface pressure"
    delta_sigma = max(0, (sigma - sigma_b) / (1 - sigma_b))   # "height factor"
    k_v = k_f * delta_sigma
    k_T = k_a + (k_s - k_a) * delta_sigma * cos(lat)^4

    T_equi = max(T_min,
                 (T_equator - deltaT * sin(lat)^2 - deltaTheta * log(sigma) * cos(lat)^2) *
                 sigma^(R / c_p))

    # project onto r
    # Make sure that r is not smaller than radius_earth
    z = max(r - EARTH_RADIUS, 0.0)
    r = z + EARTH_RADIUS
    dotprod = (u[2] * x[1] + u[3] * x[2] + u[4] * x[3]) / (r * r)

    du0 = zero(eltype(u))

    du2 = -k_v * (u[2] - dotprod * x[1])
    du3 = -k_v * (u[3] - dotprod * x[2])
    du4 = -k_v * (u[4] - dotprod * x[3])

    du5 = -k_T * u[5] / temperature * (temperature - T_equi)

    return SVector(du0, du2, du3, du4, du5, du0)
end

@inline function source_terms_held_suarez(u, x, t,
                                          equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    return source_terms_coriolis(u, x, t, equations) +
           source_terms_hs_relaxation(u, x, t, equations)
end

equations = CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = 1004,
                                                                        c_v = 717,
                                                                        gravity = EARTH_GRAVITATIONAL_ACCELERATION)

initial_condition = initial_condition_isothermal

boundary_conditions = Dict(:inside => boundary_condition_slip_wall,
                           :outside => boundary_condition_slip_wall)

polydeg = 4
surface_flux = (FluxLMARS(340), flux_zero)
volume_flux = (flux_ec, flux_nonconservative_waruszewski_etal)

solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

# these are the settings used in Souza et al.
lat_lon_trees_per_dim = 10
layers = 8

mesh = Trixi.T8codeMeshCubedSphere(lat_lon_trees_per_dim, layers, EARTH_RADIUS, 30000.0,
                                   polydeg = polydeg, initial_refinement_level = 0)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_held_suarez,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

# this test case is typically run for a long time and evaluated statistically
T = 365
tspan = (0.0, T * SECONDS_PER_DAY)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 5000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 5000,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution)

###############################################################################
# Use a Runge-Kutta method with automatic (error based) time step size control
# Enable threading of the RK method for better performance on multiple threads
sol = solve(ode,
            RDPK3SpFSAL49(thread = Trixi.True());
            abstol = 1.0e-5, reltol = 1.0e-5, ode_default_options()...,
            callback = callbacks);
