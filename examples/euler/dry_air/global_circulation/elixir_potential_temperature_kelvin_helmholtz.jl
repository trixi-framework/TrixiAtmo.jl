# An idealized mountain-triggered mesoscale flow test case (DCMIP-2025 Test Case 2).
# This setup examines the interaction between a stratified rotating flow and isolated
# topography, leading to vortex shedding in the lee of the mountain.
#
# The test-case uses a reduced-radius Earth configuration (R = Rₑ/20), scaled rotation rate,
# and a Gaussian mountain with peak height 2000 m.
#
# References:
# - DCMIP-2025 Organizing Committee (2025)
#   Test Case 2: Mountain-Triggered Meso-scale Flow Phenomena
#   https://sites.google.com/umich.edu/dcmip-2025/dcmip-2025-test-cases/test-case-2-mountain-triggered-meso-scale-flow-phenomena

using OrdinaryDiffEqLowStorageRK
using Trixi, TrixiAtmo
using LinearAlgebra: norm

function cartesian_to_sphere(x)
    r = norm(x)
    lambda = atan(x[2], x[1])
    if lambda < 0
        lambda += 2 * pi
    end
    phi = asin(x[3] / r)

    return lambda, phi, r
end

function initial_condition_kelvin_helmholtz_sphere(
    x, t,
    equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D
)
    lon, lat, r = cartesian_to_sphere(x)

    RealT = eltype(x)

    radius_earth = EARTH_RADIUS / 20
    R = equations.c_p - equations.c_v

    slope = 15
    lat_width = 0.5

    B = tanh(slope * (lat + lat_width)) -
        tanh(slope * (lat - lat_width))

    rho = 0.5 + 0.75 * B

    u_lon = 0.5 * (B - 1)

    perturb = convert(RealT, 0.1)
    mode = 2
    u_lat = perturb * sin(mode * lon)

    u_rad = zero(RealT)

    e_lon = SVector(-sin(lon), cos(lon), zero(RealT))
    e_lat = SVector(
        -sin(lat) * cos(lon),
        -sin(lat) * sin(lon),
         cos(lat)
    )
    e_rad = SVector(
        cos(lat) * cos(lon),
        cos(lat) * sin(lon),
        sin(lat)
    )

    v = u_lon * e_lon + u_lat * e_lat + u_rad * e_rad

    g = EARTH_GRAVITATIONAL_ACCELERATION
    z = max(norm(x) - radius_earth, zero(RealT))
    phi = g * z

    p = 1.0

    return prim2cons(SVector(rho, v[1], v[2], v[3], p, phi), equations)
end

@inline function source_terms_coriolis(u, x, t,
                                       equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    radius_earth = EARTH_RADIUS / 20  # a
    angular_velocity = EARTH_ROTATION_RATE * 20  # Ω

    du1 = zero(eltype(u))
    du2 = zero(eltype(u))
    du3 = zero(eltype(u))
    du4 = zero(eltype(u))
    du5 = zero(eltype(u))
    # Coriolis term, -2Ω × ρv = -2 * angular_velocity * (0, 0, 1) × u[2:4]
    du2 -= -2 * angular_velocity * u[3]
    du3 -= 2 * angular_velocity * u[2]

    return SVector(du1, du2, du3, du4, du5, zero(eltype(u)))
end
equations = CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = 1004,
                                                                        c_v = 717,
                                                                        gravity = EARTH_GRAVITATIONAL_ACCELERATION)

initial_condition = initial_condition_kelvin_helmholtz_sphere

boundary_conditions = (; inside = boundary_condition_slip_wall,
                       outside = boundary_condition_slip_wall)

polydeg = 5
surface_flux = (FluxLMARS(340), flux_zero)
volume_flux = (flux_tec, flux_nonconservative_souza_etal)

solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

trees_per_cube_face = (8, 4)

mesh = P4estMeshCubedSphere(trees_per_cube_face..., EARTH_RADIUS / 20,
                                      20000,
                                      polydeg = polydeg,
                                      initial_refinement_level = 0)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_coriolis,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.
T = 1.0 # 20 small earth days
tspan = (0.0, T * SECONDS_PER_DAY) # time in seconds for 1 standard earth day

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 5760, save_initial_solution = true,
                                     save_final_solution = true)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution)

###############################################################################
# Use a Runge-Kutta method with automatic (error based) time step size control
# Enable threading of the RK method for better performance on multiple threads
tol = 1e-6
sol = solve(ode,

            abstol = tol, reltol = tol, ode_default_options()...,
            callback = callbacks)
