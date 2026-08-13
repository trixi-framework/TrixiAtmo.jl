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

using OrdinaryDiffEqSSPRK
using Trixi, TrixiAtmo

equations = CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = 1004,
                                                                        c_v = 717,
                                                                        gravity = EARTH_GRAVITATIONAL_ACCELERATION)

initial_condition = initial_condition_vortex_shedding

boundary_conditions = (; inside = boundary_condition_slip_wall,
                       outside = boundary_condition_slip_wall)

polydeg = 5
surface_flux = (FluxLMARS(340), flux_zero)
volume_flux = (flux_tec, flux_nonconservative_souza_etal)

solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

function initial_topography_gaussian(lat, lon)
    h_0 = 2000
    h0 = 10000
    d = 12500
    lambda_c = pi
    phi_c = 20 * pi / 180

    r_earth = SMALL_EARTH_RADIUS
    dist = r_earth * acos(clamp(sin(phi_c) * sin(lat) +
                      cos(phi_c) * cos(lat) * cos(lon - lambda_c), -1, 1))

    return h_0 * exp(-(dist / d)^2)
end

trees_per_cube_face = (16, 8)

mesh = P4estMeshCubedSphereTopography(trees_per_cube_face..., SMALL_EARTH_RADIUS,
                                      20000,
                                      polydeg = polydeg,
                                      initial_refinement_level = 0,
                                      initial_topography = initial_topography_gaussian,
                                      adapt_vertical_grid = GalChen(),
                                      keep_columns_together = true)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_coriolis_small_earth,
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
            SSPRK43(thread = Trixi.Threaded());
            abstol = tol, reltol = tol, ode_default_options()...,
            callback = callbacks)
