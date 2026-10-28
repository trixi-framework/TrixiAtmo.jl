# An idealized baroclinic instability test case for the compressible Euler equations in
# vector-invariant form with the entropy-stable and well-balanced fluxes of Artiano et al.
# (2026). For optimal results consider increasing the resolution to 16x16x8 trees per cube
# face.
#
# References:
# - Paul A. Ullrich, Thomas Melvin, Christiane Jablonowski, Andrew Staniforth (2013)
#   A proposed baroclinic wave test case for deep- and shallow-atmosphere dynamical cores
#   https://doi.org/10.1002/qj.2241
# - Marco Artiano, Kieran Ricardo, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2026)
#   Entropy-Stable and Well-Balanced Discontinuous Galerkin Methods for the Compressible
#   Euler Equations in Vector-Invariant Form
#   https://arxiv.org/abs/2610.05820

using OrdinaryDiffEqSSPRK
using Trixi, TrixiAtmo

equations = CompressibleEulerVectorInvariantEquationsWithGravity3D(c_p = 1004, c_v = 717,
                                                                   gravity = EARTH_GRAVITATIONAL_ACCELERATION)

initial_condition = initial_condition_baroclinic_instability

boundary_conditions = (; inside = boundary_condition_slip_wall,
                       outside = boundary_condition_slip_wall)

polydeg = 5
# Both fluxes return the conservative and the nonconservative contributions from a single
# kernel. `FluxTurbo` selects a SIMD-vectorized implementation of the flux differencing
# volume integral, which is equivalent to `volume_flux = flux_etec`.
surface_flux = flux_lmars_es
volume_flux = FluxTurbo(flux_etec)

solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))
trees_per_cube_face = (8, 4)
mesh = P4estMeshCubedSphereTopography(trees_per_cube_face..., EARTH_RADIUS, 30000,
                                      polydeg = polydeg, initial_refinement_level = 0)

# The Coriolis source term is the same as for the momentum, since the variables contain the
# velocity
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_coriolis,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.
T = 10 # 10 days
tspan = (0.0, T * SECONDS_PER_DAY) # time in seconds for 10 days

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 100, save_initial_solution = true,
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
