# References:
# - Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
#   Structure-Preserving High-Order Methods for the Compressible Euler Equations
#   in Potential Temperature Formulation for Atmospheric Flows
#   https://arxiv.org/abs/2509.10311 (pre-print)
using OrdinaryDiffEqSSPRK
using Trixi, TrixiAtmo

equations = CompressibleEulerPotentialTemperatureEquationsWithGravity1D(c_p = 1004,
                                                                        c_v = 717,
                                                                        gravity = EARTH_GRAVITATIONAL_ACCELERATION)
polydeg = 2
basis = LobattoLegendreBasis(polydeg)
cs = 340.0
surface_flux = (FluxLMARS(cs), flux_zero)
volume_flux = (flux_tec, flux_nonconservative_artiano_etal)

volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
solver = DGSEM(basis, surface_flux, volume_integral)

coordinates_min = (0.0,)
coordinates_max = (1000.0,)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 3,
                periodicity = false)

boundary_conditions = (; x_pos = boundary_condition_slip_wall,
                       x_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall,
                       y_neg = boundary_condition_slip_wall)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_adiabatic, solver,
                                    boundary_conditions = boundary_conditions)

dt = 0.01
tspan = (0, 0.4)

summary_callback = SummaryCallback()

analysis_interval = 10
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback)

ode = semidiscretize(semi, tspan)

sol = solve(ode, SSPRK43(); dt = dt, ode_default_options()..., callback = callbacks,
            adaptive = false)
