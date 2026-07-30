# References:
# - Dale R. Durran, Joseph B. Klemp (1983)
#   A Compressible Model for the Simulation of Moist Mountain Waves
#   Monthly Weather Review, 111(12), 2341–2361
#   https://doi.org/10.1175/1520-0493(1983)111
#
# - F.X. Giraldo, M. Restelli (2008)
#   A study of spectral element and discontinuous Galerkin methods for the Navier–Stokes equations in nonhydrostatic mesoscale atmospheric modeling: Equation sets and test cases
#   Journal of Computational Physics, 227(8), 3849–3877
#   https://doi.org/10.1016/j.jcp.2007.12.009

using OrdinaryDiffEqSSPRK
using Trixi, TrixiAtmo

equations = CompressibleEulerPotentialTemperatureEquationsWithGravity2D(c_p = 1004,
                                                                        c_v = 717,
                                                                        gravity = EARTH_GRAVITATIONAL_ACCELERATION)

# The background state is isothermal at T_0, which is the particular case of a constant
# Brunt-Väisälä frequency obtained with theta_0 = T_0 and N = g / sqrt(c_p * T_0).
T_0 = 250
setup = MountainWaveSetup(theta_0 = T_0, u0 = 20,
                          brunt_vaisala_frequency = equations.g /
                                                    sqrt(equations.c_p * T_0),
                          z_damping_start = 15000, z_damping_end = 30000,
                          x_damping_start = 60000, x_damping_end = 120000,
                          damping_rate = 0.035)

boundary = BoundaryConditionDirichlet(setup)

boundary_conditions = (; x_neg = boundary,
                       x_pos = boundary,
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary)

# We have an isothermal background state with T0 = 250 K.
# The reference speed of sound can be computed as:
# cs = sqrt(gamma * R * T0)
cs = sqrt(equations.gamma * equations.R * T_0)
polydeg = 3
basis = LobattoLegendreBasis(polydeg)

surface_flux = (FluxLMARS(cs), flux_zero)
volume_flux = (flux_tec, flux_nonconservative_waruszewski_etal)
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

# The half width of 10 km gives a * N / u0 >> 1, i.e., a hydrostatic flow regime.
topography = WitchOfAgnesi(height = 1, half_width = 10000)
L = 240000
H = 30000

cells_per_dimension = (100, 60)
mesh = P4estMesh(cells_per_dimension, polydeg = polydeg,
                 faces = terrain_following_faces(topography, -L / 2, L / 2, H),
                 initial_refinement_level = 0, periodicity = (false, false))

semi = SemidiscretizationHyperbolic(mesh, equations, setup, solver,
                                    source_terms = setup,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.
T = 12.5
tspan = (0.0, T * 3600.0)  # 12.5 hours final time
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000

analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback)

###############################################################################
# run the simulation
sol = solve(ode,
            SSPRK43(thread = Trixi.Threaded());
            maxiters = 1.0e7,
            ode_default_options()..., callback = callbacks)
