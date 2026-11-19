# References:
# - Christoph Schär, Daniel Leuenberger, Oliver Fuhrer, Daniel Lüthi, Claude Girard (2002)
#   A New Terrain-Following Vertical Coordinate Formulation for Atmospheric Prediction Models
#   Monthly Weather Review, 130(10), 2459–2480
#   https://doi.org/10.1175/1520-0493(2002)130

using OrdinaryDiffEqSSPRK
using Trixi, TrixiAtmo

###############################################################################
# semidiscretization of the compressible Euler equations
equations = CompressibleEulerPotentialTemperatureEquationsWithGravity2D(c_p = 1004,
                                                                        c_v = 717,
                                                                        gravity = EARTH_GRAVITATIONAL_ACCELERATION)

setup = MountainWaveSetup(theta_0 = 280, u0 = 10, brunt_vaisala_frequency = 0.01,
                          z_damping_start = 13000, z_damping_end = 21000,
                          x_damping_start = 20000, x_damping_end = 25000,
                          damping_rate = 0.03)
aux_field = aux_field_geopotential

boundary = BoundaryConditionDirichlet(setup)
boundary_conditions = (; x_neg = boundary,
                       x_pos = boundary,
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary)

polydeg = 3
basis = LobattoLegendreBasis(polydeg)
surface_flux = (FluxLMARS(340), flux_zero)
volume_flux = (flux_etec, flux_nonconservative_artiano_etal)
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

topography = SchärMountain(height = 250, half_width = 5000, wavelength = 4000)
L = 50000
H = 21000

cells_per_dimension = (100, 50)
cells_per_dimension = (20, 12)
mesh = P4estMesh(cells_per_dimension, polydeg = polydeg,
                 faces = terrain_following_faces(topography, -L / 2, L / 2, H),
                 initial_refinement_level = 0, periodicity = (false, false))

semi = SemidiscretizationHyperbolic(mesh, equations, setup, solver,
                                    source_terms = setup,
                                    boundary_conditions = boundary_conditions,
                                    aux_field = aux_field)

###############################################################################
# ODE solvers, callbacks etc.
T = 5
tspan = (0.0, T * 3600.0) #1 hour.

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
            maxiters = 1.0e7, ode_default_options()..., callback = callbacks)
