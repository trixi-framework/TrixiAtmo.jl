# An idealized baroclinic instability test case, integrated with a Rosenbrock-W method and a
# vertically implicit Jacobian (horizontally explicit, vertically implicit time integration)
#
# This elixir takes about 2 minutes on an NVIDIA GH200 GPU.
#
# References:
# - Paul A. Ullrich, Thomas Melvin, Christiane Jablonowski, Andrew Staniforth (2013)
#   A proposed baroclinic wave test case for deep- and shallow-atmosphere dynamical cores
#   https://doi.org/10.1002/qj.2241

using Trixi, TrixiAtmo

equations = CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = 1004,
                                                                        c_v = 717,
                                                                        gravity = EARTH_GRAVITATIONAL_ACCELERATION)

initial_condition = initial_condition_baroclinic_instability

boundary_conditions = (; inside = boundary_condition_slip_wall,
                       outside = boundary_condition_slip_wall)

polydeg = 4
surface_flux = flux_es
volume_flux = FluxTurbo(flux_ec, flux_nonconservative_waruszewski_etal)

solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))
trees_per_cube_face = (20, 13)
mesh = P4estMeshCubedSphereTopography(trees_per_cube_face..., EARTH_RADIUS, 30000,
                                      polydeg = polydeg, initial_refinement_level = 0,
                                      keep_columns_together = true)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_coriolis,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.
T = 10 # 10 days
tspan = (0.0, T * SECONDS_PER_DAY) # time in seconds for 10 days

# Setting `real_type` allows to change the real number type, e.g., to `Float32`.
# This is particularly useful when changing the `storage_type` to a GPU array
# type such as `ROCArray` (AMD) or `CuArray` (NVIDIA CUDA).
# On such backends, `flux_differencing_kernel` selects the kernel used for the
# flux differencing volume integral, one of `HalfSweep()` (default), `FullSweep()`,
# or `FullSweepGlobal()`.
semi = trixi_adapt(semi; real_type = nothing, storage_type = nothing,
                   flux_differencing_kernel = nothing)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback)

###############################################################################
# Use a Rosenbrock-W method with the vertically implicit Jacobian, such that the time step
# is restricted by the horizontal dynamics only
jacobian = VerticalJacobianTwoSided(semi; c = 340)

dt = 85.0 # seconds
sol = Trixi.solve(ode, SSPRW32(); dt = dt, callback = callbacks, jacobian = jacobian)
