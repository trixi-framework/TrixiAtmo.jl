# The idealized baroclinic instability test case of
# `elixir_potential_temperature_baroclinic_instability.jl` with turbulent fluxes of constant
# horizontal and vertical eddy viscosities, integrated with a Rosenbrock-W method and a
# vertically implicit Jacobian (horizontally explicit, vertically implicit time integration),
# where the turbulent fluxes are treated explicitly
#
# References:
# - Paul A. Ullrich, Thomas Melvin, Christiane Jablonowski, Andrew Staniforth (2013)
#   A proposed baroclinic wave test case for deep- and shallow-atmosphere dynamical cores
#   https://doi.org/10.1002/qj.2241

using Trixi, TrixiAtmo

# Eddy viscosities of the momentum in m²/s, with eddy diffusivities of the potential
# temperature for a turbulent Prandtl number of one
viscosity_horizontal = 1.0e5
viscosity_vertical = 1.0
turbulence = HorizontalVerticalTurbulence(ConstantEddyViscosity(viscosity_horizontal),
                                          ConstantEddyViscosity(viscosity_vertical))
equations = CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = 1004,
                                                                        c_v = 717,
                                                                        gravity = EARTH_GRAVITATIONAL_ACCELERATION,
                                                                        turbulence = turbulence)
equations_parabolic = TurbulentDiffusion3D(equations)

initial_condition = initial_condition_baroclinic_instability
aux_field = aux_field_gravity_earth

boundary_conditions = (; inside = boundary_condition_slip_wall,
                       outside = boundary_condition_slip_wall)
# No turbulent fluxes through the surface and the top
boundary_condition_no_turbulent_flux = BoundaryConditionNeumann((x, t, equations) -> zero(SVector{5,
                                                                                                  eltype(x)}))
boundary_conditions_parabolic = (; inside = boundary_condition_no_turbulent_flux,
                                 outside = boundary_condition_no_turbulent_flux)

polydeg = 4
surface_flux = flux_lmars_es
volume_flux = FluxTurbo(flux_ec, flux_nonconservative_isothermal)

solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))
trees_per_cube_face = (20, 13)
mesh = P4estMeshCubedSphereTopography(trees_per_cube_face..., EARTH_RADIUS, 30000,
                                      polydeg = polydeg, initial_refinement_level = 0,
                                      keep_columns_together = true)

semi = SemidiscretizationHyperbolicParabolic(mesh, (equations, equations_parabolic),
                                             initial_condition, solver,
                                             source_terms = source_terms_coriolis,
                                             boundary_conditions = (boundary_conditions,
                                                                    boundary_conditions_parabolic),
                                             aux_field = aux_field)

###############################################################################
# ODE solvers, callbacks etc.
T = 10 # 10 days
tspan = (0.0, T * SECONDS_PER_DAY) # time in seconds for 10 days

# Setting `real_type` allows to change the real number type, e.g., to `Float32`.
# This is particularly useful when changing the `storage_type` to a GPU array
# type such as `ROCArray` (AMD) or `CuArray` (NVIDIA CUDA).
# On such backends, `tuning` selects performance settings of the kernels, e.g.,
# `Tuning(flux_differencing_kernel = FullSweep())` for the kernel used for the
# flux differencing volume integral.
semi = trixi_adapt(semi; real_type = nothing, storage_type = nothing,
                   tuning = nothing)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = SECONDS_PER_DAY, save_initial_solution = true,
                                     save_final_solution = true)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution)

###############################################################################
# Use a Rosenbrock-W method with the vertically implicit Jacobian, such that the time step
# is restricted by the horizontal dynamics only
jacobian = JacobianVerticalTwoSided(c = 340)

dt = 85 # seconds
sol = Trixi.solve(ode, SSPRW32(); dt = dt, callback = callbacks, jacobian = jacobian)
