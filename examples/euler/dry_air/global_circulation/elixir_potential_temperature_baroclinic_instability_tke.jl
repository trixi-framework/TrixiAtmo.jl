# The idealized baroclinic instability test case of
# `elixir_potential_temperature_baroclinic_instability.jl` with turbulent fluxes of a constant
# horizontal eddy viscosity and the turbulent kinetic energy (TKE) model of Deardorff (1980)
# in the vertical, with bulk fluxes of the momentum and the potential temperature through the
# surface as in the simple physics of DCMIP2016, integrated with a Rosenbrock-W method and a
# vertically implicit Jacobian (horizontally explicit, vertically implicit time integration),
# where the turbulent fluxes are treated explicitly and the source terms of the TKE are
# integrated exactly by Strang splitting. The positivity-preserving limiter of Zhang and Shu
# keeps the TKE nonnegative during the transport.
#
# References:
# - Paul A. Ullrich, Thomas Melvin, Christiane Jablonowski, Andrew Staniforth (2013)
#   A proposed baroclinic wave test case for deep- and shallow-atmosphere dynamical cores
#   https://doi.org/10.1002/qj.2241
# - James W. Deardorff (1980)
#   Stratocumulus-capped mixed layers derived from a three-dimensional model
#   https://doi.org/10.1007/BF00119502
# - Kevin A. Reed, Christiane Jablonowski (2012)
#   Idealized tropical cyclone simulations of intermediate complexity: A test case for AGCMs
#   https://doi.org/10.1029/2011MS000099
# - Paul A. Ullrich et al. (2016)
#   Dynamical Core Model Intercomparison Project (DCMIP2016) test case document
#   https://doi.org/10.5281/zenodo.1298671

using Trixi, TrixiAtmo

# Horizontal eddy viscosity of the momentum in m²/s, with the eddy diffusivity of the
# potential temperature for a turbulent Prandtl number of one
viscosity_horizontal = 1.0e5
turbulence = HorizontalVerticalTurbulence(ConstantEddyViscosity(viscosity_horizontal),
                                          DeardorffTKE())
equations = CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = 1004,
                                                                        c_v = 717,
                                                                        gravity = EARTH_GRAVITATIONAL_ACCELERATION,
                                                                        turbulence = turbulence)
equations_parabolic = TurbulentDiffusion3D(equations)

# Initial condition of the baroclinic instability with the constant specific TKE
# `initial_tke` in m²/s²
initial_tke = 0.01
initial_condition = let equations_dry_air = CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = 1004,
                                                                                                        c_v = 717,
                                                                                                        gravity = EARTH_GRAVITATIONAL_ACCELERATION),
    initial_tke = initial_tke

    function initial_condition_baroclinic_instability_tke(x, t, equations)
        u = initial_condition_baroclinic_instability(x, t, equations_dry_air)
        return SVector(u..., u[1] * initial_tke)
    end
end
aux_field = aux_field_gravity_earth

# The sea surface temperature is the temperature of the initial condition at the surface
@inline function sea_surface_temperature(x, equations)
    lon, lat, _ = cartesian_to_sphere(x)
    rho, _, p = TrixiAtmo.basic_state_baroclinic_instability_longitudinal_velocity(lon, lat,
                                                                                   zero(eltype(x)))
    return p / (rho * equations.R)
end

boundary_conditions = (; inside = boundary_condition_slip_wall,
                       outside = boundary_condition_slip_wall)
# Bulk fluxes through the surface and no turbulent fluxes through the top
boundary_condition_no_turbulent_flux = BoundaryConditionNeumann((x, t, equations) -> zero(SVector{6,
                                                                                                  eltype(x)}))
boundary_conditions_parabolic = (;
                                 inside = BoundaryConditionSurfaceFluxes(sea_surface_temperature),
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

# Nonnegative density of the TKE at all nodes of the stages and time steps
limiter! = PositivityPreservingLimiterZhangShu(thresholds = (0.0,),
                                               variables = (density_turbulent_kinetic_energy,))

dt = 85 # seconds
sol = Trixi.solve(ode, SSPRW32(); dt = dt, callback = callbacks, jacobian = jacobian,
                  stage_limiter! = limiter!, step_limiter! = limiter!,
                  source_splitting! = turbulence_sources_exact!)
