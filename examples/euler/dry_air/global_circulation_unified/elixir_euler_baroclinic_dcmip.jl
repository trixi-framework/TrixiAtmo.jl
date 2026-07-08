using OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiAtmo: initial_condition_dry_air_warm_bubble_generator,
                 source_terms_gravity_cartZ_generator

###############################################################################
# Parameters

RealType = Float64
earth_scale = 1.0
earth_radius = 6.371229e6

parameters = Parameters{RealType}(;
                                  earth_gravitational_acceleration = 9.81,
                                  c_dry_air_const_pressure = 1004,
                                  c_dry_air_const_volume = 717,
                                  earth_rotation_rate = 7.29212e-5 * earth_scale,
                                  earth_radius = earth_radius)

###############################################################################
# Thermodynamics

td_single = IdealGas(; parameters)
td_potT = EnergyTotal(td_single)

###############################################################################
# Equations

equations = CompressibleEulerAtmo(; n_dims = 3, n_vars_aux = 1,
                                  parameters = parameters,
                                  thermodynamic_state = td_single,
                                  thermodynamic_equation = td_potT)

###############################################################################
# Initial and boundary conditions

initial_condition_reference = initial_condition_baroclinic_instability_generator(parameters;
                                                                                 perturbation_stream_function = False(),
                                                                                 earth_scaling_factor = earth_scale)

initial_condition = transform_initial_condition(initial_condition_reference, equations)

boundary_conditions = (; inside = boundary_condition_slip_wall,
                       outside = boundary_condition_slip_wall)

###############################################################################
# Source terms

source_terms_gravity_reference = source_terms_gravity_spherical_generator(equations)

source_terms_coriolis_reference = source_terms_coriolis_generator(equations)

source_terms = transform_source_terms_sum((source_terms_gravity_reference,
                                           source_terms_coriolis_reference),
                                          equations)

###############################################################################
# Solver

polydeg = 5

surface_flux = (FluxLMARS(340), flux_zero)
volume_flux = (flux_ranocha, flux_nonconservative_waruszewski_etal)

solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Mesh

trees_per_cube_face = (8, 4)
mesh = P4estMeshCubedSphere(trees_per_cube_face..., earth_radius, 30000,
                            polydeg = polydeg, initial_refinement_level = 0)

###############################################################################
# Semidiscretization

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms,
                                    boundary_conditions = boundary_conditions)
T = 10 # 10 days
tspan = (0.0, T * SECONDS_PER_DAY)

ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 1000

analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = analysis_interval,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = "out_baroclinic")

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        stepsize_callback)

tol = 1e-6
sol = solve(ode,
            SSPRK43(thread = Trixi.Threaded());
            abstol = tol, reltol = tol, ode_default_options()...,
            callback = callbacks)
