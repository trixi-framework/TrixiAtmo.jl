using OrdinaryDiffEqSSPRK
using Trixi
using TrixiAtmo: IdealGas,
                 initial_condition_baroclinic_instability_generator,
                 source_terms_gravity_spherical_generator,
                 source_terms_coriolis_generator,
                 geopotential_spherical,
                 velocity_magnitude

###############################################################################
# Parameters

RealType = Float64
earth_scale = 1.0
earth_radius = 6.371229e6

# 1: EnergyTotal
# 2: PotentialTemperature
td_variant = 2
amr = true

parameters = Parameters{RealType}(;
                                  earth_gravitational_acceleration = 9.81,
                                  c_dry_air_const_pressure = 1004,
                                  c_dry_air_const_volume = 717,
                                  earth_rotation_rate = 7.29212e-5 * earth_scale,
                                  earth_radius = earth_radius)

###############################################################################
# Thermodynamics

td_single = IdealGas(; parameters)
if td_variant == 1
    td_eq = EnergyTotal(td_single)
else
    td_eq = PotentialTemperature(td_single)
end

###############################################################################
# Equations

equations = CompressibleEulerAtmo(; n_dims = 3, n_vars_aux = 1,
                                  parameters = parameters,
                                  thermodynamic_state = td_single,
                                  thermodynamic_equation = td_eq)

###############################################################################
# Initial and boundary conditions

initial_condition_reference = initial_condition_baroclinic_instability_generator(parameters;
                                                                                 perturbation_stream_function = true,
                                                                                 earth_scaling_factor = earth_scale)

initial_condition = transform_initial_condition(initial_condition_reference, equations)

# TODO simple!
boundary_conditions = (; inside = boundary_condition_slip_wall_simple,
                       outside = boundary_condition_slip_wall_simple)

###############################################################################
# Source terms

# NOTE: nonconservative discretizations of gravity requires aux vars with geopotential
#       gravity is then not part of source terms!

# source_terms_gravity_reference = source_terms_gravity_spherical_generator(equations)

source_terms_coriolis_reference = source_terms_coriolis_generator(equations)

source_terms = transform_source_terms(source_terms_coriolis_reference, equations)

###############################################################################
# Solver

polydeg = 5

surface_flux = (FluxLMARS(340), flux_zero)
if td_variant == 1
    volume_flux = (flux_ranocha, flux_nonconservative_waruszewski_etal)
else
    volume_flux = (flux_tec, flux_nonconservative_souza_etal)
end

solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Mesh

trees_per_cube_face = (8, 4)
initial_refinement_level = 0
if !amr
    initial_refinement_level = 1
end

mesh = P4estMeshCubedSphere(trees_per_cube_face..., earth_radius, 30000,
                            polydeg = polydeg,
                            initial_refinement_level = initial_refinement_level)

###############################################################################
# Semidiscretization

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms,
                                    boundary_conditions = boundary_conditions,
                                    aux_field = geopotential_spherical)

days = 10
tspan = (0.0, days * SECONDS_PER_DAY)

ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 10_000

analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = analysis_interval,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim,
                                     output_directory = "out_baroclinic")

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution)
if amr
    amr_indicator = IndicatorMax(semi, variable = velocity_magnitude)

    amr_controller = ControllerThreeLevel(semi, amr_indicator,
                                          base_level = 0,
                                          med_level = 1, med_threshold = 28.0,
                                          max_level = 2, max_threshold = 45.0)

    amr_callback = AMRCallback(semi, amr_controller,
                               interval = analysis_interval,
                               adapt_initial_condition = true,
                               adapt_initial_condition_only_refine = true)

    callbacks = CallbackSet(callbacks.discrete_callbacks..., amr_callback)
end

tol = 1e-6
sol = solve(ode,
            SSPRK43(thread = Trixi.Threaded());
            maxiters = 1e7,
            abstol = tol, reltol = tol, ode_default_options()...,
            callback = callbacks)
