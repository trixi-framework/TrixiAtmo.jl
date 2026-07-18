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

td_variants = ["Etot", "Tpot", "Eint"]

td_variant = 3
amr = true
tracer = true
vel_max = 45
tracer_max = 0.1

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

if td_variant == 1
    td_eq = EnergyTotal(td_single)
elseif td_variant == 2
    td_eq = PotentialTemperature(td_single)
elseif td_variant == 3
    td_eq = EnergyInternal(td_single)
end

###############################################################################
# Equations

equations = CompressibleEulerAtmo(; n_dims = 3, n_vars_aux = 1,
                                  n_vars_passive = tracer ? 1 : 0,
                                  parameters = parameters,
                                  thermodynamic_state = td_single,
                                  thermodynamic_equation = td_eq)

###############################################################################
# Initial and boundary conditions

initial_condition_reference = initial_condition_baroclinic_instability_generator(parameters;
                                                                                 perturbation_stream_function = true,
                                                                                 earth_scaling_factor = earth_scale)

if tracer
    @inline function initial_tracer(x, equations::CompressibleEulerAtmo)
        lon, lat, r = cartesian_to_spherical_coordinates(x)
        z = r - equations.parameters.earth_radius

        lon_pos = 0.75
        lon_dist = min(abs(lon - lon_pos), abs(lon + 2 * pi - lon_pos),
                   abs(lon - 2 * pi - lon_pos))

        # Initial condition for tracers: blob as a fraction of density
        return 1e-2 +
               0.1 * exp(-45 * (((lon_dist))^2 +
                    ((lat - 0.45) / 1.2)^2 +
                    ((z - 10_000) / 20_000)^2))
    end
else
    initial_tracer = (x, e) -> 0
end

initial_condition = transform_initial_condition(initial_condition_reference, equations;
                                                initial_condition_passive = initial_tracer)

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

if td_variant == 1
    surface_flux = (FluxLMARS(340), flux_zero)
    volume_flux = (flux_ranocha, flux_nonconservative_waruszewski_etal)
elseif td_variant == 2
    surface_flux = (FluxLMARS(340), flux_zero)
    volume_flux = (flux_tec, flux_nonconservative_souza_etal)
elseif td_variant == 3
    surface_flux = (flux_surface_artiano, flux_nonconservative_surface_artiano)
    volume_flux = (flux_volume_artiano, flux_nonconservative_waruszewski_etal)
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

days = 14
tspan = (0.0, days * SECONDS_PER_DAY)

ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 10_000

analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

dir_name = "out_baroclinic_$(td_variants[td_variant])"
if tracer
    dir_name *= "+t"
end
if amr
    dir_name *= "_velmax-$vel_max"
end

save_solution = SaveSolutionCallback(interval = analysis_interval,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim,
                                     output_directory = dir_name)

stepsize_callback = StepsizeCallback(cfl = 0.9)

callbacks = CallbackSet(summary_callback,
                        stepsize_callback,
                        analysis_callback,
                        alive_callback)
if amr
    @inline function velocity_tracer(u, aux, equations::CompressibleEulerAtmo)
        vel = velocity_magnitude(u, equations) / vel_max
        rho = TrixiAtmo.density_total(u, equations)
        tracer = TrixiAtmo.vars_passive(u, equations)[1] / (rho * tracer_max)
        return max(vel, tracer)
    end
    amr_indicator = IndicatorMax(semi, variable = velocity_tracer)

    amr_controller = ControllerThreeLevel(semi, amr_indicator,
                                          base_level = 0,
                                          med_level = 1, med_threshold = 0.65,
                                          # 27.8 initially -> > 27.8/max
                                          max_level = 2, max_threshold = 1.0)

    amr_callback = AMRCallback(semi, amr_controller,
                               interval = analysis_interval,
                               adapt_initial_condition = true,
                               adapt_initial_condition_only_refine = true)

    callbacks = CallbackSet(callbacks.discrete_callbacks..., amr_callback)
end

callbacks = CallbackSet(callbacks.discrete_callbacks..., save_solution)

tol = 1e-6
sol = solve(ode,
            CarpenterKennedy2N54(williamson_condition = false);
            #SSPRK43(); # thread = Trixi.Threaded());
            maxiters = 1e8,
            abstol = tol, reltol = tol, ode_default_options()...,
            callback = callbacks)
