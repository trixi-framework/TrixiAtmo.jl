using OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiAtmo
using TrixiAtmo: initial_condition_dry_air_warm_bubble_generator, source_terms_gravity_cartZ_generator
using Plots


###############################################################################
# Parameters

RealType = Float64
n_dim = 2
n_passive = 0

parameters = Parameters{RealType}(;
    earth_gravitational_acceleration = 9.81
)


###############################################################################
# Thermodynamics

td_multi = IdealGasesAndLiquids(; parameters,
    n_gas = 2, n_condens = 1, n_precip = 1)

td_totE = TotalEnergy(td_multi)

# td_single = IdealGas(; parameters)
# td_potT = PotentialTemperature(td_single)

microphysics = MicrophysicsRelaxation{RealType}()


###############################################################################
# Equations

equations = CompressibleEulerAtmo{n_dim}(; parameters = parameters,
                                        thermodynamic_state = td_multi,
                                        thermodynamic_equation = td_totE,
                                        #thermodynamic_state = td_single,
                                        #thermodynamic_equation = td_potT,
                                        microphysics = microphysics)


###############################################################################
# Initial and boundary conditions

initial_condition_reference = initial_condition_dry_air_warm_bubble_generator(
    parameters;
    perturbation_center_x = 500.0,
    perturbation_center_z = 260.0,
    perturbation_radius = 250.0
)
initial_condition = transform_initial_condition(initial_condition_reference, equations)


###############################################################################
# Source terms

source_terms_gravity_reference = source_terms_gravity_cartZ_generator(Val(2);
    g = parameters.earth_gravitational_acceleration)
source_terms_gravity = transform_source_terms(source_terms_gravity_reference, equations)


###############################################################################
# Solver

surface_flux = flux_lax_friedrichs
#surface_flux = FluxLMARS(340.0)
#volume_flux = flux_ranocha
#volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

polydeg = 3
solver = DGSEM(polydeg, surface_flux) #, volume_integral)

coordinates_min = (0.0, 0.0)
#coordinates_max = (20_000.0, 10_000.0)
coordinates_max = (1000.0, 1500.0)

trees_per_dimension = (8, 8)
mesh = P4estMesh(trees_per_dimension, polydeg = 3,
                 coordinates_min = coordinates_min, coordinates_max = coordinates_max,
                 periodicity = (true, false), initial_refinement_level = 2)


boundary_conditions = Dict(:y_neg => boundary_condition_slip_wall,
                           :y_pos => boundary_condition_slip_wall)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_gravity,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 10.0)  # 1000 seconds final time

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000

analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = analysis_interval,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = "out_warm_bubble")

visualization = VisualizationCallback(semi; interval = 10, show_mesh = false,
                                      solution_variables = cons2cons,
                                      variable_names = ["rho_v1", "rho_v2"])

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        visualization,
                        stepsize_callback)

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            dt = 0.1, # solve needs some value here but it will be overwritten by the stepsize_callback
            ode_default_options()..., callback = callbacks);
