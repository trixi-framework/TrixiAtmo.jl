using OrdinaryDiffEqSSPRK
using OrdinaryDiffEqLowStorageRK

using Trixi
using TrixiAtmo: IdealGas,
                 initial_condition_dry_air_warm_bubble_generator,
                 source_terms_gravity_cartZ_generator
using Plots

###############################################################################
# Parameters

RealType = Float64
n_dim = 2

# override to match parameters in Trixi.jl
parameters = Parameters{RealType}(;
                                  earth_gravitational_acceleration = EARTH_GRAVITATIONAL_ACCELERATION,
                                  c_dry_air_const_pressure = 1004.0,
                                  c_dry_air_const_volume = 717.0)

###############################################################################
# Thermodynamics

td_single = IdealGas(; parameters)
td_potT = PotentialTemperature(td_single)

###############################################################################
# Equations

equations = CompressibleEulerAtmo(; n_dims = 2,
                                  parameters = parameters,
                                  thermodynamic_state = td_single,
                                  thermodynamic_equation = td_potT)

###############################################################################
# Initial and boundary conditions

initial_condition_reference = initial_condition_dry_air_warm_bubble_generator(parameters;
                                                                              perturbation_center_x = 500.0,
                                                                              perturbation_center_z = 260.0,
                                                                              perturbation_radius = 250.0)
initial_condition = transform_initial_condition(initial_condition_reference, equations)

boundary_conditions = (; y_neg = boundary_condition_slip_wall_simple,
                       y_pos = boundary_condition_slip_wall_simple)

###############################################################################
# Source terms

source_terms_gravity_reference = source_terms_gravity_cartZ_generator(equations)
source_terms_gravity = transform_source_terms(source_terms_gravity_reference, equations)

###############################################################################
# Solver

surface_flux = FluxLMARS(340)
volume_flux = flux_tec
polydeg = 3
solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Mesh

coordinates_min = (0.0, 0.0)
coordinates_max = (1000.0, 1500.0)
cells_per_dimension = (64, 64)
mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                      periodicity = (true, false))

###############################################################################
# Semidiscretization
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_gravity,
                                    boundary_conditions = boundary_conditions)
tspan = (0.0, 10.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)
#extra_analysis_integrals = (entropy,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 100,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = "out_bubble_potT_new")

visualization = VisualizationCallback(semi; interval = 100, show_mesh = false,
                                      solution_variables = cons2cons,
                                      variable_names = ["rho_v1", "rho_v2"])

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        visualization, save_solution,
                        alive_callback)

sol = solve(ode,
            SSPRK43(thread = Trixi.Threaded());
            maxiters = 1.0e7, ode_default_options()..., callback = callbacks)
