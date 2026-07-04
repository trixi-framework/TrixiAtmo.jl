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
                                  earth_gravitational_acceleration = 9.81,
                                  c_dry_air_const_pressure = 1004.0,
                                  c_dry_air_const_volume = 717.0)

###############################################################################
# Thermodynamics

td_single = IdealGas(; parameters)
td_totE = TotalEnergy(td_single)

###############################################################################
# Equations

equations = CompressibleEulerAtmo(; n_dims = 2,
                                  parameters = parameters,
                                  thermodynamic_state = td_single,
                                  thermodynamic_equation = td_totE)

###############################################################################
# Initial and boundary conditions

initial_condition_reference = initial_condition_dry_air_warm_bubble_generator(parameters)
initial_condition = transform_initial_condition(initial_condition_reference, equations)

boundary_conditions = (; y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

###############################################################################
# Source terms

source_terms_gravity_reference = source_terms_gravity_cartZ_generator(equations)
source_terms_gravity = transform_source_terms(source_terms_gravity_reference, equations)

###############################################################################
# Solver

surface_flux = FluxLMARS(340.0)

volume_flux = flux_kennedy_gruber
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

polydeg = 3
solver = DGSEM(polydeg, surface_flux, volume_integral)

###############################################################################
# Mesh

coordinates_min = (0.0, 0.0)
coordinates_max = (20_000.0, 10_000.0)

cells_per_dimension = (64, 32)
mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                      periodicity = (true, false))

###############################################################################
# Semidiscretization

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_gravity,
                                    boundary_conditions = boundary_conditions)

tspan = (0.0, 10.0)  # 1000 seconds final time

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
