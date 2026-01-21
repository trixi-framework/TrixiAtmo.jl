# References:
# - André Robert (1993)
#   Bubble Convection Experiments with a Semi-implicit Formulation of the Euler Equations
#   Journal of the Atmospheric Sciences, Volume 50: Issue 13
#   https://doi.org/10.1175/1520-0469(1993)050%3C1865:BCEWAS%3E2.0.CO;2

using OrdinaryDiffEqSSPRK
using Trixi, TrixiAtmo
using TrixiAtmo: initial_condition_dry_air_warm_bubble_generator, 
                 source_terms_gravity_cartZ_generator
using Plots

RealType = Float64
n_dim = 2

parameters = Parameters{RealType}(;
    earth_gravitational_acceleration = 9.81
)
td_single = IdealGas(; parameters)
td_potT = PotentialTemperature(td_single)

microphysics = MicrophysicsRelaxation{RealType}()

equations = CompressibleEulerAtmo{n_dim}(; parameters = parameters,
                                           thermodynamic_state = td_single,
                                           thermodynamic_equation = td_potT,
                                           microphysics = microphysics)
#surface_flux = flux_lax_friedrichs
surface_flux = FluxLMARS(340)
#volume_flux = flux_tec
polydeg = 3
solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux)
               #volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

boundary_conditions = (x_neg = boundary_condition_periodic,
                       x_pos = boundary_condition_periodic,
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

coordinates_min = (0.0, 0.0)
coordinates_max = (1000.0, 1500.0)
cells_per_dimension = (16, 16)
mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                      periodicity = (true, false))

initial_condition_reference = initial_condition_dry_air_warm_bubble_generator(
    parameters;
    perturbation_center_x = 500.0,
    perturbation_center_z = 260.0,
    perturbation_radius = 250.0
)
initial_condition = transform_initial_condition(initial_condition_reference, equations)

source_terms_gravity_reference = source_terms_gravity_cartZ_generator(equations)
source_terms_gravity = transform_source_terms(source_terms_gravity_reference, equations)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_gravity,
                                    boundary_conditions = boundary_conditions)

tspan = (0.0, 10)
ode = semidiscretize(semi, tspan)

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
                        visualization,save_solution,
                        alive_callback)

sol = solve(ode,
            SSPRK43(thread = Trixi.True());
            maxiters = 1.0e7, ode_default_options()..., callback = callbacks)
