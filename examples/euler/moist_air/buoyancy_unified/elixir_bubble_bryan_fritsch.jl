using OrdinaryDiffEqLowStorageRK
using OrdinaryDiffEqSSPRK
using TrixiAtmo
using Plots


###############################################################################
# Parameters

RealType = Float64
n_dim = 2
n_condens = 0

parameters = Parameters{RealType}()


###############################################################################
# Thermodynamics

td_state = TrixiAtmo.IdealGas(; parameters)
#td_state = IdealGasesAndLiquids(; parameters,
#    n_gas = 1 + n_condens,
#    n_condens = n_condens)

td_equation = PotentialTemperature(td_state)
#td_equation = TotalEnergy(td_state)


###############################################################################
# Microphysics

microphysics = MicrophysicsRelaxation{RealType}(;
    saturation_factor = 0.0)
    #saturation_factor = 1.0)


###############################################################################
# Equations

equations = CompressibleEulerAtmo{n_dim}(; parameters = parameters,
                                        thermodynamic_state = td_state,
                                        thermodynamic_equation = td_equation,
                                        microphysics = microphysics)


###############################################################################
# Initial and boundary conditions

initial_condition_reference = initial_condition_bryan_fritsch_bubble_generator(
    parameters;
    moist = (n_condens != 0),
    perturbation_center_x = 10000.0,
    perturbation_center_z = 2000.0,
    perturbation_radius = 2000.0
)
initial_condition = transform_initial_condition(initial_condition_reference, equations)


###############################################################################
# Boundary conditions

boundary_conditions = (; x_neg = boundary_condition_slip_wall,
                         x_pos = boundary_condition_slip_wall,
                         y_neg = boundary_condition_slip_wall,
                         y_pos = boundary_condition_slip_wall)


###############################################################################
# Source terms

source_terms_gravity_reference = source_terms_gravity_cartZ_generator(equations)

if n_condens == 0
    source_terms = transform_source_terms(source_terms_gravity_reference, equations)
else
    source_terms_phase_change_reference = source_terms_phase_change_generator(microphysics)
    source_terms = transform_source_terms_sum((source_terms_gravity_reference,
                                               source_terms_phase_change_reference),
                                              equations)
end


###############################################################################
# Mesh

coordinates_min = (0.0, 0.0)
coordinates_max = (20000.0, 10000.0)

cells_per_dimension = (32, 16)

mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                      periodicity = false)


###############################################################################
# Solver

#surface_flux = flux_lax_friedrichs
surface_flux = FluxLMARS(360.0)
#volume_flux = flux_ranocha
#volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

polydeg = 3
solver = DGSEM(polydeg, surface_flux) #, volume_integral)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms,
                                    boundary_conditions = boundary_conditions)


###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 1000.0)  # 1000 seconds final time

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000

analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = analysis_interval,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = "out_warm_bubble")

visualization = VisualizationCallback(semi; interval = 100,
                                      #variable_names = ["v1", "r_vapor"],
                                      show_mesh = false)

#stepsize_callback = StepsizeCallback(cfl = 0.5)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        #save_solution,
                        #stepsize_callback,
                        visualization)
                        

sol = solve(ode,
            #CarpenterKennedy2N54(williamson_condition = false); dt = 0.1,
            SSPRK43(thread = Trixi.True());
            ode_default_options()..., callback = callbacks);
