# This test case is used to compute convergence rates via a linearized solution.
# The setup follows the approach commonly adopted in benchmark studies; therefore,
# a fixed CFL number is employed.
#
# References:
# - Michael Baldauf and Slavko Brdar (2013):
#   "An analytic solution for linear gravity waves in a channel as a test
#   for numerical models using the non-hydrostatic, compressible Euler equations"
#   Q. J. R. Meteorol. Soc., DOI: 10.1002/qj.2105
#   https://doi.org/10.1002/qj.2105
#
# - Maciej Waruszewski, Jeremy E. Kozdon, Lucas C. Wilcox, Thomas H. Gibson,
#   and Francis X. Giraldo (2022):
#   "Entropy stable discontinuous Galerkin methods for balance laws
#   in non-conservative form: Applications to the Euler equations with gravity"
#   JCP, DOI: 10.1016/j.jcp.2022.111507
#   https://doi.org/10.1016/j.jcp.2022.111507
#
# - Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025):
#   "Structure-Preserving High-Order Methods for the Compressible Euler Equations
#   in Potential Temperature Formulation for Atmospheric Flows"
#   https://arxiv.org/abs/2509.10311

using OrdinaryDiffEqSSPRK
using Trixi
using TrixiAtmo: IdealGas,
                 initial_condition_gravity_waves_generator,
                 geopotential_cartZ
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
td_totE = EnergyTotal(td_single)

###############################################################################
# Equations

equations = CompressibleEulerAtmo(; n_dims = 2, n_vars_aux = 1,
                                  parameters = parameters,
                                  thermodynamic_state = td_single,
                                  thermodynamic_equation = td_totE)

###############################################################################
# Initial and boundary conditions

temperature_ref = 250.0
initial_condition_reference = initial_condition_gravity_waves_generator(parameters;
                                                                        temperature_ref = temperature_ref)
initial_condition = transform_initial_condition(initial_condition_reference, equations)

boundary_conditions = (; y_neg = boundary_condition_slip_wall_simple,
                       y_pos = boundary_condition_slip_wall_simple)

###############################################################################
# Solver

# We have an isothermal background state with T0 = 250 K.
# The reference speed of sound can be computed as:
# cs = sqrt(gamma * R * T0)
cs = sqrt(td_single.gamma_gas[1] * td_single.R_gas[1] * temperature_ref)
surface_flux = (FluxLMARS(cs), flux_zero)

volume_flux = (flux_ranocha, flux_nonconservative_waruszewski_etal)
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

polydeg = 3
solver = DGSEM(polydeg, surface_flux, volume_integral)

###############################################################################
# Mesh

coordinates_min = (0.0, 0.0)
coordinates_max = (300_000.0, 10_000.0)
trees_per_dimension = (60, 8)

mesh = P4estMesh(trees_per_dimension, polydeg = 1,
                 coordinates_min = coordinates_min, coordinates_max = coordinates_max,
                 periodicity = (true, false))

###############################################################################
# Semidiscretization

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_conditions,
                                    aux_field = geopotential_cartZ)

tspan = (0.0, 1800.0)

ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)
#                                     extra_analysis_integrals = (entropy,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = analysis_interval,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = "out_gravity_waves")

stepsize_callback = StepsizeCallback(cfl = 1.0)

visualization = VisualizationCallback(semi; interval = 10, show_mesh = false,
                                      solution_variables = cons2cons,
                                      variable_names = ["rho_v1", "rho_v2"],
                                      aspect_ratio = 10)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        visualization,
                        stepsize_callback)

sol = solve(ode,
            SSPRK43(thread = Trixi.Threaded());
            dt = 0.1, # solve needs some value here but it will be overwritten by the stepsize_callback
            ode_default_options()..., callback = callbacks, adaptive = false);
