using OrdinaryDiffEqSSPRK
using Trixi
using TrixiAtmo
using TrixiAtmo: source_terms_rainy, initial_condition_bubble_rainy,
                 saturation_residual, saturation_residual_jacobian,
                 cons2eq_pot_temp, saturation_vapour_pressure,
                 flux_ec_rain, boundary_condition_simple_slip_wall,
                 generate_hydrostatic_residual, generate_perturbation_residual
using NLsolve: nlsolve

# domain
coordinates_min = (0.0, 0.0)
coordinates_max = (2400.0, 2400.0)

# create layers for initial condition
equations = CompressibleRainyEulerEquations2D()
atmosphere_data = AtmosphereLayersRainyBubble(equations; total_height = coordinates_max[2] + 1)

# Create the initial condition with the initial data set
function initial_condition_rainy(x, t, equations::CompressibleRainyEulerEquations2D)
    return initial_condition_bubble_rainy(x, t, equations; atmosphere_data)
end

###############################################################################
# semidiscretization of the compressible rainy Euler equations

boundary_conditions = (x_neg = boundary_condition_periodic,
                       x_pos = boundary_condition_periodic,
                       y_neg = boundary_condition_simple_slip_wall,
                       y_pos = boundary_condition_simple_slip_wall)

polydeg = 3

surface_flux = flux_lax_friedrichs
volume_integral = VolumeIntegralFluxDifferencing(flux_ec_rain)

solver = DGSEM(polydeg, surface_flux, volume_integral)

cells_per_dimension = (64, 64)
mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                      periodicity = (true, false))

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_rainy, solver,
                                    source_terms = source_terms_rainy,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 600.0)

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000

# entropy
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:entropy_conservation_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 1000,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = "out",
                                     solution_variables = cons2eq_pot_temp)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution)

stage_limiter! = NonlinearSolveDG(saturation_residual, saturation_residual_jacobian,
                                  SVector(7, 8, 9))

###############################################################################
# run the simulation
sol = solve(ode, SSPRK43(stage_limiter!); ode_default_options()...,
            maxiters = 1.0e7, save_everystep = false, callback = callbacks);
