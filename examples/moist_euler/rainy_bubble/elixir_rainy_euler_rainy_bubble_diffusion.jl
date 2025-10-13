using OrdinaryDiffEqSSPRK
using Trixi
using TrixiAtmo
using TrixiAtmo: source_terms_rainy, initial_condition_bubble_rainy,
                 saturation_residual, saturation_residual_jacobian,
                 cons2eq_pot_temp, saturation_vapour_pressure,
                 flux_chandrashekar, flux_LMARS,
                 source_terms_no_phase_change,
                 boundary_condition_laplace,
                 boundary_condition_simple_slip_wall,
                 flux_ec_rain
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

diffusivity = 0.5f0
equations_parabolic = LaplaceDiffusion2D(diffusivity, equations)

boundary_conditions = (x_neg = boundary_condition_periodic,
                       x_pos = boundary_condition_periodic,
                       y_neg = boundary_condition_simple_slip_wall,
                       y_pos = boundary_condition_simple_slip_wall)

boundary_conditions_parabolic = (x_neg = boundary_condition_periodic,
                                 x_pos = boundary_condition_periodic,
                                 y_neg = boundary_condition_laplace,
                                 y_pos = boundary_condition_laplace)

polydeg = 3

surface_flux = flux_lax_friedrichs
volume_integral = VolumeIntegralFluxDifferencing(flux_ec_rain)

solver = DGSEM(polydeg, surface_flux, volume_integral)

initial_refinement_level = 6
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = initial_refinement_level,
                periodicity = (true, false), n_cells_max = 1_000_000)

semi = SemidiscretizationHyperbolicParabolic(mesh, (equations, equations_parabolic),
                                             initial_condition_rainy, solver;
                                             source_terms_rainy,
                                             boundary_conditions = (boundary_conditions,
                                                                    boundary_conditions_parabolic))

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
