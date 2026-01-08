using OrdinaryDiffEqLowStorageRK
using Trixi, TrixiAtmo
using TrixiAtmo: cons2aeqpot, saturation_pressure, initial_condition_moist_bubble,
                 source_terms_moist_bubble
using NLsolve: nlsolve

###############################################################################
# semidiscretization of the compressible moist Euler equations

c_pd = 1004 # specific heat at constant pressure for dry air
c_vd = 717  # specific heat at constant volume for dry air
c_pv = 1885 # specific heat at constant pressure for moist air
c_vv = 1424 # specific heat at constant volume for moist air
equations = CompressibleMoistEulerEquations2D(c_pd = c_pd, c_vd = c_vd, c_pv = c_pv,
                                              c_vv = c_vv,
                                              gravity = EARTH_GRAVITATIONAL_ACCELERATION)

# Create background atmosphere data set
atmosphere_data = AtmosphereLayers(equations)

# Create the initial condition with the initial data set
function initial_condition_moist(x, t, equations)
    rho, rho_v1, rho_v2, rho_E, rho_qv, rho_ql = initial_condition_moist_bubble(x, t,
                                                                                equations_moist,
                                                                                atmosphere_data)
    return SVector(rho, rho_v1, rho_v2, rho_E, rho_qv, rho_ql)
end

initial_condition = initial_condition_moist

boundary_condition = (x_neg = boundary_condition_slip_wall,
                      x_pos = boundary_condition_slip_wall,
                      y_neg = boundary_condition_slip_wall,
                      y_pos = boundary_condition_slip_wall)

source_term = source_terms_moist_bubble

###############################################################################
# Get the DG approximation space

polydeg = 4

surface_flux = FluxLMARS(360.0)
volume_flux = flux_chandrashekar

volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(polydeg, surface_flux, volume_integral)

coordinates_min = (0.0, 0.0)
coordinates_max = (20000.0, 10000.0)

cells_per_dimension = (64, 32)

# Create curved mesh with 64 x 32 elements
mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                      periodicity = (false, false))

###############################################################################
# create the semi discretization object

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition,
                                    source_terms = source_term)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 1000.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
solution_variables = cons2aeqpot

analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:entropy_conservation_error,),
                                     extra_analysis_integrals = (entropy, energy_total,
                                                                 saturation_pressure))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 1000,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = solution_variables)

stepsize_callback = StepsizeCallback(cfl = 0.2)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks)
