
using OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiAtmo

###############################################################################
# semidiscretization of the compressible Euler equations
equations = CompressibleEulerEquationsWithGravityNoPressure2D(1.4)

function initial_condition_constant(x, t,
                                    equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    rho = exp(-x[2])
    v1 = 0.0
    v2 = 0.0
    p = exp(-x[2])
    prim = SVector(rho, v1, v2, p, x[2], 1.0)
    return prim2cons(prim, equations)
end

initial_condition = initial_condition_constant

volume_flux = (flux_kennedy_gruber, flux_nonconservative_chandrashekar_isothermal)

surface_flux = (FluxHydrostaticReconstruction(FluxPlusDissipation(flux_kennedy_gruber,
                                                                  DissipationLocalLaxFriedrichs()),
                                              hydrostatic_reconstruction),
                FluxHydrostaticReconstruction(flux_nonconservative_chandrashekar_isothermal,
                                              hydrostatic_reconstruction))

polydeg = 0
basis = LobattoLegendreBasis(polydeg)
limiter_idp = SubcellLimiterIDP(equations, basis;
                                positivity_variables_cons = ["rho"],)
#  volume_integral = VolumeIntegralSubcellLimiting(limiter_idp;
#                                                  volume_flux_dg = volume_flux,
#                                                  volume_flux_fv = surface_flux)
#
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
#volume_integral = VolumeIntegralPureLGLFiniteVolume(volume_flux_fv = surface_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

coordinates_min = (0.0, 0.0)
coordinates_max = (1.0, 1.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 4,
                n_cells_max = 10_000,
                periodicity = false)

boundary_conditions = (; x_neg = boundary_condition_slip_wall,
                       x_pos = boundary_condition_slip_wall,
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 100.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     analysis_polydeg = polydeg,
                                     save_analysis = true)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

mean_temperature_callback = MeanTemperatureCallback()

save_solution = SaveSolutionCallback(dt = 10,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim) #extra_node_variables = (:limiting_coefficient,)

stepsize_callback = StepsizeCallback(cfl = 0.5)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        mean_temperature_callback,
                        save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

stage_callbacks = () #SubcellLimiterIDPCorrection(),

sol = Trixi.solve(ode, Trixi.SimpleSSPRK33(stage_callbacks = stage_callbacks);
                  dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
                  ode_default_options()...,
                  callback = callbacks);
