using OrdinaryDiffEq
using Trixi



function initial_condition_bubble_no_rain(x, t, equations::CompressibleRainyEulerEquations2D)
    #TODO
    return
end


###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleRainyEulerEquations2D()

#TODO check boundary conditions
boundary_conditions = (x_neg = boundary_condition_periodic,
                       x_pos = boundary_condition_periodic,
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

polydeg = 3
basis = LobattoLegendreBasis(polydeg)

surface_flux = flux_lax_friedrichs

solver = DGSEM(basis, surface_flux)

#TODO adjust values
coordinates_min = (0.0, 0.0)
coordinates_max = (1.0, 1.0)

cells_per_dimension = (64, 32)
mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                      periodicity = (true, false))

#TODO import source terms from equations
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_bubble_no_rain, solver,
                                    source_terms = source_terms_no_rain,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

#TODO adjust values
tspan = (0.0, 1000.0)  # 1000 seconds final time

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000

analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:entropy_conservation_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = analysis_interval,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = "out",
                                     solution_variables = cons2prim)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            maxiters = 1.0e7,
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks);

summary_callback()