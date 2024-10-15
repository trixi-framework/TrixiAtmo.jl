using OrdinaryDiffEq
using Trixi
using TrixiAtmo
using TrixiAtmo: source_terms_rainy, saturation_residual,
                 saturation_residual_jacobian, NonlinearSolveDG,
                 cons2eq_pot_temp, saturation_vapour_pressure
using NLsolve: nlsolve


## Hydrostatic base state

function initial_constraint_residual(guess, equations::CompressibleRainyEulerEquations2D)
    # rho_dry, rho_vapour, temperature = guess
    c_pd          = equations.c_dry_air_const_pressure
    R_d           = equations.R_dry_air
    R_v           = equations.R_vapour
    eps           = equations.eps
    ref_p         = equations.ref_pressure
    surf_p        = 8.5e4                   # surface pressure
    humidity_rel0 = 0.2                     # hydrostatic relative humidity

    rho_vs = saturation_vapour_pressure(guess[3]) / (R_v * guess[3])

    res1  = surf_p   - guess[3] * (guess[1] * R_d + guess[2] * R_v)
    res2  = theta_d0 - guess[3] * (ref_p / surf_p)^(R_d / c_pd)
    res3  = rho_vs * (guess[1] + guess[2] / eps) * humidity_rel0
    res3 -= guess[2] * (guess[1] + rho_vs / eps)

    return SVector(res1, res2, 1000 * res3)
end


## Initial condition

function initial_condition_bubble_rainy(x, t, equations::CompressibleRainyEulerEquations2D)
    #TODO
    return SVector(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
end



###############################################################################
# semidiscretization of the compressible rainy Euler equations

equations = CompressibleRainyEulerEquations2D()

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

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_bubble_rainy, solver,
                                    source_terms = source_terms_rainy,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

#TODO adjust values
tspan = (0.0, 1000.0)

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000

# entropy?
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
            dt = 1.0,
            save_everystep = false, callback = callbacks);

summary_callback()