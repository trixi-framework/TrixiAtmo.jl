using OrdinaryDiffEq
using Trixi
using TrixiAtmo
using TrixiAtmo: source_terms_no_phase_change



###  Bryan-Fritsch Moist Benchmark similar to Example 4.2 of:
# Sabine Doppler, Philip L. Lederer, Joachim Schöberl, Henry von Wahl,
# A discontinuous Galerkin approach for atmospheric flows with implicit condensation,
# Journal of Computational Physics,
# Volume 499,
# 2024,
# 112713,
# ISSN 0021-9991

function initial_condition_moist_bubble_bryan_fritsch(x, t, equations::CompressibleRainyEulerEquations2D) 
    # Position of the bubble defined by center and radius:
    x_center = 10_000.0
    z_center =  2_000.0
    
    x_radius =  2_000.0
    z_radius =  2_000.0

    # Hydrostatic Base State (compare Doppler et al. Appendix B.2)
    q_w      =     0.02               # total water fraction
    theta_e  =    320.0               # wet equivalent potential temperature
    
    p_ref = equations.ref_pressure

    L = min(sqrt(((x[1] - x_center) / x_radius)^2 + ((x[2] - z_center) / z_radius)^2), 1.0)

    #TODO

    return SVector(0, 0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, 0, 0, 0)
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

coordinates_min = (     0.0,      0.0)
coordinates_max = (20_000.0, 10_000.0)

cells_per_dimension = (64, 32)
mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                      periodicity = (true, false))

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_moist_bubble_bryan_fritsch, solver,
                                    source_terms = source_terms_no_phase_change,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

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