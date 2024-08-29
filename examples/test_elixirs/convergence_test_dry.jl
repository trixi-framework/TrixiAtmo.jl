using OrdinaryDiffEq
using Trixi
using TrixiAtmo



# from compressible_euler_2d.jl
function initial_condition_convergence_test(x, t, equations::CompressibleRainyEulerEquations2D)
    RealT = eltype(x)
    c = 2
    A = convert(RealT, 0.1)
    L = 2
    f = 1.0f0 / L
    ω = 2 * convert(RealT, pi) * f
    ini = c + A * sin(ω * (x[1] + x[2] - t))

    rho = ini
    rho_v1 = ini
    rho_v2 = ini
    rho_e = ini^2

    return SVector(rho, 0.0, 0.0, rho_v1, rho_v2, rho_e, 0.0, 0.0, 0.0)
end


# from compressible_euler_2d.jl
function source_terms_convergence_test(u, x, t, equations::CompressibleRainyEulerEquations2D)
    # Same settings as in `initial_condition`
    RealT = eltype(u)
    c = 2
    A = convert(RealT, 0.1)
    L = 2
    f = 1.0f0 / L
    ω = 2 * convert(RealT, pi) * f
    γ = equations.c_dry_air_const_pressure / equations.c_dry_air_const_volume

    x1, x2 = x
    si, co = sincos(ω * (x1 + x2 - t))
    rho = c + A * si
    rho_x = ω * A * co
    # Note that d/dt rho = -d/dx rho = -d/dy rho.

    tmp = (2 * rho - 1) * (γ - 1)

    du1 = rho_x
    du2 = rho_x * (1 + tmp)
    du3 = du2
    du4 = 2 * rho_x * (rho + tmp)

    return SVector(du1, 0.0, 0.0, du2, du3, du4, 0.0, 0.0, 0.0)
end


# from elixir_euler_source_terms.jl
###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleRainyEulerEquations2D()

initial_condition = initial_condition_convergence_test

solver = DGSEM(polydeg = 3, surface_flux = flux_lax_friedrichs)

coordinates_min = (0.0, 0.0)
coordinates_max = (2.0, 2.0)

cells_per_dimension = (16, 16)

mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_convergence_test)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 2.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 100,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks);
summary_callback() # print the timer summary

# For copy-paste convenience:
#convergence_test("TrixiAtmo.jl/examples/test_elixirs/convergence_test_dry.jl", 4)