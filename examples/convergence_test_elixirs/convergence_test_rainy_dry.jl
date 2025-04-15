using OrdinaryDiffEq
using Trixi
using TrixiAtmo



# adapted from compressible_euler_2d.jl
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


# adapted from compressible_euler_2d.jl
function source_terms_convergence_test(u, x, t, equations::CompressibleRainyEulerEquations2D)
    c_pd = equations.c_dry_air_const_pressure
    c_vd = equations.c_dry_air_const_volume

    # Same settings as in `initial_condition`
    RealT = eltype(u)
    c = 2
    A = convert(RealT, 0.1)
    L = 2
    f = 1.0f0 / L
    ω = 2 * convert(RealT, pi) * f
    γ = c_pd / c_vd

    x1, x2 = x
    si, co = sincos(ω * (x1 + x2 - t))
    rho = c + A * si
    rho_x = ω * A * co
    # Notice that d/dt rho = -d/dx rho = -d/dy rho.

    tmp = (2 * rho - 1) * (γ - 1) #+ R_d * ref_temp

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

cells_per_dimension = (8, 8)

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

save_solution = SaveSolutionCallback(interval = analysis_interval,
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


#= For copy-paste convenience:
#convergence_test("TrixiAtmo.jl/examples/convergence_test_elixirs/convergence_test_rainy_dry.jl", 5)

####################################################################################################
l2
rho_dry             rho_moist           rho_rain            rho_v1              rho_v2              energy_density      rho_vapour          rho_cloud           temperature
error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     
EOC
2.80e-05  -         0.00e+00  -         0.00e+00  -         2.90e-05  -         2.90e-05  -         9.66e-05  -         0.00e+00  -         0.00e+00  -         0.00e+00  
-
9.30e-07  4.91      0.00e+00  NaN       0.00e+00  NaN       1.42e-06  4.36      1.42e-06  4.36      4.82e-06  4.32      0.00e+00  NaN       0.00e+00  NaN       0.00e+00  
NaN
7.00e-08  3.73      0.00e+00  NaN       0.00e+00  NaN       9.51e-08  3.90      9.51e-08  3.90      3.29e-07  3.87      0.00e+00  NaN       0.00e+00  NaN       0.00e+00  
NaN
4.64e-09  3.92      0.00e+00  NaN       0.00e+00  NaN       6.08e-09  3.97      6.08e-09  3.97      2.12e-08  3.95      0.00e+00  NaN       0.00e+00  NaN       0.00e+00  
NaN
2.93e-10  3.98      0.00e+00  NaN       0.00e+00  NaN       3.82e-10  3.99      3.82e-10  3.99      1.33e-09  3.99      0.00e+00  NaN       0.00e+00  NaN       0.00e+00  
NaN

mean      4.14      mean      NaN       mean      NaN       mean      4.05      mean      4.05      mean      4.04      mean      NaN       mean      NaN       mean      
NaN
----------------------------------------------------------------------------------------------------
linf
rho_dry             rho_moist           rho_rain            rho_v1              rho_v2              energy_density      rho_vapour          rho_cloud           temperature
error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     EOC       error     
EOC
1.97e-04  -         0.00e+00  -         0.00e+00  -         2.02e-04  -         2.02e-04  -         8.91e-04  -         0.00e+00  -         0.00e+00  -         0.00e+00  
-
9.62e-06  4.36      0.00e+00  NaN       0.00e+00  NaN       1.17e-05  4.10      1.17e-05  4.10      4.89e-05  4.19      0.00e+00  NaN       0.00e+00  NaN       0.00e+00  
NaN
6.24e-07  3.95      0.00e+00  NaN       0.00e+00  NaN       7.49e-07  3.97      7.49e-07  3.97      3.23e-06  3.92      0.00e+00  NaN       0.00e+00  NaN       0.00e+00  
NaN
4.06e-08  3.94      0.00e+00  NaN       0.00e+00  NaN       4.97e-08  3.91      4.97e-08  3.91      2.10e-07  3.94      0.00e+00  NaN       0.00e+00  NaN       0.00e+00  
NaN
2.55e-09  3.99      0.00e+00  NaN       0.00e+00  NaN       3.19e-09  3.96      3.19e-09  3.96      1.34e-08  3.97      0.00e+00  NaN       0.00e+00  NaN       0.00e+00  
NaN

mean      4.06      mean      NaN       mean      NaN       mean      3.99      mean      3.99      mean      4.00      mean      NaN       mean      NaN       mean      
NaN
----------------------------------------------------------------------------------------------------
Dict{Symbol, Any} with 3 entries:
  :variables => ("rho_dry", "rho_moist", "rho_rain", "rho_v1", "rho_v2", "energy_density", "rho_vapour", "rho_cloud", "temperature")
  :l2        => [4.13577, NaN, NaN, 4.05385, 4.05385, 4.03601, NaN, NaN, NaN]
  :linf      => [4.05911, NaN, NaN, 3.9867, 3.98669, 4.00476, NaN, NaN, NaN]=#