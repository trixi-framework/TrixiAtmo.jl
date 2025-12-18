using OrdinaryDiffEq
using Trixi
using TrixiAtmo: cons2pot, source_terms_gravity


##############################################################################################
# initial condition

equations = CompressibleEulerEquations2D(1004.0 / 717.0)

@inline function initial_condition_straka(x,t,equations)
    #constants
    R_d = 287.0 #gas constatnt
    C_p = 1004.0 #specific heat at const pressure
    C_v = 717.0 #specific heat at const volume
    p_0 = 100_000.0 #reference pressure 
    g = 9.81 #gravity
    T_s = 300.0 #surface temperature

    #base state values 
    T_bar = T_s - x[2] * g / C_p 
    p_bar = p_0 * (T_bar * T_s ^ -1) ^ (C_p/R_d)
    theta_bar = T_s 
    v1 = 0.0 
    v2 = 0.0 

    # temperature pertubation
    x_c = 0.0 
    x_r = 4000.0
    z_c = 3000.0
    z_r = 2000.0

    L = (((x[1] - x_c) / x_r  ) ^ 2 + ((x[2] - z_c) / z_r ) ^ 2) ^ 0.5

    if L > 1.0 
        Del_T = 0.0 
    else 
        Del_T = -15.0 * (cos(pi*L) + 1.0) / 2.0
    end 

    T = T_bar + Del_T
    rho = p_bar / (R_d * T)
    E = C_v * T  + 0.5 * (v1^2 + v2^2)
    return SVector(rho, rho * v1, rho * v2, rho * E)
end

# source terms
@inline function source(u, x, t, equations::CompressibleEulerEquations2D)
    return (source_terms_gravity(u, equations::CompressibleEulerEquations2D))
end 

# Background reference state
function hydrostatic_background(x, t, equations::CompressibleEulerEquations2D)
    #constants
    R_d = 287.0 #gas constatnt
    C_p = 1004.0 #specific heat at const pressure
    C_v = 717.0 #specific heat at const volume
    p_0 = 100_000.0 #reference pressure 
    g = 9.81 
    T_s = 300.0 #surface temperature

    #base state values 
    T_bar = T_s - x[2] * g * C_p ^ -1
    p_bar = p_0 * (T_bar * T_s ^ -1) ^ (C_p/R_d)
    rho = p_bar * (R_d * T_bar) ^ -1
    v1 = 0.0 
    v2 = 0.0 

    E = C_v * T_bar + 0.5 * (v1^2 + v2^2)
    return SVector(rho, rho * v1, rho * v2, rho * E)
end


##############################################################################################
# semidiscretization 

# mesh 
coordinates_min = (-25_600.0, 0.0)
coordinates_max = ( 25_600.0, 6_400.0)
cells_per_dimension = (512, 64)
mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                      periodicity = (false))

initial_condition = initial_condition_straka

# flux and solver 
polydeg = 3
surface_flux = flux_lax_friedrichs 
solver = DGSEM(polydeg = polydeg,
               surface_flux = surface_flux,
               volume_integral = VolumeIntegralWeakForm())

source_term = source

boundary_conditions =  (x_neg = boundary_condition_slip_wall,
                        x_pos = boundary_condition_slip_wall,
                        y_neg = boundary_condition_slip_wall,
                        y_pos = boundary_condition_slip_wall)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_term,
                                    boundary_conditions = boundary_conditions)


###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 900.0)

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
solution_variables = cons2pot

analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:entropy_conservation_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

background_state = compute_reference_state(hydrostatic_background, semi, cons2pot)

save_solution = SaveSolutionCallback(interval = analysis_interval,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = "out",
                                     solution_variables = solution_variables,
                                     reference_solution = background_state)

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