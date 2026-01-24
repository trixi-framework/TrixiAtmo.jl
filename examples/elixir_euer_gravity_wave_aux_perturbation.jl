using OrdinaryDiffEq, Trixi
using TrixiAtmo

##############################################################################################

equations = PerturbationEulerEquations2DAuxVars(1.4)


# initial condition with perturbation
@inline function initial_condition_gravity_wave(x, t, equations)
    # constants 
    g = 9.81 
    c_p = 1004.0 
    c_v = 717.0

    theta_c = 0.01 
    h_c = 10_000.0
    a_c = 5_000.0
    x_c = 100_000.0 
    theta_0 = 300.0 # constant of integration 
    bvfrequency = 0.01 # Brunt-Väisälä frequency
    p_0 = 100_000.0  # reference pressure
    R = c_p - c_v


    # exner pressure 
    exner = 1 + g^2 / (c_p * theta_0 * bvfrequency^2) *
            (exp(-bvfrequency^2 / g * x[2]) - 1)

    # potential temperature: perturbated, background and total 
    theta_prime = theta_c * sin(pi * x[2] / h_c) /
                                      (1 + ((x[1] - x_c) / a_c)^2)
    theta_mean = theta_0 * exp(bvfrequency^2 / g * x[2])
    theta = theta_prime + theta_mean


    # density: background, total and perturbation
    rho_mean = p_0 / (R * theta_mean) * exner ^ (c_v / R) 
    rho = p_0 / (R * (theta_mean + theta_prime)) * exner ^ (c_v /R) 
    rho_prime = rho - rho_mean

    # velocity 
    v1, v2 = 20.0, 0.0

    # energy: total and background  
    E = c_v * theta * exner + 0.5 * (v1^2 + v2^2)
    E_mean = c_v * theta_mean * exner + 0.5 * (v1^2 + v2^2) 


    # conservative variables 
    rhoE_prime = rho * E - rho_mean * E_mean
    rho_v1 = rho * v1 
    rho_v2 = rho * v2 

    return SVector(rho_prime, rho_v1, rho_v2, rhoE_prime)
end


# auxiliary field with background state
@inline function background_state(x)
    g = 9.81 
    c_p = 1004.0 
    c_v = 717.0

    # Exner pressure from hydrostatic balance for x[2]
    theta_0 = 300.0 # constant of integration 
    bvfrequency = 0.01 # Brunt-Väisälä frequency

    exner = 1 +
            g^2 / (c_p * theta_0 * bvfrequency^2) *
            (exp(-bvfrequency^2 / g * x[2]) - 1)

    # potential temperature
    theta_mean = theta_0 * exp(bvfrequency^2 / g * x[2])


    # pressure
    p_0 = 100_000.0  # reference pressure
    R = c_p - c_v    # gas constant (dry air)

    # density
    rho_mean = p_0 / (R * theta_mean) * exner ^ (c_v / R)
 
    v1 = 20.0
    v2 = 0.0
    E_mean = c_v * theta_mean * exner + 0.5 * (v1^2 + v2^2) 
    return SVector(rho_mean, rho_mean * E_mean)
end 


# Source terms  
@inline function source_terms(u, x, t, equations::PerturbationEulerEquations2DAuxVars)
    g = 9.81
    rho_prime = u[1]
    return SVector(zero(eltype(u)), zero(eltype(u)), -g * rho_prime, zero(eltype(u)))
end



##############################################################################################
# semidiscretization 

polydeg = 4 
coordinates_min = (0.0, 0.0)
coordinates_max = (300_000.0, 10_000.0)

cells_per_dimension = (300, 10)
mesh = P4estMesh(cells_per_dimension; polydeg = 2, coordinates_min, coordinates_max)

boundary_conditions = (x_neg = boundary_condition_periodic,
                       x_pos = boundary_condition_periodic,
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)
              
initial_condition = initial_condition_gravity_wave

basis = LobattoLegendreBasis(polydeg)


surface_flux = flux_lax_friedrichs
volume_integral = VolumeIntegralWeakForm()

solver = DGSEM(basis, surface_flux, volume_integral)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms,
                                    boundary_conditions = boundary_conditions,
                                    aux_field = background_state)


###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 3000.0)  # 3000 seconds final time

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
solution_variables = cons2prim

analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = analysis_interval,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = "out",
                                     solution_variables = solution_variables)

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

