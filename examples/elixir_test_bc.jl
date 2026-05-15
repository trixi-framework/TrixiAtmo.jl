using OrdinaryDiffEq, Trixi
using TrixiAtmo: cons2temppert
using TrixiAtmo


##############################################################################################

equations = PerturbationEulerEquations2DAuxVars(1004.0 / 717.0)

# initial condition with perturbation
@inline function initial_condition_gravity_wave(x, t, equations)
    # constants 
    g = 9.81 
    c_p = 1004.0 
    c_v = 717.0

    theta_c = 0.01 
    x_c = 100_000.0 # x-coordinate of center point of the perturbated area 
    h_c = 10_000.0 # y-coordinate
    a_c = 5_000.0 # width 
    theta_0 = 300.0 # constant of integration 
    bvfrequency = 0.01 # Brunt-Väisälä frequency
    p_0 = 100_000.0  # reference pressure
    R = c_p - c_v


    # exner pressure from hydrostatic balance
    exner = 1 + g^2 / (c_p * theta_0 * bvfrequency^2) *
            (exp(-bvfrequency^2 / g * x[2]) - 1)

    # potential temperature: perturbated, background and total 
    theta_prime = theta_c * (sin(pi * x[2] / h_c)) /
                                      (1 + ((x[1] - x_c) / a_c)^2)
    theta_mean = theta_0 * exp(bvfrequency^2 / g * x[2])
    theta = theta_prime + theta_mean


    # density: background, total and perturbation
    rho_mean = p_0 / (R * theta_mean) * exner ^ (c_v / R) 
    rho = p_0 / (R * theta) * exner ^ (c_v /R) 
    rho_prime = rho - rho_mean

    # total and mean velocities  
    v1, v2 = 20.0, 0.0
    v1_mean, v2_mean = 20.0, 0.0

    rhov1_prime = rho * v1 - rho_mean * v1_mean 
    rhov2_prime = rho * v2 - rho_mean * v2_mean


    # energy: total and background  
    e = c_v * theta * exner + 0.5 * (v1^2 + v2^2)
    e_mean = c_v * theta_mean * exner + 0.5 * (v1^2 + v2^2) #v1_mean and v2_mean, but is identical here

    rhoe_prime = rho * e - rho_mean * e_mean

    return SVector(rho_prime, rhov1_prime, rhov2_prime, rhoe_prime)
end


# auxiliary field with background state
@inline function background_state(x)
    g = 9.81 
    c_p = 1004.0 
    c_v = 717.0

    # Exner pressure from hydrostatic balance for x[2]
    theta_0 = 300.0 # constant of integration 
    bvfrequency = 0.01 # Brunt-Väisälä frequency
    # pressure
    p_0 = 100_000.0  # reference pressure
    R = c_p - c_v    # gas constant (dry air)

    exner = 1 +
            g^2 / (c_p * theta_0 * bvfrequency^2) *
            (exp(-bvfrequency^2 / g * x[2]) - 1)

    # potential temperature
    theta_mean = theta_0 * exp(bvfrequency^2 / g * x[2])

    # density
    rho_mean = p_0 / (R * theta_mean) * exner ^ (c_v / R)
 
    # velocity 
    v1_mean, v2_mean = 20.0, 0.0

    # energy  
    e_mean = c_v * theta_mean * exner + 0.5 * (v1_mean^2 + v2_mean^2) 

    return SVector(rho_mean, rho_mean*v1_mean, rho_mean*v2_mean, e_mean)
end 


# Source terms   
@inline function source_terms(u, aux, x, t, equations::PerturbationEulerEquations2DAuxVars)
    g = 9.81
    rho = u[1] #+ aux[1]
    rho_v2 = u[3] + aux[3]
    return SVector(zero(eltype(u)), zero(eltype(u)), -g * rho, -g * rho_v2)
end



##############################################################################################
# semidiscretization 

polydeg = 3
coordinates_min = (0.0, -10_000.0)
coordinates_max = (300_000.0, 20_000.0)

cells_per_dimension = (300, 10)
mesh = P4estMesh(cells_per_dimension; polydeg = 2, coordinates_min, coordinates_max,
                 periodicity = (true, false))

boundary_conditions = Dict(:y_neg => boundary_condition_slip_wall_aux,
                       :y_pos => boundary_condition_slip_wall_aux)

              
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

tspan = (0.0, 1000.0)  # 3000 seconds final time

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
solution_variables = cons2all

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

