using OrdinaryDiffEq
using Trixi
using TrixiAtmo: cons2pot


#Inertia-gravity waves test case 

struct InteriaWaveSetup
    # Physical constants
    g::Float64       # gravity of earth
    c_p::Float64     # heat capacity for constant pressure (dry air)
    c_v::Float64     # heat capacity for constant volume (dry air)
    gamma::Float64   # heat capacity ratio (dry air)

    function InteriaWaveSetup(; g = 9.81, c_p = 1004.0, c_v = 717.0, gamma = c_p / c_v)
        new(g, c_p, c_v, gamma)
    end
end

# Initial condition
function (setup::InteriaWaveSetup)(x, t,
                                   equations::CompressibleEulerEquations2D)
    @unpack g, c_p, c_v = setup

    # Exner pressure from hydrostatic balance for x[2]
    theta_0 = 300.0 #constant of integration 
    bvfrequency = 0.01 #Brunt-Väisälä frequency

    exner = 1 +
            g^2 / (c_p * theta_0 * bvfrequency^2) *
            (exp(-bvfrequency^2 / g * x[2]) - 1)

    # mean potential temperature
    theta_mean = theta_0 * exp(bvfrequency^2 / g * x[2])

    # perturbed potential temperature
    theta_c = 0.01
    h_c = 10_000.0
    a_c = 5_000.0
    x_c = 100_000.0

    theta_prime = theta_c * sin(pi * x[2] / h_c) /
                                      (1 + ((x[1] - x_c) / a_c)^2)

    # potential temperature
    theta = theta_mean + theta_prime 

    # temperature
    T = theta * exner

    # pressurey    
    p_0 = 100_000.0  # reference pressure
    R = c_p - c_v    # gas constant (dry air)
    #< p = p_0 * exner^(c_p / R)

    # density
    rho = p_0 / (R * theta) * exner ^ (c_v / R)

    v1 = 20.0
    v2 = 0.0
    E = c_v * T + 0.5 * (v1^2 + v2^2) 
    return SVector(rho, rho * v1, rho * v2, rho * E)
end


# Source terms  
@inline function (setup::InteriaWaveSetup)(u, x, t, equations::CompressibleEulerEquations2D)
    @unpack g = setup
    rho, _, rho_v2, _ = u
    return SVector(zero(eltype(u)), zero(eltype(u)), -g * rho, -g * rho_v2)

end

# Background reference state
function hydrostatic_background(x, t, equations::CompressibleEulerEquations2D)
    g = 9.81
    c_p = 1004.0f0
    c_v = 717.0f0    
    theta_0 = 300.0 #constant of integration 
    bvfrequency = 0.01 #Brunt-Väisälä frequency

    exner = 1 +
            g^2 / (c_p * theta_0 * bvfrequency^2) *
            (exp(-bvfrequency^2 / g * x[2]) - 1)

    # mean potential temperature
    theta = theta_0 * exp(bvfrequency^2 / g * x[2])

    # pressure
    p_0 = 100_000.0  # reference pressure
    R = c_p - c_v    # gas constant (dry air)
    #p = p_0 * exner^(c_p / R)

    # temperature
    T = theta * exner

   # density
    rho = p_0 / (R * theta) * exner ^ (c_v / R)

    v1 = 20.0
    v2 = 0.0
    E = c_v * T + 0.5 * (v1^2 + v2^2) 
    return SVector(rho, rho * v1, rho * v2, rho * E)
end

###############################################################################
# semidiscretization of the compressible Euler equations
interia_wave_setup = InteriaWaveSetup()

equations = CompressibleEulerEquations2D(interia_wave_setup.gamma)

boundary_conditions = (x_neg = boundary_condition_periodic,
                       x_pos = boundary_condition_periodic,
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

polydeg = 4
basis = LobattoLegendreBasis(polydeg)

surface_flux = FluxLMARS(340.0) 

volume_flux = flux_kennedy_gruber
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)
 
coordinates_min = (0.0, 0.0)
coordinates_max = (300_000.0, 10_000.0)

cells_per_dimension = (300, 10)
mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                      periodicity = (true, false))

semi = SemidiscretizationHyperbolic(mesh, equations, interia_wave_setup, solver,
                                    source_terms = interia_wave_setup,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 3000.0)  # 3000 seconds final time

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