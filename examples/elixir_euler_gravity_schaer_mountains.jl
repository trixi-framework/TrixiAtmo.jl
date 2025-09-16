using OrdinaryDiffEq
using Trixi

#schaer mountain test case 

# Initial condition
function initial_condition_schaer_mountain(x, t,
                                           equations::CompressibleEulerEquationsWithGravity2D)
    g = 9.81
    c_p = 1004.0
    c_v = 717.0
    gamma = c_p / c_v

    # Exner pressure from hydrostatic balance for x[2]
    potential_temperature_int = 280.0 #constant of integration 
    bvfrequency = 0.01 #Brunt-Väisälä frequency

    exner = 1 +
            g^2 / (c_p * potential_temperature_int * bvfrequency^2) *
            (exp(-bvfrequency^2 / g * x[2]) - 1)

    # mean potential temperature
    potential_temperature = potential_temperature_int * exp(bvfrequency^2 / g * x[2])

    # temperature
    T = potential_temperature * exner

    # pressure
    p_0 = 100_000.0  # reference pressure
    R = c_p - c_v    # gas constant (dry air)
    p = p_0 * exner^(c_p / R)

    # density
    rho = p / (R * T)

    #Geopotential
    phi = g * x[2]

    v1 = 10.0
    v2 = 0.0
    E = c_v * T + 0.5 * (v1^2 + v2^2) + phi
    return SVector(rho, rho * v1, rho * v2, rho * E, phi)
end

###############################################################################

function mapping(xi_, eta_)
    xi = xi_ * 25_000  #xi_ * 10_000.0 
    eta = eta_ * 10_500 + 10_500# eta_ * 500.0 + 500.0 
    #upper boundary 
    H = 21_000.0

    #topography
    h_c = 250.0
    lambda_c = 4000.0
    a_c = 5000.0

    topo = -h_c * exp(-(xi / a_c)^2) * cos(pi * xi / lambda_c)^2

    x = xi
    y = H * (eta - topo) / (H - topo)
    return SVector(x, y)
end

# Create curved mesh with 200 x 100 elements
polydeg = 3
cells_per_dimension = (60, 30)
mesh = P4estMesh(cells_per_dimension; polydeg = polydeg, mapping = mapping,
                 periodicity = false)

###############################################################################
# semidiscretization of the compressible Euler equations
equations = CompressibleEulerEquationsWithGravity2D(1004.0 / 717.0)

initial_condition = initial_condition_schaer_mountain

boundary_conditions_dirichlet = Dict(:x_neg => BoundaryConditionDirichlet(initial_condition_schaer_mountain),
                                     :x_pos => BoundaryConditionDirichlet(initial_condition_schaer_mountain),
                                     :y_neg => boundary_condition_slip_wall,
                                     :y_pos => boundary_condition_slip_wall)

basis = LobattoLegendreBasis(polydeg)

volume_flux = (flux_kennedy_gruber, flux_nonconservative_waruszewski)
surface_flux = (FluxLMARS(340.0), flux_nonconservative_waruszewski)

volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

coordinates_min = (-25_000.0, 0.0)
coordinates_max = (25_000.0, 21_000.0)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_conditions_dirichlet)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 60 * 60)# * 10)  # 10h = 36000 s

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
solution_variables = cons2prim

analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:entropy_conservation_error,))

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
