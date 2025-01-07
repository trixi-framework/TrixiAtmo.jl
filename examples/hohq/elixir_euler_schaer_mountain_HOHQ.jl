using OrdinaryDiffEq
using Trixi
using TrixiAtmo: cons2pot, source_terms_rayleigh_sponge, source_terms_rayleigh_sponge_left, source_terms_rayleigh_sponge_right, source_terms_gravity

#schaer mountain test case 


# Initial condition
function initial_condition_schaer_mountain(x, t,
                                           equations::CompressibleEulerEquations2D)
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

    v1 = 10.0
    v2 = 0.0
    E = c_v * T + 0.5 * (v1^2 + v2^2)
    return SVector(rho, rho * v1, rho * v2, rho * E)
end

# Source terms 
@inline function source(u, x, t, equations::CompressibleEulerEquations2D)
    return (source_terms_gravity(u, equations::CompressibleEulerEquations2D))
end


###############################################################################
# semidiscretization of the compressible Euler equations
equations = CompressibleEulerEquations2D(1004.0 / 717.0)

initial_condition = initial_condition_schaer_mountain
source_term = source

boundary_conditions_dirichlet = Dict(:left => BoundaryConditionDirichlet(initial_condition_schaer_mountain),
                       :right => BoundaryConditionDirichlet(initial_condition_schaer_mountain),
                       :bottom => boundary_condition_slip_wall,
                       :top => boundary_condition_slip_wall)

boundary_conditions_periodic = (x_neg = boundary_condition_periodic, 
                       x_pos = boundary_condition_periodic, 
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

polydeg = 3
basis = LobattoLegendreBasis(polydeg)

surface_flux = FluxLMARS(340.0)

volume_flux = flux_kennedy_gruber
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

mesh_file = joinpath("mesh", "schaer_mountain_250.inp")
mesh_file_unstructured = joinpath("mesh", "schaer_mountain_250.mesh")
mesh = P4estMesh{2}(mesh_file, polydeg = polydeg)
#mesh = UnstructuredMesh2D(mesh_file_unstructured)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_term,
                                    boundary_conditions = boundary_conditions_dirichlet)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 60 * 60 * 10)  # 10h = 36000 s

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