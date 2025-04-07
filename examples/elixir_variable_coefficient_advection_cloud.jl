using OrdinaryDiffEq, Trixi
using TrixiAtmo

##############################################################################################
# initial condition

equations = VariableCoefficientAdvectionEquation2D()

# initial condition, round density cloud is defined
@inline function initial_condition_schaer_mountain_cloud(x, t, equations)
    RealT = eltype(x)
    x_0, z_0 = -50000.0f0, 9000.0f0
    rho_0 = 1.0f0
    A_x, A_z = 25000.0f0, 3000.0f0

    r = sqrt(((x[1] - x_0) / A_x)^2 + ((x[2] - z_0) / A_z)^2)

    if r <= 1
        rho = convert(RealT, rho_0 * (cos((pi * r) / 2))^2)
    else
        rho = 0.0f0
    end

    return SVector(rho)
end

# wind profile in horizontal direktion
@inline function velocity_schaer_mountain(x)
    RealT = eltype(x)
    u_0 = 10.0f0
    z_1 = 4000.0f0
    z_2 = 5000.0f0

    if x[2] <= z_1
        u_1 = 0.0f0
    elseif z_2 <= x[2]
        u_1 = 1.0f0
    else
        u_1 = convert(RealT, u_0 * (sin((pi / 2) * (x[2] - z_1) / (z_2 - z_1))))^2
    end
    return SVector(u_1, 0.0f0)
end

##############################################################################################
# semidiscretization

polydeg = 3

# P4est HOHQ mesh 
mesh_file = joinpath("src/meshes", "schaer_mountain_advection.inp")
mesh = P4estMesh{2}(mesh_file, polydeg = polydeg)

initial_condition = initial_condition_schaer_mountain_cloud

# flux and solver 
surface_flux = flux_central 

solver = DGSEM(polydeg = polydeg,
               surface_flux = surface_flux,
               volume_integral = VolumeIntegralWeakForm())
               
# boundary conditions 

boundary_conditions_dirichlet = Dict(:left => BoundaryConditionDirichletAux(initial_condition_schaer_mountain_cloud, velocity_schaer_mountain),
                                     :right => BoundaryConditionDirichletAux(initial_condition_schaer_mountain_cloud, velocity_schaer_mountain),
                                     :bottom => BoundaryConditionDirichletAux(initial_condition_schaer_mountain_cloud, velocity_schaer_mountain),#boundary_condition_slip_wall,
                                     :top => BoundaryConditionDirichletAux(initial_condition_schaer_mountain_cloud, velocity_schaer_mountain), #boundary_condition_slip_wall,
                                     :bottom_left => BoundaryConditionDirichletAux(initial_condition_schaer_mountain_cloud, velocity_schaer_mountain), 
                                     :bottom_right => BoundaryConditionDirichletAux(initial_condition_schaer_mountain_cloud, velocity_schaer_mountain))
                                     #:bottom_left_connection => BoundaryConditionDirichletAux(initial_condition_schaer_mountain_cloud, velocity_schaer_mountain),
                                     #:bottom_right_connection => BoundaryConditionDirichletAux(initial_condition_schaer_mountain_cloud, velocity_schaer_mountain))


# the velocity is passed as auxiliary_field into the cache
semi = SemidiscretizationHyperbolic(mesh,
                                    equations,
                                    initial_condition,
                                    solver,
                                    #source_terms = source_term,
                                    boundary_conditions = boundary_conditions_dirichlet,
                                    auxiliary_field = velocity_schaer_mountain)

##############################################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 5000.0)

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
solution_variables = cons2prim

analysis_callback = AnalysisCallback(semi,
                                     interval = analysis_interval,
                                     extra_analysis_errors = (:entropy_conservation_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = analysis_interval,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = "out",
                                     solution_variables = solution_variables)

#stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution)
                        #stepsize_callback)

###############################################################################
# run the simulation
sol = solve(ode,
            CarpenterKennedy2N54(williamson_condition = false),
            maxiters = 1.0e7,
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false,
            callback = callbacks);

summary_callback()
