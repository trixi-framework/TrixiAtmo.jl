using OrdinaryDiffEq, Trixi
using TrixiAtmo

##############################################################################################
# initial condition

equations = VariableCoefficientAdvectionEquation2D()

# initial condition, round density cloud is defined
@inline function initial_condition_schaer_mountain_cloud(x, t, equations)
    RealT = eltype(x)
    x_0, z_0 = -50000.0f0, 9000.0f0 #x_0, z_0 = -50000.0f0, 9000.0f0 #center location of the cloud at t=0
    rho_0 = 1.0f0
    A_x, A_z = 25000.0f0, 3000.0f0 #halfwidths in both directions

    r = sqrt(((x[1] - x_0) / A_x)^2 + ((x[2] - z_0) / A_z)^2)

    if r <= 1
        rho = convert(RealT, rho_0 * (cos((pi * r) / 2))^2)
    else
        rho = 0.0f0
    end

    return SVector(rho)
end

# wind profile in horizontal direction
@inline function velocity_schaer_mountain(x)
    RealT = eltype(x)
    u_0 = 10.0f0
    z_1 = 4000.0f0#3000.0f0#4000.0f0
    z_2 = 5000.0f0#6000.0f0#5000.0f0

    if x[2] <= z_1
        u_1 = 0.0f0
    elseif z_2 <= x[2]
        u_1 = 1.0f0
    else
        u_1 = convert(RealT, (sin((pi / 2) * (x[2] - z_1) / (z_2 - z_1))))^2
    end
    return SVector(u_0 *u_1, 0.0f0)
end

# constant wind profile for testing 
@inline function velocity_schaer_mountain_constant(x)
    u_0 = 10.0f0
    return SVector(u_0, 0.0f0)
end


##############################################################################################
# semidiscretization

polydeg = 2


initial_condition = initial_condition_schaer_mountain_cloud

# P4est transformation mesh
function mapping(xi_, eta_) 
    # transformation from [-1,1]x[-1,1] to [-75000,75000]x[0,15000]
    xi = xi_ * 75_000  
    eta = eta_ * 7_500 + 7_500
    
    H = 18_000.0 #upper boundary 
    #Y_cut = 15_000.0 #cut-off height for y

    #topography
    h_0 = 3_000.0
    lambda_c = 8_000.0
    a = 25_000.0

    topo_ = -h_0 * cos((pi * xi/(2 * a)))^2 * cos(pi * xi / lambda_c)^2

    if -a <= xi <= a
        topo = topo_
    else 
        topo = 0
    end 

    x = xi
    y = H * (eta - topo) / (H - topo)

    return SVector(x, y)
end


function mapping2(xi_, eta_) #try to refine the grid in the velocity transition area  
     # transformation from [-1,1]x[-1,1] to [-75000,75000]x[0,15000]
    xi = xi_ * 75_000  
    eta = eta_ * 7_500 + 7_500
    
    H = 18_000.0 #upper boundary 
    #Y_cut = 15_000.0 #cut-off height for y

    #topography
    h_c = 3_000.0
    lambda_c = 8_000.0
    a_c = 9_511.0

    topo = -h_c * cos(-(xi * x/ a_c)^2) * cos(pi * xi / lambda_c)^2

    x = xi
    y_ = H * (eta - topo) / (H - topo)

    z_1 = 4000.0f0#3000.0f0#4000.0f0
    z_2 = 5000.0f0#6000.0f0#5000.0f0

    if y_ < z_1
        y = y_
    elseif y_ < z_2
        # compress [z1, z2] by factor 2
        y = z_1 + (y_ - z_1) / 5
    else
        # above z2: shift upward to preserve continuity
        y = z_1 + (z_2 - z_1)/5 + (y_ - z_2)
    end

    return SVector(x, y)
end



cells_per_dimension = (50,10)
mesh_transform = P4estMesh(cells_per_dimension; polydeg=2, mapping)

mesh = mesh_transform
# flux and solver 
surface_flux = flux_lax_friedrichs 

solver = DGSEM(polydeg = polydeg,
               surface_flux = surface_flux,
               volume_integral = VolumeIntegralWeakForm())
               
# boundary conditions 
boundary_conditions_dirichlet_transform = (x_neg = BoundaryConditionDirichletAux(initial_condition_schaer_mountain_cloud, velocity_schaer_mountain),
                       x_pos = BoundaryConditionDirichletAux(initial_condition_schaer_mountain_cloud, velocity_schaer_mountain),
                       y_neg = BoundaryConditionDirichletAux(initial_condition_schaer_mountain_cloud, velocity_schaer_mountain),
                       y_pos = BoundaryConditionDirichletAux(initial_condition_schaer_mountain_cloud, velocity_schaer_mountain))

# the velocity is passed as auxiliary_field into the cache

semi = SemidiscretizationHyperbolic(mesh_transform,
                                    equations,
                                    initial_condition,
                                    solver,
                                    #source_terms = source_term,
                                    boundary_conditions = boundary_conditions_dirichlet_transform,
                                    aux_field = velocity_schaer_mountain)

##############################################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 10000.0) #10000.0

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
solution_variables = cons2prim#cons2prim_and_aux

analysis_callback = AnalysisCallback(semi,
                                     interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 10,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = "out_transform",
                                     solution_variables = solution_variables)

stepsize_callback = StepsizeCallback(cfl = 0.1)

visualization = VisualizationCallback(interval = 100)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode,
            CarpenterKennedy2N54(williamson_condition = false),
            maxiters = 1.0e7,
            dt = 1, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false,
            callback = callbacks);
summary_callback()
