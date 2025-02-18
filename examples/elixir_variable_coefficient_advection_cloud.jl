using OrdinaryDiffEq, Trixi, TrixiAtmo

##############################################################################################
# initial condition

equations = VariableCoefficientAdvectionEquation2D()

function initial_condition_schaer_mountain_cloud(x, t, equations)
    x_0, z_0 = -50_000.0, 9_000.0
    rho_0 = 1.0
    A_x, A_z = 25_000.0, 3_000.0

    r = sqrt(((x[1] - x_0) / A_x)^2 + ((x[2] - z_0) / A_z)^2)

    if r <= 1
        rho = rho_0 * (cos((pi*r)/2))^2
    else 
        rho = 0
    end 
    
    return SVector(rho)
end

##############################################################################################
# semidiscretization

polydeg = 3

# mesh 
mesh_file = joinpath("mesh", "schaer_mountain_1000.inp")
mesh = P4estMesh{2}(mesh_file, polydeg = polydeg)



initial_condition = initial_condition_schaer_mountain_cloud

# flux and solver 
surface_flux = FluxLMARS(340.0)
volume_flux = flux_kennedy_gruber

solver = DGSEM(
    polydeg = polydeg,
    surface_flux = surface_flux,
    volume_integral = VolumeIntegralFluxDifferencing(volume_flux),
)

# source terms 
@inline function source(u, x, t, equations::VariableCoefficientAdvectionEquation2D)
    return (source_terms_gravity(u, equations::VariableCoefficientAdvectionEquation2D))
end

source_term = source

# boundary conditions 

boundary_conditions_dirichlet = Dict(
    :left => BoundaryConditionDirichlet(initial_condition_schaer_mountain_cloud),
    :right => BoundaryConditionDirichlet(initial_condition_schaer_mountain_cloud),
    :bottom => boundary_condition_slip_wall,
    :top => boundary_condition_slip_wall,
)


semi = SemidiscretizationHyperbolic(
    mesh,
    equations,
    initial_condition,
    solver,
    source_terms = source_term,
    boundary_conditions = boundary_conditions_dirichlet,
    auxiliary_field = velocity_schaer_mountain,
)

##############################################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 5000.0)

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
solution_variables = cons2prim

analysis_callback = AnalysisCallback(
    semi,
    interval = analysis_interval,
    extra_analysis_errors = (:entropy_conservation_error,),
)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(
    interval = analysis_interval,
    save_initial_solution = true,
    save_final_solution = true,
    output_directory = "out",
    solution_variables = solution_variables,
)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(
    summary_callback,
    analysis_callback,
    alive_callback,
    save_solution,
    stepsize_callback,
)

###############################################################################
# run the simulation
sol = solve(
    ode,
    CarpenterKennedy2N54(williamson_condition = false),
    maxiters = 1.0e7,
    dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
    save_everystep = false,
    callback = callbacks,
);

summary_callback()
