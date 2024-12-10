using OrdinaryDiffEq, Trixi

equations = CompressibleEulerEquations2D(1.4) # set gas gamma = 1.4

# freestream flow state with Ma_inf = 0.3
@inline function uniform_flow_state(x, t, equations::CompressibleEulerEquations2D)

    # set the freestream flow parameters
    rho_freestream = 1.0
    u_freestream = 0.3
    p_freestream = inv(equations.gamma)

    theta = 0.0 # zero angle of attack
    si, co = sincos(theta)
    v1 = u_freestream * co
    v2 = u_freestream * si

    prim = SVector(rho_freestream, v1, v2, p_freestream)
    return prim2cons(prim, equations)
end

# initial condition
initial_condition = uniform_flow_state

# boundary condition types
boundary_condition_uniform_flow = BoundaryConditionDirichlet(uniform_flow_state)

# boundary condition dictionary
boundary_conditions = Dict(:Bottom => boundary_condition_uniform_flow,
                           :Top => boundary_condition_uniform_flow,
                           :Right => boundary_condition_uniform_flow,
                           :Left => boundary_condition_uniform_flow,
                           :LeftSlant => boundary_condition_slip_wall,
                           :RightSlant => boundary_condition_slip_wall,
                           :IceCream => boundary_condition_slip_wall);

# DGSEM solver.
#    1) polydeg must be >= the polynomial order set in the HOHQMesh control file to guarantee
#       freestream preservation. As a extra task try setting polydeg=3
#    2) VolumeIntegralFluxDifferencing with central volume flux is activated
#       for dealiasing
volume_flux = flux_ranocha
solver = DGSEM(polydeg = 4, surface_flux = flux_hll,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

# create the unstructured mesh from your mesh file
mesh_file = joinpath("out", "ice_cream_straight_sides.mesh")
mesh = UnstructuredMesh2D(mesh_file)

# Create semidiscretization with all spatial discretization-related components
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_conditions)

# Create ODE problem from semidiscretization with time span from 0.0 to 2.0
tspan = (0.0, 2.0)
ode = semidiscretize(semi, tspan)

# Create the callbacks to output solution files and adapt the time step
summary_callback = SummaryCallback()
save_solution = SaveSolutionCallback(interval = 10,
                                     save_initial_solution = true,
                                     save_final_solution = true)
stepsize_callback = StepsizeCallback(cfl = 1.0)

amr_indicator = IndicatorLöhner(semi, variable=density)

amr_controller = ControllerThreeLevel(semi, amr_indicator,
                                      base_level=0,
                                      med_level =1, med_threshold=0.05,
                                      max_level =3, max_threshold=0.1)

amr_callback = AMRCallback(semi, amr_controller,
                           interval=5,
                           adapt_initial_condition=true,
                           adapt_initial_condition_only_refine=true)

callbacks = CallbackSet(summary_callback, save_solution, stepsize_callback, amr_callback)

    # Evolve ODE problem in time using `solve` from OrdinaryDiffEq
    sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
                dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
                save_everystep = false, callback = callbacks)
    # print the timer summary
    summary_callback()