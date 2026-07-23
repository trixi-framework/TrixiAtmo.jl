###############################################################################
# DGSEM for the linear advection equation on a prismed icosahedral grid
###############################################################################

using OrdinaryDiffEqLowStorageRK, Trixi, TrixiAtmo

###############################################################################
# Spatial discretization

function initial_condition_gravity_waves(x, t,
                                         equations)
    g = 9.81
    c_p, c_v = 1004, 717
    # center of perturbation
    x_c = 100_000.0
    a = 5_000
    H = 10_000
    R = c_p - c_v    # gas constant (dry air)
    T0 = 250
    delta = g / (R * T0)
    DeltaT = 0.001
    Tb = DeltaT * sinpi(x[2] / H) * exp(-(x[1] - x_c)^2 / a^2)
    ps = 100_000  # reference pressure
    rhos = ps / (T0 * R)
    rho_b = rhos * (-Tb / T0)
    p = ps * exp(-delta * x[2])
    rho = rhos * exp(-delta * x[2]) + rho_b * exp(-0.5 * delta * x[2])
    v1 = 20
    v2 = 0
    if p <= 0 || rho <= 0
        @show p rho
        error("Non-physical initial condition: pressure and density must be positive.")
    end
    return SVector(rho, v1, v2, p)
end

function geopotential(x)
    return 9.81 * x[2]
end

initial_condition = initial_condition_gravity_waves

equations = CovariantEulerEnergyEquationsWithGravity2D(1.4,
                                      global_coordinate_system = GlobalCartesianCoordinates())

###############################################################################
# Build DG solver.

polydeg = 4

dg = DGMulti(element_type = Quad(),
             approximation_type = SBP(),
             surface_flux = flux_lax_friedrichs,
             polydeg = polydeg)

###############################################################################
# Build mesh.

coordinates_min = (0.0, 0.0)
coordinates_max = (300_000.0, 10_000.0)
cells_per_dimension = (60, 8)

mesh = DGMultiMesh(dg,cells_per_dimension;
                   periodicity = (true, false),
                   coordinates_min = coordinates_min,
                   coordinates_max = coordinates_max)


# Transform the initial condition to the proper set of conservative variables
initial_condition_transformed = transform_initial_condition(initial_condition, equations)

boundary_conditions = (; entire_boundary = boundary_condition_slip_wall)


# A semidiscretization collects data structures and functions for the spatial discretization
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_transformed, dg,
                                    metric_terms = MetricTermsCovariant(manifold = FlatManifold()),
                                    auxiliary_field = geopotential,
                                    source_terms = source_terms_gravity,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0 to T
tspan = (0.0, 1800.0)
ode = semidiscretize(semi, tspan)

# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation 
# setup and resets the timers
summary_callback = SummaryCallback()

# The AnalysisCallback allows to analyse the solution in regular intervals and prints the 
# results
analysis_callback = AnalysisCallback(semi, interval = 100,
                                     save_analysis = true,
                                     extra_analysis_errors = (:conservation_error,),
                                     uEltype = real(dg))

# The SaveSolutionCallback allows to save the solution to a file in regular intervals
save_solution = SaveSolutionCallback(interval = 100,
                                     solution_variables = contravariant_cons2global_prim)

# The StepsizeCallback handles the re-calculation of the maximum Δt after each time step
stepsize_callback = StepsizeCallback(cfl = 0.7)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE 
# solver
callbacks = CallbackSet(summary_callback, analysis_callback, save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

# OrdinaryDiffEq's `solve` method evolves the solution in time and executes the passed 
# callbacks
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 1.0, save_everystep = false, callback = callbacks)