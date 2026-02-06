using OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiAtmo

# To run a convergence test, we have two options:
# 1. Use the p4est variable initial_refinement_level to refine the grid:
#    - To do this, line 46 ("initial_refinement_level = 0") must NOT be a comment
#    - Call convergence_test("../examples/advection/cartesian/elixir_quad_icosahedron.jl", 4, initial_refinement_level = 0)
#    - NOT OPTIMAL: Good convergence the first iterations, but then it stagnates. Reason: The geometry does not improve with refinement.
# 2. Use the variable trees_per_face_dimension of P4estMeshQuadIcosahedron2D
#    - To do this, line 46 ("initial_refinement_level = 0") MUST BE commented/removed
#    - Call convergence_test("../examples/advection/cartesian/ elixir_quad_icosahedron.jl", 4, cells_per_dimension = (1,1))
#    - OPTIMAL convergence of polydeg + 1. Reason: The geometry improves with refinement.

###############################################################################
# semidiscretization of the linear advection equation
initial_condition = initial_condition_gaussian
polydeg = 3
cells_per_dimension = (2, 2)

# We use the ShallowWaterEquations3D equations structure but modify the rhs! function to
# convert it to a variable-coefficient advection equation
equations = ShallowWaterEquations3D(gravity = 0)

# Create DG solver with polynomial degree = 3 and (local) Lax-Friedrichs/Rusanov flux as surface flux
solver = DGSEM(polydeg = polydeg,
               surface_flux = (flux_lax_friedrichs, flux_nonconservative_wintermeyer_etal))

# Source term function to transform the Euler equations into a linear advection equation with variable advection velocity
function source_terms_convert_to_linear_advection(u, du, x, t,
                                                  equations::ShallowWaterEquations3D,
                                                  normal_direction)
    v1 = u[2] / u[1]
    v2 = u[3] / u[1]
    v3 = u[4] / u[1]

    s2 = du[1] * v1 - du[2]
    s3 = du[1] * v2 - du[3]
    s4 = du[1] * v3 - du[4]

    return SVector(0, s2, s3, s4, 0)
end

# Hack to use the weak form kernel with ShallowWaterEquations3D (a non-conservative equation).
# This works only because we have a constant bottom topography, so the equations are effectively
# conservative. Note that the weak form kernel is NOT equal to the flux differencing kernel
# with central fluxes because of the curved geometry!
@inline function Trixi.weak_form_kernel!(du, u,
                                         element,
                                         mesh::Union{StructuredMesh{2},
                                                     UnstructuredMesh2D,
                                                     P4estMesh{2}, T8codeMesh{2}},
                                         nonconservative_terms::Trixi.True,
                                         equations::Trixi.AbstractEquations{3},
                                         dg::DGSEM, cache, alpha = true)
    Trixi.weak_form_kernel!(du, u, element, mesh, Trixi.False(), equations, dg, cache)
end

# Create a 2D quad-based icosahedral mesh the size of the Earth
mesh = P4estMeshQuadIcosahedron2D(cells_per_dimension[1], EARTH_RADIUS,
                                  #initial_refinement_level = 0,
                                  polydeg = polydeg)

initial_condition_transformed = transform_initial_condition(initial_condition, equations)

# A semidiscretization collects data structures and functions for the spatial discretization
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_transformed, solver,
                                    source_terms = source_terms_convert_to_linear_advection)

###############################################################################
# ODE solvers, callbacks etc.

ode = semidiscretize(semi, (0.0, 12 * SECONDS_PER_DAY))

# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation setup
# and resets the timers
summary_callback = SummaryCallback()

# The AnalysisCallback allows to analyse the solution in regular intervals and prints the results
analysis_callback = AnalysisCallback(semi, interval = 100,
                                     save_analysis = true,
                                     extra_analysis_errors = (:conservation_error,))

# The SaveSolutionCallback allows to save the solution to a file in regular intervals
save_solution = SaveSolutionCallback(interval = 100,
                                     solution_variables = cons2prim_and_vorticity)

# The StepsizeCallback handles the re-calculation of the maximum Δt after each time step
stepsize_callback = StepsizeCallback(cfl = 0.7)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE solver
callbacks = CallbackSet(summary_callback, analysis_callback, save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

# OrdinaryDiffEq's `solve` method evolves the solution in time and executes the passed callbacks
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks)
