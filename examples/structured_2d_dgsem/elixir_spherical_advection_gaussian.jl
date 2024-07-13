###############################################################################
# DGSEM for the linear advection equation on the cubed sphere
###############################################################################

using OrdinaryDiffEq, Trixi, TrixiAtmo
using Rotations  # to compute solid body rotation

###############################################################################
# Problem definition

function initial_condition_gaussian_bell(x, t, equations)
    
    (; manifold) = equations
    RealT = eltype(x)

    # Get Earth's radius
    (; a) = manifold

    # speed of propagation
    V = convert(RealT, 2π) * a / (12 * SECONDS_PER_DAY)

    # axis of rotation
    α = convert(RealT, π/4)

    # initial height of the bump
    h_ref = 1000f0

    # initial location of the bump
    (λ_0, θ_0) = (convert(RealT, 3π/2), convert(RealT, 0f0))

    # width of the bump
    b_0 = 5f0 / (a^2)

    # get zonal and meridional components of the velocity
    λ, θ = face2global(x, manifold) 
    v_λ = V * (cos(θ)*cos(α) + sin(θ)*cos(λ)*sin(α))
    v_θ = -V * sin(λ)*sin(α)

    # compute exact solution through solid body rotation
    x_cartesian = global2cartesian(λ, θ, manifold)
    x_0_cartesian = SVector(global2cartesian(λ_0, θ_0, manifold))
    x_0_rotated = AngleAxis(V*t/a, -sin(α), 0.0, cos(α)) * x_0_cartesian

    # evaluate the Gaussian bump function
    h = h_ref * exp(-b_0 * ((x_cartesian[1] - x_0_rotated[1])^2 +
                        (x_cartesian[2] - x_0_rotated[2])^2 +
                        (x_cartesian[3] - x_0_rotated[3])^2))

    return global2cons(h, v_λ, v_θ, λ, θ, equations)
end

# End time of simulation
days = 12
T = days * SECONDS_PER_DAY

###############################################################################
# Spatial discretization

p = 3
cells_per_dimension = 2

# Create DG solver with polynomial degree = p and a local Lax-Friedrichs flux
solver = DGSEM(polydeg = p, surface_flux = flux_lax_friedrichs)

# Get geometric and topological information for a cubed sphere of radius EARTH_RADIUS
mesh = StructuredMesh((cells_per_dimension, cells_per_dimension), 
    (-π/4, -π/4), (π/4, π/4), periodicity = false)
faces = [CubedSphereFace2D{i}(EARTH_RADIUS) for i in 1:6]
connectivity = cubed_sphere_connectivity(spherical_coupling)

# Define an equation system for each face of the cubed sphere
equations = [SphericalLinearAdvectionEquation2D(face) for face in faces]

# Define a SemidiscretizationHyperbolic for each of the six faces of the cube
semis = [SemidiscretizationHyperbolic(mesh, equations[i],
    initial_condition_gaussian_bell, solver, source_terms = nothing,
    boundary_conditions = connectivity[i]) for i in eachindex(faces)]

# Create a SemidiscretizationCoupled that bundles all the semidiscretizations
semi = SemidiscretizationCoupled(semis...)

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0 to T
ode = semidiscretize(semi, (0.0, T))

summary_callback = SummaryCallback()
alive_callback = AliveCallback(alive_interval = 10)
stepsize_callback = StepsizeCallback(cfl = 0.4)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE solver
callbacks = CallbackSet(summary_callback, alive_callback, stepsize_callback)

###############################################################################
# Run the simulation

# OrdinaryDiffEq's `solve` method evolves the solution in time and executes the passed callbacks
sol = solve(ode, CarpenterKennedy2N54(), dt = 1.0, save_everystep = false, callback = callbacks)
summary_callback()  # print the timer summary