# # Rising bubble (moist Euler + gravity)
#
# This tutorial sets up a warm rising-bubble test case in `TrixiAtmo.jl`.
#
# The setup is based on a standard atmospheric benchmark where a localized warm
# perturbation is placed in a background atmosphere. Since the perturbation is
# warmer than the air around it, it should rise as the solution evolves.
#
# This example uses the moist Euler equations, but the moisture variables are
# initialized to zero in our initial condition.
#
# In this tutorial we go through:
#
# - the governing equations
# - the warm-bubble initial condition
# - the mesh and DG discretization
# - the source term and boundary conditions
# - the time integration setup
#
# ## Load packages
#
# We start by loading the time integrator, the core `Trixi.jl` tools, and the
# atmospheric extensions from `TrixiAtmo.jl`.

using OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiAtmo
using TrixiAtmo: source_terms_geopotential, cons2drypot

# ## Governing equations
#
# Here we use the two-dimensional compressible moist Euler equations with
# gravity. The parameters below specify the thermodynamic constants for dry air
# and water vapor.

c_pd = 1004 # specific heat at constant pressure for dry air
c_vd = 717  # specific heat at constant volume for dry air
c_pv = 1885 # specific heat at constant pressure for moist air
c_vv = 1424 # specific heat at constant volume for moist air

equations = CompressibleMoistEulerEquations2D(c_pd = c_pd, c_vd = c_vd,
                                              c_pv = c_pv, c_vv = c_vv,
                                              gravity = EARTH_GRAVITATIONAL_ACCELERATION)

equations

# ## Initial condition: warm bubble
#
# This follows the warm bubble setup from Wicker and Skamarock (1998).
#
# The background potential temperature is
# $\theta_{ref} = 300$ K, and a localized perturbation is added inside a bubble
# centered at $(x_c, z_c) = (10000, 2000)$ with radius $r_c = 2000$.
#
# Inside the bubble, the perturbation is
#
# $$\Delta \theta = 2 \cos^2\!\left(\frac{\pi}{2}\frac{r}{r_c}\right),$$
#
# and outside the bubble it is zero.
#
# The pressure is chosen from a hydrostatic neutral background state. Then
# density and temperature are reconstructed from pressure and potential
# temperature.
#
# The moisture variables are initialized to zero.

function initial_condition_warm_bubble(x, t, equations::CompressibleMoistEulerEquations2D)
    @unpack p_0, kappa, g, c_pd, c_vd, R_d, R_v = equations

    xc = 10000
    zc = 2000
    r = sqrt((x[1] - xc)^2 + (x[2] - zc)^2)
    rc = 2000
    θ_ref = 300
    Δθ = 0

    if r <= rc
        Δθ = 2 * cospi(0.5f0 * r / rc)^2
    end

    ## Potential temperature with warm perturbation
    θ = θ_ref + Δθ

    ## Background hydrostatic pressure
    p = p_0 * (1 - kappa * g * x[2] / (R_d * θ_ref))^(c_pd / R_d)

    ## Recover density and temperature
    rho = p / ((p / p_0)^kappa * R_d * θ)
    T = p / (R_d * rho)

    v1 = 20
    v2 = 0
    rho_v1 = rho * v1
    rho_v2 = rho * v2
    rho_E = rho * c_vd * T + 1 / 2 * rho * (v1^2 + v2^2)

    return SVector(rho, rho_v1, rho_v2, rho_E,
                   zero(eltype(g)), zero(eltype(g)))
end

initial_condition = initial_condition_warm_bubble

# ## Boundary conditions and source term
#
# The mesh below is periodic in the horizontal direction. At the bottom and top
# boundaries we use slip-wall boundary conditions.
#
# Gravity is added through `source_terms_geopotential`.

boundary_conditions = (; y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

source_term = source_terms_geopotential

# ## DG discretization
#
# We use a DGSEM discretization with polynomial degree `polydeg = 4`. The
# surface flux is LMARS, and the volume flux is the Chandrashekar flux in
# flux-differencing form.

polydeg = 4

surface_flux = FluxLMARS(360.0)
volume_flux = flux_chandrashekar
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(polydeg, surface_flux, volume_integral)

solver

# ## Mesh
#
# The computational domain is
#
# - $x \in [0, 20000]$
# - $z \in [-5000, 15000]$
#
# We use a `P4estMesh` with horizontal periodicity and non-periodic vertical
# boundaries.

coordinates_min = (0.0, -5000.0)
coordinates_max = (20000.0, 15000.0)

trees_per_dimension = (2, 2)

mesh = P4estMesh(trees_per_dimension, polydeg = 1,
                 coordinates_min = coordinates_min,
                 coordinates_max = coordinates_max,
                 initial_refinement_level = 5,
                 periodicity = (true, false))

mesh

# ## Semidiscretization
#
# Next we combine the mesh, equations, initial condition, solver, boundary
# conditions, and source term into a semidiscretization.

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_conditions,
                                    source_terms = source_term)

semi

# ## Time interval and ODE problem
#
# The problem is evolved on the time interval `tspan = (0.0, 1000.0)`.

tspan = (0.0, 1000.0)
ode = semidiscretize(semi, tspan)

# ## Callbacks and diagnostics
#
# We use a standard collection of callbacks:
#
# - `SummaryCallback()` for a summary of the run
# - `AnalysisCallback(...)` for diagnostics
# - `AliveCallback(...)` for progress information
# - `SaveSolutionCallback(...)` to write output
# - `StepsizeCallback(...)` to control the time step with a CFL condition
#
# The saved variables are given by `cons2drypot`.

summary_callback = SummaryCallback()

analysis_interval = 1000
solution_variables = cons2drypot

analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:entropy_conservation_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 1000,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = solution_variables)

stepsize_callback = StepsizeCallback(cfl = 0.2)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        stepsize_callback)

callbacks

# ## Run the simulation
#
# The full simulation can be run with:
#
# ```julia
# sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
#             maxiters = 1.0e7,
#             dt = 1.0,
#             save_everystep = false,
#             callback = callbacks)
# ```
#
# The value `dt = 1.0` is only a placeholder here. In practice, the actual time
# step is adjusted by `StepsizeCallback`.

# ## Expected behavior
#
# Since the perturbation is warmer than the surrounding air, it is less dense
# and should rise under gravity. As the solution develops, the bubble should
# move upward and begin to deform.
#
# This kind of test is a standard atmospheric benchmark because it checks how
# well the method handles buoyancy, gravity, compressibility, and the spatial
# discretization all at once.