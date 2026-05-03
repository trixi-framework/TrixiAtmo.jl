# # Rising bubble (moist Euler + gravity)
#
# This tutorial sets up a warm rising-bubble test case in `TrixiAtmo.jl`.
#
# The setup follows the warm-bubble benchmark described by Bryan and Fritsch,
# which includes both dry and moist atmospheric variants:
#
# Bryan, G. H., and Fritsch, J. M. (2002). A benchmark simulation for moist
# nonhydrostatic numerical models. Monthly Weather Review, 130(12), 2917-2928.
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

# ```julia
# using OrdinaryDiffEqLowStorageRK
# using Trixi
# using TrixiAtmo
# using TrixiAtmo: cons2drypot
# ```

using OrdinaryDiffEqLowStorageRK # hide
using Trixi # hide
using TrixiAtmo # hide
using TrixiAtmo: cons2drypot # hide

# ## Governing equations
#
# Here we use the two-dimensional compressible moist Euler equations with
# gravity. The parameters below specify the specific heat values in
# $\mathrm{J}\,\mathrm{kg}^{-1}\,\mathrm{K}^{-1}$ for dry air and water vapor.

c_pd = 1004 # specific heat at constant pressure for dry air
c_vd = 717  # specific heat at constant volume for dry air
c_pv = 1885 # specific heat at constant pressure for water vapor
c_vv = 1424 # specific heat at constant volume for water vapor

equations = CompressibleMoistEulerEquations2D(c_pd = c_pd, c_vd = c_vd,
                                              c_pv = c_pv, c_vv = c_vv,
                                              gravity = EARTH_GRAVITATIONAL_ACCELERATION)

equations

# ## Initial condition: warm bubble
#
# This follows the warm bubble setup from Bryan and Fritsch (2002).
#
# The background potential temperature is
# $\theta_{ref} = 300\,\mathrm{K}$, and a localized perturbation is added inside
# a bubble centered at $(x_c, z_c) = (10000\,\mathrm{m}, 2000\,\mathrm{m})$
# with radius $r_c = 2000\,\mathrm{m}$.
#
# Inside the bubble, the perturbation is
#
# $$\Delta \theta = 2\,\mathrm{K}
# \cos^2\!\left(\frac{\pi}{2}\frac{r}{r_c}\right),$$
#
# and outside the bubble it is zero.
#
# The pressure is chosen from a hydrostatic neutral background state. Then
# density and temperature are reconstructed from pressure and potential
# temperature.
#
# The moisture variables are initialized to zero.
#
# For reusable examples, this initial condition could be moved to the equations
# source file. Here we keep it local so that the tutorial is self contained and
# readers can see exactly how the warm bubble state is constructed.

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
# The mesh is periodic in the horizontal direction. Therefore, no boundary
# condition is needed in the x-direction. The vertical direction is not periodic,
# so we prescribe slip-wall boundary conditions at the bottom and top boundaries.
#
# This must match the mesh periodicity setting below:
# 
# ```julia
# periodicity = (true, false)
# ```
#
# where the first entry corresponds to the horizontal x-direction and the second
# entry corresponds to the vertical z-direction.

boundary_conditions = (; y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

# Gravity acts in the vertical direction. For the conservative variables
#
#   u = (rho, rho_v1, rho_v2, rho_E, rho_qv, rho_qc),
#
# the gravitational source contributes to the vertical momentum equation and
# to the total energy equation through the vertical velocity.

function source_terms_gravity(u, x, t, equations::CompressibleMoistEulerEquations2D)
    @unpack g = equations

    rho = u[1]
    rho_v2 = u[3]

    return SVector(zero(eltype(u)),
                   zero(eltype(u)),
                   -g * rho,
                   -g * rho_v2,
                   zero(eltype(u)),
                   zero(eltype(u)))
end

source_term = source_terms_gravity

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
#```math 
#   x \in [0, 20000]\,\mathrm{m}, \quad z \in [-5000, 15000]\,\mathrm{m}
#```
#
# We use a `P4estMesh` with horizontal periodicity and non-periodic vertical
# boundaries.
#
# The initial coarse p4est mesh has 2 x 2 macro cells. During construction, each
# macro cell is uniformly refined `initial_refinement_level` times, so the actual
# computational mesh is much finer than the initial 2 x 2 layout.

coordinates_min = (0.0, -5000.0)
coordinates_max = (20000.0, 15000.0)

trees_per_dimension = (2, 2)
initial_refinement_level = 4

mesh = P4estMesh(trees_per_dimension, polydeg = 1,
                 coordinates_min = coordinates_min,
                 coordinates_max = coordinates_max,
                 initial_refinement_level = initial_refinement_level,
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
# The problem is evolved on the shortened time interval `tspan = (0.0, 400.0)`
# so that the tutorial remains inexpensive to run. For a longer benchmark-style
# simulation, increase the final time to `1000.0` and optionally increase
# `initial_refinement_level`.

tspan = (0.0, 400.0)
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

# The entropy conservation error is a standard diagnostic for entropy-conservative
# or entropy-stable discretizations. It is useful here because the volume flux is
# the Chandrashekar flux in flux-differencing form.
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:entropy_conservation_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 1000,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = solution_variables)

# The time step is selected dynamically from a CFL condition. A moderate CFL
# number keeps the tutorial run stable while avoiding an unnecessarily small
# time step.
stepsize_callback = StepsizeCallback(cfl = 0.5)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        stepsize_callback)

callbacks

# ## Run the simulation
#
# We now run the simulation. The value `dt = 1.0` is only an initial time-step
# guess; the actual time step is adjusted by `StepsizeCallback`.

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            maxiters = 1.0e7,
            dt = 1.0,
            save_everystep = false,
            callback = callbacks)

summary_callback()

# ## Expected behavior
#
# Since the perturbation is warmer than the surrounding air, it is less dense
# and should rise under gravity. As the solution develops, the bubble should
# move upward and begin to deform.
#
# This kind of test is a standard atmospheric benchmark because it checks how
# well the method handles buoyancy, gravity, compressibility, and the spatial
# discretization all at once.
