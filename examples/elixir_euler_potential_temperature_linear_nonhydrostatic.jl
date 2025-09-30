# References:
# - Dale R. Durran, Joseph B. Klemp (1983)
#   A Compressible Model for the Simulation of Moist Mountain Waves
#   Monthly Weather Review, 111(12), 2341–2361
#   https://doi.org/10.1175/1520-0493(1983)111
#
# - F.X. Giraldo, M. Restelli (2008)
#   A study of spectral element and discontinuous Galerkin methods for the Navier–Stokes equations in nonhydrostatic mesoscale atmospheric modeling: Equation sets and test cases
#   Journal of Computational Physics, 227(8), 3849–3877
#   https://doi.org/10.1016/j.jcp.2007.12.009

using OrdinaryDiffEqSSPRK
using Trixi, TrixiAtmo

struct NonHydrostaticSetup
    theta_0::Float64  # background potential temperature in K
    u0::Float64       # background horizontal velocity in m/s
    Nf::Float64       # Brunt–väisälä frequency in 1/s
    z_B::Float64      # start of vertical damping layer in m
    z_T::Float64      # end of vertical damping layer in m
    alpha::Float64    # Rayleigh damping coefficient
    xr_B::Float64     # horizontal extent where damping starts in m (symmetric with respect to y-axis)
    function NonHydrostaticSetup(alpha, xr_B; Nf = 0.01, theta_0 = 280, u0 = 10,
                                 z_B = 15000, z_T = 30000)
        new(theta_0, u0, Nf, z_B, z_T, alpha, xr_B)
    end
end

@inline function rayleigh_damping(x, z_B, z_T, alpha, xr_B)
    xr_T = 72000

    if x[2] <= z_B
        S_v = 0.0
    else
        S_v = -alpha * sinpi(0.5f0 * (x[2] - z_B) / (z_T - z_B))^2
    end
    if x[1] < xr_B
        S_h1 = 0.0
    else
        S_h1 = -alpha * sinpi(0.5f0 * (x[1] - xr_B) / (xr_T - xr_B))^2
    end

    if x[1] > -xr_B
        S_h2 = 0.0
    else
        S_h2 = -alpha * sinpi(0.5f0 * (x[1] + xr_B) / (-xr_T + xr_B))^2
    end
    return S_v, S_h1, S_h2
end

# This signature is used for source terms, adding Rayleigh damping
# to avoid reflections at the boundaries.
@inline function (setup::NonHydrostaticSetup)(u, x, t,
                                              equations::CompressibleEulerPotentialTemperatureEquationsWithGravity2D)
    @unpack theta_0, z_B, z_T, Nf, u0, alpha, xr_B = setup

    rho, rho_v1, rho_v2, rho_theta, _ = u

    g = equations.g
    v1 = rho_v1 / rho
    theta = rho_theta / rho

    S_v, S_h1, S_h2 = rayleigh_damping(x, z_B, z_T, alpha, xr_B)

    theta_b = theta_0 * exp(Nf^2 / g * x[2])
    du2 = rho * (v1 - u0) * (S_v + S_h1 + S_h2)
    du3 = rho_v2 * (S_v + S_h1 + S_h2)
    du4 = rho * (theta - theta_b) * (S_v + S_h1 + S_h2)

    return SVector(zero(eltype(u)), du2, du3, du4, zero(eltype(u)))
end

# This signature is used for the initial condition.
@inline function (setup::NonHydrostaticSetup)(x, t,
                                              equations::CompressibleEulerPotentialTemperatureEquationsWithGravity2D)
    @unpack theta_0, u0, Nf = setup
    g = equations.g
    # Exner pressure, solves hydrostatic equation for x[2]
    exner = 1 + g^2 / (equations.c_p * theta_0 * Nf^2) * (exp(-Nf^2 / g * x[2]) - 1)
    # pressure
    p_0 = 100_000.0  # reference pressure
    R = equations.c_p - equations.c_v    # gas constant (dry air)
    p = p_0 * exner^(equations.c_p / R)
    potential_temperature = theta_0 * exp(Nf^2 / g * x[2])
    T = potential_temperature * exner
    # density
    rho = p / (R * T)
    v1 = u0
    v2 = 0.0

    return prim2cons(SVector(rho, v1, v2, p, g * x[2]), equations)
end

equations = CompressibleEulerPotentialTemperatureEquationsWithGravity2D()
alpha = 0.03
xr_B = 40000.0

linear_hydrostatic_setup = NonHydrostaticSetup(alpha, xr_B)

boundary = BoundaryConditionDirichlet(linear_hydrostatic_setup)

boundary_conditions = Dict(:x_neg => boundary,
                           :x_pos => boundary,
                           :y_neg => boundary_condition_slip_wall,
                           :y_pos => boundary)
polydeg = 3
basis = LobattoLegendreBasis(polydeg)
surface_flux = (FluxLMARS(340.0), flux_zero)
volume_flux = (flux_etec, flux_nonconservative_artiano_etal)
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

a = 1000.0
L = 144000.0
H = 30000.0
peak = 1.0
y_b = peak / (1 + (L / 2 / a)^2)
alpha = (H - y_b) * 0.5

f1(s) = SVector(-L / 2, y_b + alpha * (s + 1))
f2(s) = SVector(L / 2, y_b + alpha * (s + 1))
f3(s) = SVector((s + 1 - 1) * L / 2, peak / (1 + ((s + 1 - 1) * L / 2)^2 / a^2))
f4(s) = SVector((s + 1 - 1) * L / 2, H)
cells_per_dimension = (200, 50)
cells_per_dimension = (20, 12)
mesh = P4estMesh(cells_per_dimension, polydeg = polydeg,
                 faces = (f1, f2, f3, f4),
                 initial_refinement_level = 0, periodicity = (false, false))

semi = SemidiscretizationHyperbolic(mesh, equations, linear_hydrostatic_setup, solver,
                                    source_terms = linear_hydrostatic_setup,
                                    boundary_conditions = boundary_conditions)
T = 8
###############################################################################
# ODE solvers, callbacks etc.
tspan = (0.0, T * 3600.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000

analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback)

###############################################################################
# run the simulation
sol = solve(ode,
            SSPRK43(thread = Trixi.True());
            maxiters = 1.0e7, ode_default_options()..., callback = callbacks)
