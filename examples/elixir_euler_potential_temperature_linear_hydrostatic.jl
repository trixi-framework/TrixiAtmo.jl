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

struct HydrostaticSetup
    T_0::Float64     # background temperature in K
    u0::Float64      # background horizontal velocity in m/s
    Nf::Float64      # Brunt–Väisälä frequency in 1/s
    z_B::Float64     # bottom boundary of damping layer in m
    z_T::Float64     # top boundary of damping layer in m
    alpha::Float64   # Rayleigh damping coefficient
    xr_B::Float64    # horizontal extent where damping starts in m (symmetric with respect to y-axis)
    function HydrostaticSetup(alpha, xr_B, equations; T_0 = 250.0, u0 = 20.0, z_B = 15000.0,
                              z_T = 30000.0)
        Nf = equations.g / sqrt(equations.c_p * T_0)
        new(T_0, u0, Nf, z_B, z_T, alpha, xr_B)
    end
end

@inline function rayleigh_damping(x, z_B, z_T, alpha, xr_B)
    xr_T = 120000.0

    if x[2] <= z_B
        S_v = 0.0
    else
        S_v = -alpha * sinpi(0.5 * (x[2] - z_B) / (z_T - z_B))^2
    end
    if x[1] < xr_B
        S_h1 = 0.0
    else
        S_h1 = -alpha * sinpi(0.5 * (x[1] - xr_B) / (xr_T - xr_B))^2
    end

    if x[1] > -xr_B
        S_h2 = 0.0
    else
        S_h2 = -alpha * sinpi(0.5 * (x[1] + xr_B) / (-xr_T + xr_B))^2
    end
    return S_v, S_h1, S_h2
end

# This signature is used for source terms, adding Rayleigh damping
# to avoid reflections at the boundaries.
@inline function (setup::HydrostaticSetup)(u, x, t,
                                           equations::CompressibleEulerPotentialTemperatureEquationsWithGravity2D)
    @unpack T_0, z_B, z_T, Nf, u0, alpha, xr_B = setup
    g = equations.g

    rho, rho_v1, rho_v2, rho_theta, _ = u

    v1 = rho_v1 / rho
    theta = rho_theta / rho

    S_v, S_h1, S_h2 = rayleigh_damping(x, z_B, z_T, alpha, xr_B)

    exner = exp(-Nf^2 / g * x[2])

    theta_0 = T_0 / exner

    du2 = rho * (v1 - u0) * (S_v + S_h1 + S_h2)
    du3 = rho_v2 * (S_v + S_h1 + S_h2)
    du4 = rho * (theta - theta_0) * (S_v + S_h1 + S_h2)

    return SVector(zero(eltype(u)), du2, du3, du4, zero(eltype(u)))
end

# This signature is used for the initial condition.
@inline function (setup::HydrostaticSetup)(x, t,
                                           equations::CompressibleEulerPotentialTemperatureEquationsWithGravity2D)
    @unpack T_0, u0, Nf = setup
    g = equations.g

    # Exner pressure, solves hydrostatic equation for x[2]
    exner = exp(-Nf^2 / g * x[2])
    # pressure
    p_0 = 100_000.0  # reference pressure
    R = equations.c_p - equations.c_v    # gas constant (dry air)
    p = p_0 * exner^(equations.c_p / equations.R)

    # density
    rho = p / (R * T_0)
    v1 = u0
    v2 = 0.0

    return prim2cons(SVector(rho, v1, v2, p, g * x[2]), equations)
end

equations = CompressibleEulerPotentialTemperatureEquationsWithGravity2D(1004, 717, 9.81)
alpha = 0.035
xr_B = 60000.0
linear_hydrostatic_setup = HydrostaticSetup(alpha, xr_B, equations)

boundary = BoundaryConditionDirichlet(linear_hydrostatic_setup)

boundary_conditions = Dict(:x_neg => boundary,
                           :x_pos => boundary,
                           :y_neg => boundary_condition_slip_wall,
                           :y_pos => boundary)

# We have an isothermal background state with T0 = 250 K. 
# The reference speed of sound can be computed as:
# cs = sqrt(gamma * R * T0)
cs = sqrt(equations.gamma * equations.R * linear_hydrostatic_setup.T_0)
polydeg = 3
basis = LobattoLegendreBasis(polydeg)

surface_flux = (FluxLMARS(cs), flux_zero)
volume_flux = (flux_tec, flux_nonconservative_waruzewski_etal)
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

a = 10000.0
L = 240000.0
H = 30000.0
peak = 1.0
y_b = peak / (1 + (L / 2 / a)^2)
alpha_b = (H - y_b) * 0.5

f1(s) = SVector(-L / 2, y_b + alpha_b * (s + 1))
f2(s) = SVector(L / 2, y_b + alpha_b * (s + 1))
f3(s) = SVector((s + 1 - 1) * L / 2, peak / (1 + ((s + 1 - 1) * L / 2)^2 / a^2))
f4(s) = SVector((s + 1 - 1) * L / 2, H)
cells_per_dimension = (100, 60)
mesh = P4estMesh(cells_per_dimension, polydeg = polydeg,
                 faces = (f1, f2, f3, f4),
                 initial_refinement_level = 0, periodicity = (false, false))

semi = SemidiscretizationHyperbolic(mesh, equations, linear_hydrostatic_setup, solver,
                                    source_terms = linear_hydrostatic_setup,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.
T = 12.5
tspan = (0.0, T * 3600.0)  # 12.5 hours final time
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
            maxiters = 1.0e7,
            ode_default_options()..., callback = callbacks)
