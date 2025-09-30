# References:
# - Christoph Schär, Daniel Leuenberger, Oliver Fuhrer, Daniel Lüthi, Claude Girard (2002)
#   A New Terrain-Following Vertical Coordinate Formulation for Atmospheric Prediction Models
#   Monthly Weather Review, 130(10), 2459–2480
#   https://doi.org/10.1175/1520-0493(2002)130

using OrdinaryDiffEqSSPRK
using Trixi, TrixiAtmo

struct SchärSetup
    theta_0::Float64  # background potential temperature in K
    u0::Float64       # background horizontal velocity in m/s
    Nf::Float64       # Brunt–väisälä frequency in 1/s
    z_B::Float64      # start of vertical damping layer in m
    z_T::Float64      # end of vertical damping layer in m
    alpha::Float64    # Rayleigh damping coefficient
    xr_B::Float64     # horizontal extent where damping starts in m (symmetric with respect to y-axis)
    function SchärSetup(alpha, xr_B; theta_0 = 280, u0 = 10, z_B = 13000,
                        z_T = 21000)
        Nf = 0.01
        new(theta_0, u0, Nf, z_B, z_T, alpha, xr_B)
    end
end

@inline function rayleigh_damping(x, z_B, z_T, alpha, xr_B)
    xr_T = 25000

    if x[2] <= z_B
        S_v = 0
    else
        S_v = -alpha * sinpi(0.5f0 * (x[2] - z_B) / (z_T - z_B))^2
    end
    if x[1] < xr_B
        S_h1 = 0
    else
        S_h1 = -alpha * sinpi(0.5f0 * (x[1] - xr_B) / (xr_T - xr_B))^2
    end

    if x[1] > -xr_B
        S_h2 = 0
    else
        S_h2 = -alpha * sinpi(0.5f0 * (x[1] + xr_B) / (-xr_T + xr_B))^2
    end
    return S_v, S_h1, S_h2
end

# This signature is used for source terms, adding Rayleigh damping
# to avoid reflections at the boundaries.
function (setup::SchärSetup)(u, x, t,
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
function (setup::SchärSetup)(x, t,
                             equations::CompressibleEulerPotentialTemperatureEquationsWithGravity2D)
    @unpack theta_0, u0, Nf = setup
    g = equations.g
    c_p = equations.c_p
    c_v = equations.c_v
    # Exner pressure, solves hydrostatic equation for x[2]
    exner = 1 + g^2 / (c_p * theta_0 * Nf^2) * (exp(-Nf^2 / g * x[2]) - 1)
    # pressure
    p_0 = 100_000  # reference pressure
    R = c_p - c_v    # gas constant (dry air)
    p = p_0 * exner^(c_p / R)
    potential_temperature = theta_0 * exp(Nf^2 / g * x[2])
    T = potential_temperature * exner
    # density
    rho = p / (R * T)
    v1 = u0
    v2 = 0

    return prim2cons(SVector(rho, v1, v2, p, g * x[2]), equations)
end

###############################################################################
# semidiscretization of the compressible Euler equations
equations = CompressibleEulerPotentialTemperatureEquationsWithGravity2D()
alpha = 0.03
xr_B = 20000
schär_setup = SchärSetup(alpha, xr_B)
boundary = BoundaryConditionDirichlet(schär_setup)
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

a = 5000.0
L = 50000.0
H = 21000.0
lambda_c = 4000.0
hc = 250.0
y_b = hc * exp(-(L / 2 / a)^2) * cospi(L / 2 / lambda_c)^2
alpha = (H - y_b) * 0.5

f1(s) = SVector(-L / 2, y_b + alpha * (s + 1))
f2(s) = SVector(L / 2, y_b + alpha * (s + 1))
f3(s) = SVector(s * L / 2, hc * exp(-(s * L / 2 / a)^2) * cospi(s * L / 2 / lambda_c)^2)
f4(s) = SVector(s * L / 2, H)
cells_per_dimension = (100, 50)
cells_per_dimension = (20, 12)
mesh = P4estMesh(cells_per_dimension, polydeg = polydeg,
                 faces = (f1, f2, f3, f4),
                 initial_refinement_level = 0, periodicity = (false, false))

semi = SemidiscretizationHyperbolic(mesh, equations, schär_setup, solver,
                                    source_terms = schär_setup,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.
T = 5
tspan = (0.0, T * 3600.0) #1 hour.
tspan = (0.0, 360)

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
            SSPRK43(thread = Trixi.True()),
            maxiters = 1.0e7, ode_default_options()..., callback = callbacks)
