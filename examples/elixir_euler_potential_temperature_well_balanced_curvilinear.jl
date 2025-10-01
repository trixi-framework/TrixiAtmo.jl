# References:
# - Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
#   Structure-Preserving High-Order Methods for the Compressible Euler Equations 
#   in Potential Temperature Formulation for Atmospheric Flows
#   https://arxiv.org/abs/2509.10311 (pre-print)
using OrdinaryDiffEqSSPRK
using Trixi, TrixiAtmo

function initial_condition_adiabatic(x, t,
                                     equations::CompressibleEulerPotentialTemperatureEquationsWithGravity2D)
    g = equations.g
    c_p = equations.c_p
    c_v = equations.c_v
    # center of perturbation
    p0 = 100_000
    # perturbation in potential temperature
    R = c_p - c_v    # gas constant (dry air)
    potential_temperature = 300.0

    # Exner pressure, solves hydrostatic equation for x[2]
    exner = 1 - g / (c_p * potential_temperature) * x[2]

    # pressure
    R = c_p - c_v    # gas constant (dry air)
    p = p0 * exner^(c_p / R)

    # temperature
    T = potential_temperature * exner
    # density
    rho = p / (R * T)
    v1 = 0
    v2 = 0
    return prim2cons(SVector(rho, v1, v2, p, g * x[2]), equations)
end

equations = CompressibleEulerPotentialTemperatureEquationsWithGravity2D()
polydeg = 2
basis = LobattoLegendreBasis(polydeg)
cs = 340
surface_flux = (FluxLMARS(cs), flux_zero)
volume_flux = (flux_tec, flux_nonconservative_artiano_etal)
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
solver = DGSEM(basis, surface_flux, volume_integral)

function mapping(xi, eta)
    x = xi + 0.1 * sinpi(xi) * sinpi(eta)
    y = eta + 0.1 * sinpi(xi) * sinpi(eta)
    return SVector(1000 * 0.5 * (1 + x), 1000 * 0.5 * (1 + y))
end
trees_per_dimension = (16, 16)

mesh = P4estMesh(trees_per_dimension, polydeg = polydeg,
                 mapping = mapping, periodicity = (false, false),
                 initial_refinement_level = 0)

boundary_conditions = Dict(:x_pos => boundary_condition_slip_wall,
                           :x_neg => boundary_condition_slip_wall,
                           :y_pos => boundary_condition_slip_wall,
                           :y_neg => boundary_condition_slip_wall)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_adiabatic, solver,
                                    boundary_conditions = boundary_conditions)

dt = 0.01
tspan = (0, 0.4)

summary_callback = SummaryCallback()

analysis_interval = 1
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback)

ode = semidiscretize(semi, tspan)

sol = solve(ode, SSPRK43(); dt = dt, ode_default_options()..., callback = callbacks,
            adaptive = false)
