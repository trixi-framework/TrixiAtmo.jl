using OrdinaryDiffEqSSPRK
using Trixi, TrixiAtmo

function initial_condition_robert_bubble(x, t,
                                         equations::CompressibleEulerPotentialTemperatureEquations2D)
    # center of perturbation
    center_x = 500.0
    center_z = 10.0
    g = equations.g
    # radius of perturbation
    radius = 0.0
    # distance of current x to center of perturbation
    r = sqrt((x[1] - center_x)^2 + (x[2] - center_z)^2)

    # perturbation in potential temperature
    potential_temperature_ref = 300.0
    potential_temperature_perturbation = 0.0
    if r <= radius
        potential_temperature_perturbation = 5.0
    end
    potential_temperature = potential_temperature_ref + potential_temperature_perturbation

    # Exner pressure, solves hydrostatic equation for x[2]
    exner = 1 - g / (equations.c_p * potential_temperature) * x[2]

    # pressure
    p_0 = 100_000.0  # reference pressure
    R = equations.c_p - equations.c_v    # gas constant (dry air)
    p = p_0 * exner^(equations.c_p / R)

    # temperature
    T = potential_temperature * exner

    # density
    rho = p / (R * T)

    v1 = 20.0
    v2 = 0.0
    return TrixiAtmo.prim2cons(SVector(rho, v1, v2, p), equations)
end

equations = CompressibleEulerPotentialTemperatureEquations2D()

surface_flux = FluxLMARS(340)
volume_flux = flux_tec
polydeg = 3
solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

boundary_conditions = (x_neg = boundary_condition_periodic,
                       x_pos = boundary_condition_periodic,
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

coordinates_min = (0.0, 0.0)
coordinates_max = (1000.0, 1000.0)
cells_per_dimension = (32, 16)
mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                      periodicity = (true, false))

@inline function source_terms_gravity(u, x, t,
                                      equations::CompressibleEulerPotentialTemperatureEquations2D)
    rho, _, _, _ = u
    return SVector(zero(eltype(u)), zero(eltype(u)), -equations.g * rho, zero(eltype(u)))
end
source_terms = source_terms_gravity
initial_condition = initial_condition_robert_bubble
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms,
                                    boundary_conditions = boundary_conditions)
tspan = (0.0, 1800.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_integrals = (entropy,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback)

sol = solve(ode,
            SSPRK43(thread = Trixi.True()),
            maxiters = 1.0e7,
            dt = 1e-1, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks)
