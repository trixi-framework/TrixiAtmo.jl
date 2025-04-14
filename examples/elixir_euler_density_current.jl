using Trixi
using OrdinaryDiffEq

function initial_condition_density_current(x, t,
                                           equations::CompressibleEulerEquations2D)
    g = 9.81
    c_p = 1004.0
    c_v = 717.0
    # center of perturbation
    center_x = -2000.0
    center_z = 3000.0
    # radius of perturbation
    radius = 1.0
    x_r = 4000.0
    z_r = 2000.0

    # distance of current x to center of perturbation
    r = sqrt((x[1] - center_x)^2 / x_r^2 + (x[2] - center_z)^2 / z_r^2)

    # perturbation in potential temperature
    potential_temperature_ref = 300.0
    potential_temperature_perturbation = 0.0
    if r <= radius
        potential_temperature_perturbation = -15 / 2 * (1 + cospi(r))
    end
    potential_temperature = potential_temperature_ref + potential_temperature_perturbation

    # Exner pressure, solves hydrostatic equation for x[2]
    exner = 1 - g / (c_p * potential_temperature_ref) * x[2]

    # pressure
    p_0 = 100_000.0  # reference pressure
    R = c_p - c_v    # gas constant (dry air)
    p = p_0 * exner^(c_p / R)
    T = potential_temperature * exner

    # density
    rho = p / (R * T)
    v1 = 20.0
    v2 = 0.0

    return prim2cons(SVector(rho, v1, v2, p), equations)
end

@inline function source_terms_gravity(u, x, t,
                                      equations::CompressibleEulerEquations2D)
    rho, _, rho_v2, _ = u
    return SVector(zero(eltype(u)), zero(eltype(u)), -9.81 * rho,
                   -9.81 * rho_v2)
end

equations = CompressibleEulerEquations2D(1.4)

polydeg = 4

surface_flux = FluxLMARS(340.0)

volume_flux = flux_kennedy_gruber

solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

## Agnesi Profile
L = 25_600  # Length of the box domain
H = 6_400   # Height of the box domain
peak = 1000 # Peak of the mountain
a = 1000    # "width" of the mountain
y_b = peak / (1 + (L / 2 / a)^2)
alfa = (H - y_b) * 0.5
f1(s) = SVector(-L / 2, y_b + alfa * (s + 1))
f2(s) = SVector(L / 2, y_b + alfa * (s + 1))
f3(s) = SVector((s + 1 - 1) * L / 2, peak / (1 + ((s + 1 - 1) * L / 2)^2 / a^2))
f4(s) = SVector((s + 1 - 1) * L / 2, H)

trees_per_dimension = (32, 32)
mesh = P4estMesh(trees_per_dimension, polydeg = 3,
                 faces = (f1, f2, f3, f4),
                 initial_refinement_level = 1, periodicity = (true, false))

boundary_conditions = Dict(:y_neg => boundary_condition_slip_wall,
                           :y_pos => boundary_condition_slip_wall)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_density_current,
                                    solver, source_terms = source_terms_gravity,
                                    boundary_conditions = boundary_conditions)

tspan = (0.0, 900.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback, visualization_callback)

sol = solve(ode,
            SSPRK43(thread = OrdinaryDiffEq.True()),
            maxiters = 1.0e7,
            dt = 1e-1,
            save_everystep = false, callback = callbacks)