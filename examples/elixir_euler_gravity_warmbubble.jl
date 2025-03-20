
using OrdinaryDiffEqLowStorageRK
using Trixi, TrixiAtmo

# Warm bubble test case from
# - Wicker, L. J., and Skamarock, W. C. (1998)
#   A time-splitting scheme for the elastic equations incorporating
#   second-order Runge–Kutta time differencing
#   [DOI: 10.1175/1520-0493(1998)126%3C1992:ATSSFT%3E2.0.CO;2](https://doi.org/10.1175/1520-0493(1998)126%3C1992:ATSSFT%3E2.0.CO;2)
# See also
# - Bryan and Fritsch (2002)
#   A Benchmark Simulation for Moist Nonhydrostatic Numerical Models
#   [DOI: 10.1175/1520-0493(2002)130<2917:ABSFMN>2.0.CO;2](https://doi.org/10.1175/1520-0493(2002)130<2917:ABSFMN>2.0.CO;2)
# - Carpenter, Droegemeier, Woodward, Hane (1990)
#   Application of the Piecewise Parabolic Method (PPM) to
#   Meteorological Modeling
#   [DOI: 10.1175/1520-0493(1990)118<0586:AOTPPM>2.0.CO;2](https://doi.org/10.1175/1520-0493(1990)118<0586:AOTPPM>2.0.CO;2)
struct WarmBubbleSetup
    # Physical constants
    g::Float64       # gravity of earth
    c_p::Float64     # heat capacity for constant pressure (dry air)
    c_v::Float64     # heat capacity for constant volume (dry air)
    gamma::Float64   # heat capacity ratio (dry air)

    function WarmBubbleSetup(; g = 9.81, c_p = 1004.0, c_v = 717.0, gamma = c_p / c_v)
        new(g, c_p, c_v, gamma)
    end
end

# Initial condition
function (setup::WarmBubbleSetup)(x, t, equations::CompressibleEulerEquationsWithGravity2D)
    RealT = eltype(x)
    @unpack g, c_p, c_v = setup

    # center of perturbation
    center_x = 10000
    center_z = 2000
    # radius of perturbation
    radius = 2000
    # distance of current x to center of perturbation
    r = sqrt((x[1] - center_x)^2 + (x[2] - center_z)^2)

    # perturbation in potential temperature
    potential_temperature_ref = 300
    potential_temperature_perturbation = zero(RealT)
    if r <= radius
        potential_temperature_perturbation = 2 * cospi(0.5f0 * r / radius)^2
    end
    potential_temperature = potential_temperature_ref + potential_temperature_perturbation

    # Exner pressure, solves hydrostatic equation for x[2]
    exner = 1 - g / (c_p * potential_temperature) * x[2]

    # pressure
    p_0 = 100_000  # reference pressure
    R = c_p - c_v    # gas constant (dry air)
    p = p_0 * exner^(c_p / R)

    # temperature
    T = potential_temperature * exner
    # T = potential_temperature - g / (c_p) * x[2]

    # density
    rho = p / (R * T)

    # Geopotential
    phi = g * x[2]

    v1 = 20
    v2 = 0
    E = c_v * T + 0.5f0 * (v1^2 + v2^2) + phi
    return SVector(rho, rho * v1, rho * v2, rho * E, phi)
end

###############################################################################
# semidiscretization of the compressible Euler equations
warm_bubble_setup = WarmBubbleSetup()

equations = CompressibleEulerEquationsWithGravity2D(warm_bubble_setup.gamma)

initial_condition = warm_bubble_setup

volume_flux = (flux_kennedy_gruber, flux_nonconservative_waruszewski)
surface_flux = (FluxLMARS(340.0), flux_nonconservative_waruszewski)

polydeg = 3
solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

coordinates_min = (0.0, -5000.0)
coordinates_max = (20_000.0, 15_000.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 6,
                n_cells_max = 10_000,
                periodicity = (true, false))

boundary_conditions = (; x_neg = boundary_condition_periodic,
                       x_pos = boundary_condition_periodic,
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 1000.0)  # 1000 seconds final time
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     analysis_polydeg = polydeg)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 10.0, #interval = 1, #dt = 10.0,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks);
