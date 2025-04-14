using Trixi
using OrdinaryDiffEq

struct SchärSetup
    # Physical constants
    g::Float64       # gravity of earth
    c_p::Float64     # heat capacity for constant pressure (dry air)
    c_v::Float64     # heat capacity for constant volume (dry air)
    gamma::Float64   # heat capacity ratio (dry air)
    p_0::Float64     # atmospheric pressure
    theta_0::Float64 # 
    u0::Float64      #
    Nf::Float64      # 
    z_B::Float64     # start damping layer
    z_T::Float64     # end damping layer
    function SchärSetup(; g = 9.81, c_p = 1004.0, c_v = 717.0, gamma = c_p / c_v,
                        p_0 = 100_000.0, theta_0 = 280.0, u0 = 10.0, z_B = 16000.0,
                        z_T = 21000.0)
        Nf = 0.01
        new(g, c_p, c_v, gamma, p_0, theta_0, u0, Nf, z_B, z_T)
    end
end

function (setup::SchärSetup)(u, x, t, equations::CompressibleEulerEquations2D)
    @unpack g, c_p, c_v, gamma, p_0, theta_0, z_B, z_T, Nf, u0 = setup

    rho, rho_v1, rho_v2, rho_e = u

    R = c_p - c_v

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2))
    rho_theta = (p / p_0)^(c_v / c_p) * p_0 / R
    theta = rho_theta / rho

    alfa = 0.1

    xr_B = 20000.0
    xr_T = 25000.0

    if x[2] <= z_B
        S_v = 0.0
    else
        S_v = -alfa * sinpi(0.5 * (x[2] - z_B) / (z_T - z_B))^2
    end

    if x[1] < xr_B
        S_h1 = 0.0
    else
        S_h1 = -alfa * sinpi(0.5 * (x[1] - xr_B) / (xr_T - xr_B))^2
    end

    if x[1] > -xr_B
        S_h2 = 0.0
    else
        S_h2 = -alfa * sinpi(0.5 * (x[1] + xr_B) / (-xr_T + xr_B))^2
    end

    theta_b = theta_0 * exp(Nf^2 / g * x[2])

    K = p_0 * (R / p_0)^gamma
    du2 = rho * (v1 - u0) * (S_v + S_h1 + S_h2)
    du3 = rho_v2 * (S_v + S_h1 + S_h2)
    du4 = rho * (theta - theta_b) * (S_v + S_h1 + S_h2) * K * gamma / (gamma - 1.0) *
          (rho_theta)^(gamma - 1.0) + du2 * v1 + du3 * v2

    return SVector(zero(eltype(u)), du2, du3 - g * rho, du4 - g * rho_v2)
end

function (setup::SchärSetup)(x, t,
                             equations::Union{CompressibleEulerEquations2D,
                                              CompressibleEulerPotentialTemperatureEquations2D})
    @unpack g, c_p, c_v, p_0, theta_0, u0, Nf = setup

    # Exner pressure, solves hydrostatic equation for x[2]
    exner = 1 + g^2 / (c_p * theta_0 * Nf^2) * (exp(-Nf^2 / g * x[2]) - 1)
    # pressure
    p_0 = 100_000.0  # reference pressure
    R = c_p - c_v    # gas constant (dry air)
    p = p_0 * exner^(c_p / R)
    potential_temperature = theta_0 * exp(Nf^2 / g * x[2])
    T = potential_temperature * exner
    # density
    rho = p / (R * T)
    v1 = u0
    v2 = 0.0

    return prim2cons(SVector(rho, v1, v2, p), equations)
end

schär_setup = SchärSetup()

boundary_conditions = (x_neg = boundary_condition_periodic,
                       x_pos = boundary_condition_periodic,
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

polydeg = 3
basis = LobattoLegendreBasis(polydeg)
surface_flux = FluxLMARS(340.0)
volume_flux = flux_kennedy_gruber
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

a = 5000.0
L = 50000.0
H = 21000.0
lambda_c = 4000.0
hc = 250.0
y_b = hc * exp(-(L / 2 / a)^2) * cospi(L / 2 / lambda_c)^2
alfa = (H - y_b) * 0.5

f1(s) = SVector(-L / 2, y_b + alfa * (s + 1))
f2(s) = SVector(L / 2, y_b + alfa * (s + 1))
f3(s) = SVector(s * L / 2, hc * exp(-(s * L / 2 / a)^2) * cospi(s * L / 2 / lambda_c)^2)
f4(s) = SVector(s * L / 2, H)

mesh = StructuredMesh(cells_per_dimension, (f1, f2, f3, f4), periodicity = (true, false))

semi = SemidiscretizationHyperbolic(mesh, equations, schär_setup, solver,
                                    source_terms = schär_setup,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.
T = 5
tspan = (0.0, T * 3600.0) #1 hour.

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000

analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:entropy_conservation_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback)

###############################################################################
# run the simulation
sol = solve(ode,
            SSPRK43(thread = OrdinaryDiffEq.True()),
            maxiters = 1.0e7,
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks)

summary_callback()