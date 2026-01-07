using OrdinaryDiffEqLowStorageRK
using Trixi, TrixiAtmo
using TrixiAtmo: cons2aeqpot, saturation_pressure, source_terms_moist_bubble
using NLsolve: nlsolve

###############################################################################
# semidiscretization of the compressible moist Euler equations

c_pd = 1004 # specific heat at constant pressure for dry air
c_vd = 717  # specific heat at constant volume for dry air
c_pv = 1885 # specific heat at constant pressure for moist air
c_vv = 1424 # specific heat at constant volume for moist air
equations = CompressibleMoistEulerEquations2D(c_pd = c_pd, c_vd = c_vd, c_pv = c_pv,
                                              c_vv = c_vv,
                                              gravity = EARTH_GRAVITATIONAL_ACCELERATION)

# Moist bubble test case from paper:
# G.H. Bryan, J.M. Fritsch, A Benchmark Simulation for Moist Nonhydrostatic Numerical
# Models, MonthlyWeather Review Vol.130, 2917–2928, 2002,
# https://journals.ametsoc.org/view/journals/mwre/130/12/1520-0493_2002_130_2917_absfmn_2.0.co_2.xml.
function initial_condition_moist_bubble(x, t, equations::CompressibleMoistEulerEquations2D,
                                        atmosphere_layers::AtmosphereLayers)
    @unpack layer_data, preciseness, total_height = atmosphere_layers
    dz = preciseness
    z = x[2]
    if (z > total_height && !(isapprox(z, total_height)))
        error("The atmosphere does not match the simulation domain")
    end
    n = convert(Int, floor((z + eps()) / dz)) + 1
    z_l = (n - 1) * dz
    (rho_l, rho_theta_l, rho_qv_l, rho_ql_l) = layer_data[n, :]
    z_r = n * dz
    if (z_l == total_height)
        z_r = z_l + dz
        n = n - 1
    end
    (rho_r, rho_theta_r, rho_qv_r, rho_ql_r) = layer_data[n + 1, :]
    rho = (rho_r * (z - z_l) + rho_l * (z_r - z)) / dz
    rho_theta = rho * (rho_theta_r / rho_r * (z - z_l) + rho_theta_l / rho_l * (z_r - z)) /
                dz
    rho_qv = rho * (rho_qv_r / rho_r * (z - z_l) + rho_qv_l / rho_l * (z_r - z)) / dz
    rho_ql = rho * (rho_ql_r / rho_r * (z - z_l) + rho_ql_l / rho_l * (z_r - z)) / dz

    rho, rho_e, rho_qv, rho_ql = perturb_moist_profile!(x, rho, rho_theta, rho_qv, rho_ql,
                                                        equations::CompressibleMoistEulerEquations2D)

    v1 = 0
    v2 = 0
    rho_v1 = rho * v1
    rho_v2 = rho * v2
    rho_E = rho_e + 1 / 2 * rho * (v1^2 + v2^2)

    return SVector(rho, rho_v1, rho_v2, rho_E, rho_qv, rho_ql)
end

function perturb_moist_profile!(x, rho, rho_theta, rho_qv, rho_ql,
                                equations::CompressibleMoistEulerEquations2D{RealT}) where {RealT}
    @unpack kappa, p_0, c_pd, c_vd, c_pv, c_vv, R_d, R_v, c_pl, L_00 = equations
    xc = 10000
    zc = 2000
    rc = 2000
    Δθ = 2

    r = sqrt((x[1] - xc)^2 + (x[2] - zc)^2)
    rho_d = rho - rho_qv - rho_ql
    kappa_M = (R_d * rho_d + R_v * rho_qv) / (c_pd * rho_d + c_pv * rho_qv + c_pl * rho_ql)
    p_loc = p_0 * (R_d * rho_theta / p_0)^(1 / (1 - kappa_M))
    T_loc = p_loc / (R_d * rho_d + R_v * rho_qv)
    rho_e = (c_vd * rho_d + c_vv * rho_qv + c_pl * rho_ql) * T_loc + L_00 * rho_qv

    # Assume pressure stays constant
    if (r < rc && Δθ > 0)
        # Calculate background density potential temperature
        θ_dens = rho_theta / rho * (p_loc / p_0)^(kappa_M - kappa)
        # Calculate perturbed density potential temperature
        θ_dens_new = θ_dens * (1 + Δθ * cospi(0.5f0 * r / rc)^2 / 300)
        rt = (rho_qv + rho_ql) / rho_d
        rv = rho_qv / rho_d
        # Calculate moist potential temperature
        θ_loc = θ_dens_new * (1 + rt) / (1 + (R_v / R_d) * rv)
        # Adjust varuables until the temperature is met
        if rt > 0
            while true
                T_loc = θ_loc * (p_loc / p_0)^kappa
                T_C = T_loc - convert(RealT, 273.15)
                # SaturVapor
                pvs = convert(RealT, 611.2) *
                      exp(convert(RealT, 17.62) * T_C / (convert(RealT, 243.12) + T_C))
                rho_d_new = (p_loc - pvs) / (R_d * T_loc)
                rvs = pvs / (R_v * rho_d_new * T_loc)
                θ_new = θ_dens_new * (1 + rt) / (1 + (R_v / R_d) * rvs)
                if abs(θ_new - θ_loc) <= θ_loc * 1.0e-12
                    break
                else
                    θ_loc = θ_new
                end
            end
        else
            rvs = 0
            T_loc = θ_loc * (p_loc / p_0)^kappa
            rho_d_new = p_loc / (R_d * T_loc)
            θ_new = θ_dens_new * (1 + rt) / (1 + (R_v / R_d) * rvs)
        end
        rho_qv = rvs * rho_d_new
        rho_ql = (rt - rvs) * rho_d_new
        rho = rho_d_new * (1 + rt)
        rho_d = rho - rho_qv - rho_ql
        kappa_M = (R_d * rho_d + R_v * rho_qv) /
                  (c_pd * rho_d + c_pv * rho_qv + c_pl * rho_ql)
        rho_theta = rho * θ_dens_new * (p_loc / p_0)^(kappa - kappa_M)
        rho_e = (c_vd * rho_d + c_vv * rho_qv + c_pl * rho_ql) * T_loc + L_00 * rho_qv
    end
    return SVector(rho, rho_e, rho_qv, rho_ql)
end

# Create background atmosphere data set
atmosphere_data = AtmosphereLayers(equations)

# Create the initial condition with the initial data set
function initial_condition_moist(x, t, equations)
    return initial_condition_moist_bubble(x, t, equations, atmosphere_data)
end

initial_condition = initial_condition_moist

boundary_condition = (x_neg = boundary_condition_slip_wall,
                      x_pos = boundary_condition_slip_wall,
                      y_neg = boundary_condition_slip_wall,
                      y_pos = boundary_condition_slip_wall)

source_term = source_terms_moist_bubble

###############################################################################
# Get the DG approximation space

polydeg = 4

surface_flux = FluxLMARS(360.0)
volume_flux = flux_chandrashekar

volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(polydeg, surface_flux, volume_integral)

coordinates_min = (0.0, 0.0)
coordinates_max = (20000.0, 10000.0)

cells_per_dimension = (64, 32)

# Create curved mesh with 64 x 32 elements
mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                      periodicity = (false, false))

###############################################################################
# create the semi discretization object

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition,
                                    source_terms = source_term)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 1000.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
solution_variables = cons2aeqpot

analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:entropy_conservation_error,),
                                     extra_analysis_integrals = (entropy, energy_total,
                                                                 saturation_pressure))

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

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks)
