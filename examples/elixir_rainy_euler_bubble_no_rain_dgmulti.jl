using OrdinaryDiffEq
using Trixi
using TrixiAtmo
using TrixiAtmo: source_terms_no_phase_change, saturation_residual,
                 saturation_residual_jacobian, NonlinearSolveDG,
                 cons2eq_pot_temp, flux_LMARS, flux_chandrashekar
using NLsolve: nlsolve
using Plots


# Initial condition from elixir_moist_euler_bubble.jl
function moist_state(y, dz, y0, r_t0, theta_e0, equations::CompressibleMoistEulerEquations2D)

    @unpack p_0, g, c_pd, c_pv, c_vd, c_vv, R_d, R_v, c_pl, L_00 = equations

    (p, rho, T, r_t, r_v, rho_qv, theta_e) = y
    p0 = y0[1]

    F     = zeros(7, 1)
    rho_d = rho / (1 + r_t)
    p_d   = R_d * rho_d * T
    T_C   = T - 273.15
    p_vs  = 611.2 * exp(17.62 * T_C / (243.12 + T_C))
    L     = L_00 - (c_pl - c_pv) * T

    F[1] = (p - p0) / dz + g * rho
    F[2] = p - (R_d * rho_d + R_v * rho_qv) * T
    # H = 1 is assumed
    F[3] = (theta_e -
            T * (p_d / p_0)^(-R_d / (c_pd + c_pl * r_t)) *
            exp(L * r_v / ((c_pd + c_pl * r_t) * T)))
    F[4] = r_t - r_t0
    F[5] = rho_qv - rho_d * r_v
    F[6] = theta_e - theta_e0
    a    = p_vs / (R_v * T) - rho_qv
    b    = rho - rho_qv - rho_d
    # H=1 => phi=0
    F[7] = a + b - sqrt(a * a + b * b)

    return F
end

function generate_function_of_y(dz, y0, r_t0, theta_e0, equations::CompressibleMoistEulerEquations2D)
    function function_of_y(y)
        return moist_state(y, dz, y0, r_t0, theta_e0, equations)
    end
end

struct AtmosphereLayers{RealT <: Real}
    equations::CompressibleMoistEulerEquations2D
    # structure:  1--> i-layer (z = total_height/precision *(i-1)),  2--> rho, rho_theta, rho_qv, rho_ql
    layer_data::Matrix{RealT}
    total_height::RealT
    preciseness::Int
    layers::Int
    ground_state::NTuple{2, RealT}
    equivalent_potential_temperature::RealT
    mixing_ratios::NTuple{2, RealT}
end

function AtmosphereLayers(equations; total_height = 10010.0, preciseness = 10,
                          ground_state = (1.4, 100000.0),
                          equivalent_potential_temperature = 320,
                          mixing_ratios = (0.02, 0.02), RealT = Float64)

    @unpack kappa, p_0, c_pd, c_vd, c_pv, c_vv, R_d, R_v, c_pl = equations
    rho0, p0 = ground_state
    r_t0, r_v0 = mixing_ratios
    theta_e0 = equivalent_potential_temperature

    rho_qv0 = rho0 * r_v0
    T0 = theta_e0
    y0 = [p0, rho0, T0, r_t0, r_v0, rho_qv0, theta_e0]

    n = convert(Int, total_height / preciseness)
    dz = 0.01
    layer_data = zeros(RealT, n + 1, 4)

    F = generate_function_of_y(dz, y0, r_t0, theta_e0, equations)
    sol = nlsolve(F, y0)
    p, rho, T, r_t, r_v, rho_qv, theta_e = sol.zero

    rho_d = rho / (1 + r_t)
    rho_ql = rho - rho_d - rho_qv
    kappa_M = (R_d * rho_d + R_v * rho_qv) / (c_pd * rho_d + c_pv * rho_qv + c_pl * rho_ql)
    rho_theta = rho * (p0 / p)^kappa_M * T * (1 + (R_v / R_d) * r_v) / (1 + r_t)

    layer_data[1, :] = [rho, rho_theta, rho_qv, rho_ql]
    for i in (1:n)
        y0 = deepcopy(sol.zero)
        dz = preciseness
        F = generate_function_of_y(dz, y0, r_t0, theta_e0, equations)
        sol = nlsolve(F, y0)
        p, rho, T, r_t, r_v, rho_qv, theta_e = sol.zero

        rho_d = rho / (1 + r_t)
        rho_ql = rho - rho_d - rho_qv
        kappa_M = (R_d * rho_d + R_v * rho_qv) /
                  (c_pd * rho_d + c_pv * rho_qv + c_pl * rho_ql)
        rho_theta = rho * (p0 / p)^kappa_M * T * (1 + (R_v / R_d) * r_v) / (1 + r_t)

        layer_data[i + 1, :] = [rho, rho_theta, rho_qv, rho_ql]
    end

    return AtmosphereLayers{RealT}(equations, layer_data, total_height, dz, n, ground_state,
                                   theta_e0, mixing_ratios)
end

# Moist bubble test case from paper:
# G.H. Bryan, J.M. Fritsch, A Benchmark Simulation for Moist Nonhydrostatic Numerical
# Models, MonthlyWeather Review Vol.130, 2917–2928, 2002,
# https://journals.ametsoc.org/view/journals/mwre/130/12/1520-0493_2002_130_2917_absfmn_2.0.co_2.xml.
function initial_condition_moist_bubble(x, t, equations::CompressibleMoistEulerEquations2D, atmosphere_layers::AtmosphereLayers)
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

    rho, rho_e, rho_qv, rho_ql, T_loc = perturb_moist_profile!(x, rho, rho_theta, rho_qv, rho_ql,
                                                        equations::CompressibleMoistEulerEquations2D)

    v1 = 0.0
    v2 = 0.0
    rho_v1 = rho * v1
    rho_v2 = rho * v2
    rho_E = rho_e + 1 / 2 * rho * (v1^2 + v2^2)

    return SVector(rho - rho_qv - rho_ql, rho_qv + rho_ql, 0.0, rho_v1, rho_v2, rho_E, rho_qv, rho_ql, T_loc)
end

function perturb_moist_profile!(x, rho, rho_theta, rho_qv, rho_ql,
                                equations::CompressibleMoistEulerEquations2D)
    @unpack kappa, p_0, c_pd, c_vd, c_pv, c_vv, R_d, R_v, c_pl, L_00 = equations
    xc = 10000.0
    zc = 2000.0
    rc = 2000.0
    Δθ = 2.0

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
        θ_dens_new = θ_dens * (1 + Δθ * cospi(0.5 * r / rc)^2 / 300)
        rt = (rho_qv + rho_ql) / rho_d
        rv = rho_qv / rho_d
        # Calculate moist potential temperature
        θ_loc = θ_dens_new * (1 + rt) / (1 + (R_v / R_d) * rv)
        # Adjust varuables until the temperature is met
        if rt > 0
            while true
                T_loc = θ_loc * (p_loc / p_0)^kappa
                T_C = T_loc - 273.15
                # SaturVapor
                pvs = 611.2 * exp(17.62 * T_C / (243.12 + T_C))
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
        # test
        #T_loc = θ_new * (p_loc / p_0)^kappa
        #
        rho_qv = rvs * rho_d_new
        rho_ql = (rt - rvs) * rho_d_new
        rho = rho_d_new * (1 + rt)
        rho_d = rho - rho_qv - rho_ql
        kappa_M = (R_d * rho_d + R_v * rho_qv) /
                  (c_pd * rho_d + c_pv * rho_qv + c_pl * rho_ql)
        rho_theta = rho * θ_dens_new * (p_loc / p_0)^(kappa - kappa_M)
        rho_e = (c_vd * rho_d + c_vv * rho_qv + c_pl * rho_ql) * T_loc + L_00 * rho_qv
    end
    return SVector(rho, rho_e, rho_qv, rho_ql, T_loc)
end

# Create background atmosphere data set
atmosphere_data = AtmosphereLayers(CompressibleMoistEulerEquations2D())

# Create the initial condition with the initial data set
function initial_condition_moist(x, t, equations::CompressibleRainyEulerEquations2D)
    return initial_condition_moist_bubble(x, t, CompressibleMoistEulerEquations2D(), atmosphere_data)
end

###############################################################################
# semidiscretization of the compressible moist Euler equations

equations = CompressibleRainyEulerEquations2D()

initial_condition = initial_condition_moist

# tag different boundary segments
left(x, tol = 50 * eps())   = abs(x[1] - coordinates_min[1]) < tol
right(x, tol = 50 * eps())  = abs(x[1] - coordinates_max[1]) < tol
bottom(x, tol = 50 * eps()) = abs(x[2] - coordinates_min[2]) < tol
top(x, tol = 50 * eps())    = abs(x[2] - coordinates_max[2]) < tol

is_on_boundary = Dict(:left => left, :right => right, :top => top, :bottom => bottom)

boundary_conditions = (; :left   => boundary_condition_periodic,
                         :top    => boundary_condition_slip_wall,
                         :bottom => boundary_condition_slip_wall,
                         :right  => boundary_condition_periodic)

#volume_flux  = flux_chandrashekar
#volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGMulti(polydeg = 1, element_type = Quad(), approximation_type = GaussSBP(),
                 surface_integral = SurfaceIntegralWeakForm(flux_lax_friedrichs),
                 volume_integral = VolumeIntegralWeakForm())

coordinates_min = (     0.0,      0.0)
coordinates_max = (20_000.0, 10_000.0)

cells_per_dimension = (200, 100)
mesh = DGMultiMesh(solver, cells_per_dimension; coordinates_min, coordinates_max, is_on_boundary, periodicity = (true, false))

semi = SemidiscretizationHyperbolic(mesh, equations,
                                    initial_condition, solver;
                                    source_terms = source_terms_no_phase_change,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 1000.0)

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000

analysis_callback = AnalysisCallback(semi, interval = analysis_interval, uEltype = real(solver))

alive_callback = AliveCallback(analysis_interval = 1000)

#save_solution = SaveSolutionCallback(interval = 1000,
                                     #save_initial_solution = true,
                                     #save_final_solution = true,
                                     #output_directory = "out",
                                     #solution_variables = cons2eq_pot_temp)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        #save_solution,
                        stepsize_callback)

stage_limiter! = NonlinearSolveDG(saturation_residual, saturation_residual_jacobian, SVector(7, 8, 9), 1e-9)

###############################################################################
# run the simulation
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false, stage_limiter!),
            maxiters = 1.0e7,
            dt = 1.0,
            save_everystep = false, callback = callbacks);

summary_callback()

pd = PlotData2D(sol; solution_variables = cons2eq_pot_temp);
plot(pd["eq_pot_temp"], c = :vik, dpi = 1000)