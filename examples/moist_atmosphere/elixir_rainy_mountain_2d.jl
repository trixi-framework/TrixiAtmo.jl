using TrixiAtmo
using TrixiAtmo: source_terms_no_phase_change, source_terms_rainy,
                 boundary_condition_simple_slip_wall,
                 saturation_vapour_pressure,
                 saturation_residual, saturation_residual_jacobian
using Trixi
using OrdinaryDiffEq
using Plots

###############################################################################
# Initial conditions
###############################################################################
# hydrostatic equilibrium
# near-surface temperature Ts = 283.15 K
# constant mean flow U = 10 m s−1
# a constant dry static stability Nd = 11 × 10−3 s−1
# a relative humidity profile, which is constant up to zm = 5 km
# and rapidly decreases above this level according to arctan
# with the near-surface humidity RHS = 95 % (RH95 case)
v1_0::Float64 = 10.0
v2_0::Float64 = 0.0

function initial_condition_rainy_mountain(x, t,
                                          equations::CompressibleRainyEulerEquations2D)
    g = equations.gravity
    c_pd = equations.c_dry_air_const_pressure
    c_vd = equations.c_dry_air_const_volume
    c_pv = equations.c_vapour_const_pressure
    c_vv = equations.c_vapour_const_volume
    c_pl = equations.c_liquid_water
    R_d = equations.R_dry_air
    R_v = equations.R_vapour
    ref_L = equations.ref_latent_heat_vap_temp

    p_0 = 100_000.0  # surface pressure
    theta_0 = 283.15     # surface potential temperature = surface temperature
    N = 1.1e-2     # Brunt-Väisälä frequency, constant dry static stability

    theta = theta_0 * exp(N^2 * x[2] / g)
    p_d = p_0 * (1 + g^2 / (c_pd * theta_0 * N^2) * (exp(-x[2] * N^2 / g) - 1))^(c_pd / R_d)
    T = (p_d / p_0)^(R_d / c_pd) * theta
    rho_d = p_d / (R_d * T)

    RH = 0.95
    if (x[2] > 5_000)
        RH = RH * (1 + 2 * inv(pi) * atan((5_000 - x[2]) / 500))
    end

    e_s = saturation_vapour_pressure(T, equations)
    p_v = RH * e_s
    rho_v = p_v / (R_v * T)

    rho_l = 0.0  # no cloud water
    rho_r = 0.0  # no rain

    rho_m = rho_v + rho_l
    rho = rho_d + rho_v + rho_l + rho_r

    rho_e = (c_vd * rho_d + c_vv * rho_v + c_pl * rho_l) * T + ref_L * rho_v
    rho_E = rho_e + 1 / 2 * rho * (v1_0^2 + v2_0^2)

    rho_v1 = rho * v1_0
    rho_v2 = rho * v2_0

    return SVector(rho_d, rho_m, rho_r, rho_v1, rho_v2, rho_E, rho_v, rho_l, T)
end

###############################################################################
# Source terms 
###############################################################################
# A Rayleigh damping layer above 11 km is applied
@inline function source_terms_rayleigh_sponge(u, x, t,
                                              equations::CompressibleRainyEulerEquations2D)
    rho_d, rho_m, rho_r, v1, v2, energy, rho_v, rho_c, T = cons2prim(u, equations)
    rho = rho_d + rho_m + rho_r

    # height of domain
    domain_height = 14_000
    # damping height
    damping_height = 11_000
    # relaxation coefficient > 0
    alpha = 0.5

    tau_s = zero(eltype(u))
    if x[2] > damping_height
        tau_s = alpha *
                sin(0.5 * (x[2] - damping_height) / (domain_height - damping_height))^2
    end

    return SVector(0.0, 0.0, 0.0,
                   -tau_s * rho * (v1 - v1_0),
                   -tau_s * rho * (v2 - v2_0),
                   0.0, 0.0, 0.0, 0.0)
end

@inline function source_terms_rainy_mountain(u, x, t,
                                             equations::CompressibleRainyEulerEquations2D)
    return (source_terms_rayleigh_sponge(u, x, t,
                                         equations::CompressibleRainyEulerEquations2D) +
                                         source_terms_rainy(u, x, t,
                                         equations::CompressibleRainyEulerEquations2D))
end

###############################################################################
# Solver
###############################################################################
polydeg = 3
basis = LobattoLegendreBasis(polydeg)

surface_flux = flux_lax_friedrichs
#surface_flux = flux_LMARS
#volume_flux  = flux_chandrashekar
#volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux)

###############################################################################
# Mesh
###############################################################################
function mapping(xi_, eta_)
    # transform input variables between -1 and 1

    # horizontal width = 553 km
    width = 553_000

    # domain height = 14 km
    height = 14_000

    # mountain peak height H = 1 km
    H = 1_000

    # half-width length a = 11 km
    a = 11_000

    x = xi_ * width / 2

    # bell-shaped mountain
    topo = H / ((x / a)^2 + 1)^1.5

    above = height - topo
    y = topo + (eta_ + 1) * above / 2

    return SVector(x, y)
end

# horizontal grid spacing of 2.765 km = 200 cells
# vertical spacing 200 m = 70 cells
cells_per_dimension = (200, 70)
mesh = StructuredMesh(cells_per_dimension, mapping,
                      periodicity = (true,false))
#trees_per_dimension = (100, 35)
#mesh = P4estMesh(trees_per_dimension,
#                 polydeg = polydeg,
#                 mapping = mapping_mountain,
#                 initial_refinement_level = 1,
#                 periodicity = false)

###############################################################################
# semidiscretization of the compressible rainy Euler equations
###############################################################################
equations = CompressibleRainyEulerEquations2D()

boundary_conditions_Dirichlet = BoundaryConditionDirichlet(initial_condition_rainy_mountain)
#boundary_conditions = Dict(:x_neg => boundary_conditions_Dirichlet,
#                           :x_pos => boundary_conditions_Dirichlet,
#                           :y_neg => boundary_condition_slip_wall,
#                           :y_pos => boundary_condition_slip_wall)
boundary_conditions = (x_neg = boundary_condition_periodic,
                       x_pos = boundary_condition_periodic,
                       y_neg = boundary_condition_simple_slip_wall,
                       y_pos = boundary_condition_simple_slip_wall)

semi = SemidiscretizationHyperbolic(mesh,
                                    equations,
                                    initial_condition_rainy_mountain,
                                    solver,
                                    source_terms = source_terms_rainy_mountain,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 10 * 60 * 60)

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000

analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:entropy_conservation_error,))

alive_callback = AliveCallback(analysis_interval = 1000)

save_solution = SaveSolutionCallback(interval = 2000,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = "out_rainy_mountain",
                                     solution_variables = cons2prim)

stepsize_callback = StepsizeCallback(cfl = 1.0)

#visualization = VisualizationCallback(interval = 100, show_mesh = false,
#                                      #suspend = true,
#                                      #plot_creator = Trixi.show_plot_makie,
#                                      solution_variables = cons2prim,
#                                      variable_names = ["rho_vapour", "rho_cloud", "rho_rain"],
#                                      aspect_ratio = 20)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
#                        visualization,
                        save_solution,
                        stepsize_callback)

stage_limiter! = NonlinearSolveDG(saturation_residual, saturation_residual_jacobian,
                                  SVector(7, 8, 9), 1e-9)

###############################################################################
# run the simulation
sol = solve(ode,
            CarpenterKennedy2N54(williamson_condition = false, stage_limiter!),
            maxiters = 1.0e7,
            dt = 1.0,
            save_everystep = false, callback = callbacks);

summary_callback()
