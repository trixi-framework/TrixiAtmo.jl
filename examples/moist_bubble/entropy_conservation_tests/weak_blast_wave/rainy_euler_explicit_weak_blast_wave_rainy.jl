using OrdinaryDiffEq
using Trixi
using TrixiAtmo
using TrixiAtmo: saturation_residual, saturation_residual_jacobian,
                 flux_ec_rain, saturation_vapour_pressure, flux_chandrashekar,
                 RainLimiterDG, cons2eq_pot_temp, entropy

function initial_condition_weak_blast_wave_rainy(x, t,
                                                 equations::CompressibleRainyEulerExplicitEquations2D)
    # constants
    c_l = equations.c_liquid_water
    c_vd = equations.c_dry_air_const_volume
    c_vv = equations.c_vapour_const_volume
    R_d = equations.R_dry_air
    R_v = equations.R_vapour
    L_ref = equations.ref_latent_heat_vap_temp

    if sqrt(x[1]^2 + x[2]^2) > 0.5
        rho = 1.0e-4
        v1 = 0.0
        v2 = 0.0
        p = 1.0
    else
        phi = atan(x[2], x[1])
        rho = 1.1691e-4
        v1 = 0.1882 * cos(phi)
        v2 = 0.1882 * sin(phi)
        p = 1.245
    end

    rho_dry = 0.25 * rho
    rho_vapour = 0.25 * rho
    rho_cloud = 0.25 * rho
    rho_rain = 0.25 * rho

    T = p / (rho_dry * R_d + rho_vapour * R_v)

    energy_density = (c_vd * rho_dry + c_vv * rho_vapour + c_l * (rho_cloud + rho_rain)) * T
    energy_density += L_ref * rho_vapour + 0.5 * rho * (v1^2 + v2^2)

    return SVector(rho_dry, rho_vapour, rho_cloud, rho_rain, rho * v1, rho * v2,
                   energy_density)
end

###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleRainyEulerExplicitEquations2D()

initial_condition = initial_condition_weak_blast_wave_rainy

polydeg = 3
basis = LobattoLegendreBasis(polydeg)

volume_flux = flux_ec_rain
surface_flux = volume_flux
solver = DGSEM(basis, surface_flux, VolumeIntegralFluxDifferencing(volume_flux))

coordinates_min = (-2.0, -2.0)
coordinates_max = (2.0, 2.0)

mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 5,
                n_cells_max = 1_000_000,
                periodicity = true)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition_periodic)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 0.4)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     save_analysis = true,
                                     extra_analysis_integrals = (entropy,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 1000,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = "out",
                                     solution_variables = cons2eq_pot_temp)

stepsize_callback = StepsizeCallback(cfl = 0.1)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks);
summary_callback() # print the timer summary
