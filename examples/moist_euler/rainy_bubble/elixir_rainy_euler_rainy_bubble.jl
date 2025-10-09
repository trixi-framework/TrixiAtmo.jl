using OrdinaryDiffEqSSPRK
using Trixi
using TrixiAtmo
using TrixiAtmo: source_terms_rainy, saturation_residual,
                 saturation_residual_jacobian, NonlinearSolveDG,
                 cons2eq_pot_temp, saturation_vapour_pressure,
                 flux_ec_rain, boundary_condition_simple_slip_wall,
                 generate_hydrostatic_residual, generate_perturbation_residual
using NLsolve: nlsolve

# domain
coordinates_min = (0.0, 0.0)
coordinates_max = (2400.0, 2400.0)

# create layers for initial condition
equations = CompressibleRainyEulerEquations2D()
layers = AtmosphereLayersRainyBubble(equations; total_height = coordinates_max[2] + 1)

function initial_condition_bubble_rainy(x, t,
                                        equations::CompressibleRainyEulerEquations2D{RealT};
                                        atmosphere_layers = layers) where {RealT}
    # equations constants
    c_vd = equations.c_dry_air_const_volume
    c_vv = equations.c_vapour_const_volume
    ref_L = equations.ref_latent_heat_vap_temp

    # problem specific constants
    humidity_rel_bar = convert(RealT, 0.2)  # background relative humidity field
    humidity_max = 1

    # bubble parameters
    radius_outer, radius_inner = 300, 200  # radii of humidity bubble
    x_center, z_center = 1200, 800  # center of humidity bubble

    # radius relative to bubble center
    r = sqrt((x[1] - x_center)^2 + (x[2] - z_center)^2)

    # humidity definition
    if (r > radius_outer)
        # outside the bubble
        humidity = humidity_rel_bar
    elseif (r > radius_inner)
        # outer layers of the bubble
        humidity = humidity_rel_bar +
                   (humidity_max - humidity_rel_bar) *
                   cos(pi * (r - radius_inner) / (2 * (radius_outer - radius_inner)))^2
    else
        # inner layer
        humidity = humidity_max
    end

    # get atmosphere layer and height information
    @unpack layer_data, total_height, precision = atmosphere_layers
    dz = precision
    z = x[2]
    n = convert(Int, floor((z + eps()) / dz)) + 1
    z_lower = (n - 1) * dz
    z_upper = n * dz

    if (z_lower == total_height)
        z_upper = z_lower + dz
        n = n - 1
    end

    # check height consistency
    if (z > total_height && !(isapprox(z, total_height)))
        error("The atmosphere does not match the simulation domain")
    end

    # get hydrostatic pressures and approximate between lower and upper data point
    pressure_hydrostatic_lower = layer_data[n, 1]
    pressure_hydrostatic_upper = layer_data[n + 1, 1]
    pressure_hydrostatic = (pressure_hydrostatic_upper * (z - z_lower) +
                            pressure_hydrostatic_lower * (z_upper - z)) / dz

    # solve perturbation
    residual_function! = generate_perturbation_residual(pressure_hydrostatic, humidity, z,
                                                        equations)
    rho_dry, rho_vapour, temperature = nlsolve(residual_function!, layer_data[n, 2:4],
                                               ftol = 1e-9, iterations = 20).zero

    energy_density = (c_vd * rho_dry + c_vv * rho_vapour) * temperature + rho_vapour * ref_L

    return SVector(rho_dry, rho_vapour, 0, 0, 0, energy_density, rho_vapour, 0,
                   temperature)
end

###############################################################################
# semidiscretization of the compressible rainy Euler equations

boundary_conditions = (x_neg = boundary_condition_periodic,
                       x_pos = boundary_condition_periodic,
                       y_neg = boundary_condition_simple_slip_wall,
                       y_pos = boundary_condition_simple_slip_wall)

polydeg = 3

surface_flux = flux_lax_friedrichs
volume_integral = VolumeIntegralFluxDifferencing(flux_ec_rain)

solver = DGSEM(polydeg, surface_flux, volume_integral)

cells_per_dimension = (64, 64)
mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                      periodicity = (true, false))

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_bubble_rainy, solver,
                                    source_terms = source_terms_rainy,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 600.0)

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000

# entropy
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:entropy_conservation_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 1000,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = "out",
                                     solution_variables = cons2eq_pot_temp)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution)

stage_limiter! = NonlinearSolveDG(saturation_residual, saturation_residual_jacobian,
                                  SVector(7, 8, 9))

###############################################################################
# run the simulation
sol = solve(ode, SSPRK43(stage_limiter!); ode_default_options()...,
            maxiters = 1.0e7, save_everystep = false, callback = callbacks);
