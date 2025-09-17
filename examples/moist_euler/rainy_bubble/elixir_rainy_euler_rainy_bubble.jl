using OrdinaryDiffEqSSPRK
using Trixi
using TrixiAtmo
using TrixiAtmo: source_terms_rainy, saturation_residual,
                 saturation_residual_jacobian, NonlinearSolveDG,
                 cons2eq_pot_temp, saturation_vapour_pressure,
                 flux_ec_rain, boundary_condition_simple_slip_wall
using NLsolve: nlsolve

# domain
coordinates_min = (0.0, 0.0)
coordinates_max = (2400.0, 2400.0)

# hydrostatic dry potential temperature
function theta_d(z, equations::CompressibleRainyEulerEquations2D)
    # constants
    c_pd = equations.c_dry_air_const_pressure
    R_d = equations.R_dry_air
    ref_pressure = equations.ref_pressure

    # problem specific constants
    surface_temperature = 283.0
    surface_pressure = 8.5e4
    stratification = 1.3e-5

    # dry potential temperature at surface
    Theta0 = surface_temperature * (ref_pressure / surface_pressure)^(R_d / c_pd)
    # at height z
    theta_d = Theta0 * exp(stratification * z)

    return theta_d
end

# hydrostatic base state residual
function generate_hydrostatic_residual(pressure_lower, humidity_rel0, z, dz,
                                       equations::CompressibleRainyEulerEquations2D)
    # equations constants
    c_pd = equations.c_dry_air_const_pressure
    R_d = equations.R_dry_air
    R_v = equations.R_vapour
    eps = equations.eps
    ref_pressure = equations.ref_pressure
    g = equations.gravity

    function hydrostatic_residual!(residual, guess)
        # variables
        pressure, rho_dry, rho_vapour, temperature = guess

        rho_vs = saturation_vapour_pressure(temperature, equations) / (R_v * temperature)

        # pressure derivative residual approximation
        residual[1] = (pressure - pressure_lower) / dz + (rho_dry + rho_vapour) * g

        # pressure residual
        residual[2] = pressure - temperature * (rho_dry * R_d + rho_vapour * R_v)

        # hydrostatic dry potential temperature residual
        residual[3] = theta_d(z, equations) -
                      temperature * (ref_pressure / pressure)^(R_d / c_pd)

        # humidity residual
        residual[4] = rho_vs * (rho_dry + rho_vapour / eps) * humidity_rel0
        residual[4] -= rho_vapour * (rho_dry + rho_vs / eps)
        residual[4] *= 1000.0
    end

    return hydrostatic_residual!
end

function generate_perturbation_residual(pressure_hydrostatic, H_init, z,
                                        equations::CompressibleRainyEulerEquations2D)
    # equations constants
    c_pd = equations.c_dry_air_const_pressure
    R_d = equations.R_dry_air
    R_v = equations.R_vapour
    eps = equations.eps
    ref_pressure = equations.ref_pressure

    function perturbation_residual!(residual, guess)
        # variables
        rho_dry, rho_vapour, temperature = guess

        rho_vs = saturation_vapour_pressure(temperature, equations) / (R_v * temperature)
        pressure = (rho_dry * R_d + rho_vapour * R_v) * temperature

        # humidity residual
        residual[1] = rho_vs * (rho_dry + rho_vapour / eps) * H_init
        residual[1] -= rho_vapour * (rho_dry + rho_vs / eps)
        residual[1] *= 30.0

        # hydrostatic dry potential temperature residual
        residual[2] = theta_d(z, equations) -
                      temperature * (ref_pressure / pressure_hydrostatic)^(R_d / c_pd)

        # pressure residual
        residual[3] = pressure_hydrostatic - pressure
    end

    return perturbation_residual!
end

# for approximating the dz pressure gradient
struct AtmosphereLayersRainyBubble{RealT <: Real}
    layer_data   :: Matrix{RealT}
    total_height :: RealT
    precision    :: RealT
end

function AtmosphereLayersRainyBubble(equations::CompressibleRainyEulerEquations2D;
                                     total_height = coordinates_max[2] + 1.0,
                                     precision = 1.0, RealT = Float64)
    # constants
    humidity_rel0 = 0.2      # hydrostatic relative humidity
    surface_pressure = 8.5e4

    # surface layer with initial guesses for rho_dry, rho_vapour and temperature
    surface_layer = [surface_pressure, 1.4, 0.04, 300.0]

    # allocate layer_data
    n = convert(Int, total_height / precision)
    layer_data = zeros(RealT, n + 1, 4)

    # solve (slightly above) surface layer
    dz = 0.01
    z = 0.01
    residual_function! = generate_hydrostatic_residual(surface_pressure, humidity_rel0, z,
                                                       dz, equations)
    layer_data[1, :] .= nlsolve(residual_function!, surface_layer).zero

    # adjust to chosen precision
    dz = precision

    # iterate up the atmosphere
    for i in (1:n)
        z += dz
        residual_function! = generate_hydrostatic_residual(layer_data[i, 1], humidity_rel0,
                                                           z, dz, equations)
        guess = deepcopy(layer_data[i, :])
        layer_data[i + 1, :] .= nlsolve(residual_function!, guess, ftol = 1e-10,
                                        iterations = 20).zero
    end

    return AtmosphereLayersRainyBubble{RealT}(layer_data, total_height, precision)
end

# create layers for initial condition
equations = CompressibleRainyEulerEquations2D()
layers = AtmosphereLayersRainyBubble(equations)

function initial_condition_bubble_rainy(x, t, equations::CompressibleRainyEulerEquations2D;
                                        atmosphere_layers = layers)
    # equations constants
    c_vd = equations.c_dry_air_const_volume
    c_vv = equations.c_vapour_const_volume
    ref_L = equations.ref_latent_heat_vap_temp

    # problem specific constants
    humidity_rel_bar = 0.2                # background relative humidity field
    humidity_max = 1.0

    # bubble parameters
    radius_outer, radius_inner = 300.0, 200.0      # radii of humidity bubble
    x_center, z_center = 1200.0, 800.0      # center of humidity bubble

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
                   cos(pi * (r - radius_inner) / (2.0 * (radius_outer - radius_inner)))^2
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

    return SVector(rho_dry, rho_vapour, 0.0, 0.0, 0.0, energy_density, rho_vapour, 0.0,
                   temperature)
end

###############################################################################
# semidiscretization of the compressible rainy Euler equations

boundary_conditions = (x_neg = boundary_condition_periodic,
                       x_pos = boundary_condition_periodic,
                       y_neg = boundary_condition_simple_slip_wall,
                       y_pos = boundary_condition_simple_slip_wall)

polydeg = 3
basis = LobattoLegendreBasis(polydeg)

surface_flux = flux_lax_friedrichs
volume_integral = VolumeIntegralFluxDifferencing(flux_ec_rain)

solver = DGSEM(basis, surface_flux, volume_integral)

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
