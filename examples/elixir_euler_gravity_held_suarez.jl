# Held-Suarez test case
# Following Souza et al.:
# The Flux-Differencing Discontinuous Galerkin Method Applied to an Idealized Fully
# Compressible Nonhydrostatic Dry Atmosphere

using OrdinaryDiffEq
using Trixi, TrixiAtmo
using LinearAlgebra


###############################################################################
# semidiscretization of the compressible Euler equations
gamma = 1.4
equations = CompressibleEulerEquationsWithGravity3D(gamma)

function initial_condition_isothermal(x, t, equations::CompressibleEulerEquationsWithGravity3D)
    # equation (60) in the paper
    temperature = 285                     # T_I
    gas_constant = 287                    # R_d
    surface_pressure = 1e5                # p_0
    radius_earth = 6.371229e6             # r_planet
    gravitational_acceleration = 9.80616  # g
    c_v = 717.5                 # Specific heat capacity of dry air at constant volume

    r = norm(x)
    # Make sure that r is not smaller than radius_earth
    z = max(r - radius_earth, 0.0)
    r = z + radius_earth

    # pressure
    # geopotential formulation?
    p = surface_pressure *
        exp(gravitational_acceleration *
            (radius_earth^2 / r - radius_earth) /
            (gas_constant * temperature))

    # density (via ideal gas law)
    rho = p / (gas_constant * temperature)

    # geopotential
    phi = gravitational_acceleration * (radius_earth - radius_earth^2 / r)

    E = c_v * temperature + phi

    return SVector(rho, 0, 0, 0, rho * E, phi)
end

@inline function source_terms_coriolis(u, x, t,
                                       equations::CompressibleEulerEquationsWithGravity3D)
    radius_earth = 6.371229e6             # r_planet
    angular_velocity = 7.29212e-5         # Ω
    #                  7.27220521664e-05  (Giraldo)

    r = norm(x)
    # Make sure that r is not smaller than radius_earth
    z = max(r - radius_earth, 0.0)
    r = z + radius_earth

    # Coriolis term, -2Ω × ρv = -2 * angular_velocity * (0, 0, 1) × u[2:4]
    du2 =  2 * angular_velocity * u[3]
    du3 = -2 * angular_velocity * u[2]

    return SVector(0, du2, du3, 0, 0, 0)
end

function cartesian_to_sphere(x)
    r = norm(x)
    lambda = atan(x[2], x[1])
    if lambda < 0
        lambda += 2 * pi
    end
    phi = asin(x[3] / r)

    return lambda, phi, r
end

@inline function source_terms_hs_relaxation(u, x, t,
                                            equations::CompressibleEulerEquationsWithGravity3D)
    # equations (55)-(58) in the paper
    secondsperday = 60*60*24
    radius_earth = 6.371229e6   # r_planet
    k_f = 1/secondsperday       # Damping scale for momentum
    k_a = 1/(40*secondsperday)  # Polar relaxation scale
    k_s = 1/(4*secondsperday)   # Equatorial relaxation scale
    T_min = 200                 # Minimum equilibrium temperature
    T_equator = 315             # Equatorial equilibrium temperature
    surface_pressure = 1e5      # p_0
    deltaT = 60                 # Latitudinal temperature difference
    deltaTheta = 10             # Vertical temperature difference
    sigma_b = 0.7               # Dimensionless damping height
    gas_constant = 287          # R_d
    c_v = 717.5                 # Specific heat capacity of dry air at constant volume
    c_p = 1004.5                # Specific heat capacity of dry air at constant pressur
 
    _, _, _, _, pressure = cons2prim(u, equations)
    lon, lat, r = cartesian_to_sphere(x)
    temperature = pressure / (u[1] * gas_constant)

    sigma = pressure / surface_pressure   # "p_0 instead of instantaneous surface pressure"
    delta_sigma = max(0, (sigma-sigma_b)/(1-sigma_b))   # "height factor"
    k_v = k_f * delta_sigma
    k_T = k_a + (k_s - k_a) * delta_sigma * cos(lat)^4

    T_equi = max(T_min,
                 (T_equator - deltaT * sin(lat)^2 - deltaTheta * log(sigma) * cos(lat)^2) *
                  sigma^(gas_constant/c_p))

    # project onto r, normalize! @. Yₜ.c.uₕ -= k_v * Y.c.uₕ
    # Make sure that r is not smaller than radius_earth
    z = max(r - radius_earth, 0.0)
    r = z + radius_earth
    dotprod = (u[2] * x[1] + u[3] * x[2] + u[4] * x[3]) / (r*r)
    
    du2 = -k_v * (u[2] - dotprod * x[1])
    du3 = -k_v * (u[3] - dotprod * x[2])
    du4 = -k_v * (u[4] - dotprod * x[3])

    du5 = -k_T * u[1] * c_v * (temperature - T_equi)
 
    return SVector(0, du2, du3, du4, du5, 0)
end

@inline function source_terms_held_suarez(u, x, t,
                                          equations::CompressibleEulerEquationsWithGravity3D)
    return source_terms_coriolis(u,x,t,equations) +
           source_terms_hs_relaxation(u,x,t,equations)
end

initial_condition = initial_condition_isothermal

boundary_conditions = Dict(:inside => boundary_condition_slip_wall,
                           :outside => boundary_condition_slip_wall)

volume_flux = (flux_kennedy_gruber, flux_nonconservative_waruszewski)
surface_flux = (FluxLMARS(340.0), flux_nonconservative_waruszewski)

# Giraldo: (10,8), polydeg 4
# baroclinic instability: (16, 8), polydeg 5
lat_lon_trees_per_dim = 8
layers = 4
polydeg = 3

solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

mesh = Trixi.T8codeMeshCubedSphere(lat_lon_trees_per_dim, layers, 6.371229e6, 30000.0,
                                   polydeg = polydeg, initial_refinement_level = 0)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_held_suarez,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 365 * 24 * 60 * 60) # time in seconds for 1 year
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 5000,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim,
                                     output_directory = "out_heldsuarez")

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution)

###############################################################################
# run the simulation

# Use a Runge-Kutta method with automatic (error based) time step size control
# Enable threading of the RK method for better performance on multiple threads
sol = solve(ode,
            RDPK3SpFSAL49(thread = Trixi.True()); abstol = 1.0e-8, reltol = 1.0e-8,
            ode_default_options()..., maxiters=1e8,
            callback = callbacks);
