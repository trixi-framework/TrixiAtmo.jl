###############################################################################
# DGSEM for the linear advection equation on a prismed icosahedral grid
###############################################################################

using OrdinaryDiffEqLowStorageRK, Trixi, TrixiAtmo, LinearAlgebra


# Unperturbed balanced steady-state.
# Returns primitive variables with only the velocity in longitudinal direction (rho, u, p).
# The other velocity components are zero.
function basic_state_baroclinic_instability_longitudinal_velocity(lon, lat, z)
    # Parameters from Table 1 in the paper
    # Corresponding names in the paper are commented
    radius_earth = 6.371229e6  # a
    half_width_parameter = 2           # b
    gravitational_acceleration = 9.81     # g
    k = 3           # k
    surface_pressure = 1.0f5         # p₀
    gas_constant = 287         # R
    surface_equatorial_temperature = 310       # T₀ᴱ
    surface_polar_temperature = 240       # T₀ᴾ
    lapse_rate = 0.005       # Γ
    angular_velocity = 7.29212e-5  # Ω

    # Distance to the center of the Earth
    r = z + radius_earth

    # In the paper: T₀
    temperature0 = 0.5 * (surface_equatorial_temperature + surface_polar_temperature)
    # In the paper: A, B, C, H
    const_a = 1 / lapse_rate
    const_b = (temperature0 - surface_polar_temperature) /
              (temperature0 * surface_polar_temperature)
    const_c = 0.5 * (k + 2) * (surface_equatorial_temperature - surface_polar_temperature) /
              (surface_equatorial_temperature * surface_polar_temperature)
    const_h = gas_constant * temperature0 / gravitational_acceleration

    # In the paper: (r - a) / bH
    scaled_z = z / (half_width_parameter * const_h)

    # Temporary variables
    temp1 = exp(lapse_rate / temperature0 * z)
    temp2 = exp(-scaled_z^2)

    # In the paper: ̃τ₁, ̃τ₂
    tau1 = const_a * lapse_rate / temperature0 * temp1 +
           const_b * (1 - 2 * scaled_z^2) * temp2
    tau2 = const_c * (1 - 2 * scaled_z^2) * temp2

    # In the paper: ∫τ₁(r') dr', ∫τ₂(r') dr'
    inttau1 = const_a * (temp1 - 1) + const_b * z * temp2
    inttau2 = const_c * z * temp2

    # Temporary variables
    temp3 = r / radius_earth * cos(lat)
    temp4 = temp3^k - k / (k + 2) * temp3^(k + 2)

    # In the paper: T
    temperature = 1 / ((r / radius_earth)^2 * (tau1 - tau2 * temp4))

    # In the paper: U, u (zonal wind, first component of spherical velocity)
    big_u = gravitational_acceleration / radius_earth * k * temperature * inttau2 *
            (temp3^(k - 1) - temp3^(k + 1))
    temp5 = radius_earth * cos(lat)
    u = -angular_velocity * temp5 + sqrt(angular_velocity^2 * temp5^2 + temp5 * big_u)

    # Hydrostatic pressure
    p = surface_pressure *
        exp(-gravitational_acceleration / gas_constant * (inttau1 - inttau2 * temp4))

    # Density (via ideal gas law)
    rho = p / (gas_constant * temperature)
    return rho, u, p
end

# Perturbation as in Equations 25 and 26 of the paper (analytical derivative)
function perturbation_stream_function(lon, lat, z)
    # Parameters from Table 1 in the paper
    # Corresponding names in the paper are commented
    perturbation_radius = 1 / 6      # d₀ / a
    perturbed_wind_amplitude = 1      # Vₚ
    perturbation_lon = pi / 9     # Longitude of perturbation location
    perturbation_lat = 2 * pi / 9 # Latitude of perturbation location
    pertz = 15000    # Perturbation height cap

    # Great circle distance (d in the paper) divided by a (radius of the Earth)
    # because we never actually need d without dividing by a
    great_circle_distance_by_a = acos(sin(perturbation_lat) * sin(lat) +
                                      cos(perturbation_lat) * cos(lat) *
                                      cos(lon - perturbation_lon))

    # In the first case, the vertical taper function is per definition zero.
    # In the second case, the stream function is per definition zero.
    if z > pertz || great_circle_distance_by_a > perturbation_radius
        return 0, 0
    end

    # Vertical tapering of stream function
    perttaper = 1 - 3 * z^2 / pertz^2 + 2 * z^3 / pertz^3

    # sin/cos(pi * d / (2 * d_0)) in the paper
    sin_, cos_ = sincos(0.5f0 * pi * great_circle_distance_by_a / perturbation_radius)

    # Common factor for both u and v
    factor = 16 / (3 * sqrt(3)) * perturbed_wind_amplitude * perttaper * cos_^3 * sin_

    u_perturbation = -factor * (-sin(perturbation_lat) * cos(lat) +
                      cos(perturbation_lat) * sin(lat) * cos(lon - perturbation_lon)) /
                     sin(great_circle_distance_by_a)

    v_perturbation = factor * cos(perturbation_lat) * sin(lon - perturbation_lon) /
                     sin(great_circle_distance_by_a)
    return u_perturbation, v_perturbation
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

# Initial condition for an idealized baroclinic instability test
# https://doi.org/10.1002/qj.2241, Section 3.2 and Appendix A
function initial_condition_baroclinic_instability(x, t, equations)
    lon, lat, r = cartesian_to_sphere(x)
    radius_earth = EARTH_RADIUS
    # Make sure that the r is not smaller than radius_earth
    z = max(r - radius_earth, 0.0)

    # Unperturbed basic state
    rho, u, p = basic_state_baroclinic_instability_longitudinal_velocity(lon, lat, z)

    # Stream function type perturbation
    u_perturbation, v_perturbation = perturbation_stream_function(lon, lat, z)

    u += u_perturbation
    v = v_perturbation

    # Convert spherical velocity to Cartesian
    v1 = -sin(lon) * u - sin(lat) * cos(lon) * v
    v2 = cos(lon) * u - sin(lat) * sin(lon) * v
    v3 = cos(lat) * v

    return SVector(rho, v1, v2, v3, p)
end

# Geopotential function to pass as auxiliary_field keyword argument in constructor for 
# SemidiscretizationHyperbolic
@inline function geopotential_earth(x)
    r = norm(x)
    radius_earth = EARTH_RADIUS
    GM = radius_earth^2 * EARTH_GRAVITATIONAL_ACCELERATION
    # Geopotential relative to the surface value -> O(g*z), not O(GM/a)
    return GM * (1 / radius_earth - 1 / r)
end

###############################################################################
# Spatial discretization

initial_condition = initial_condition_baroclinic_instability

gamma = 1.4

equations = CovariantEulerEnergyEquations3D(c_p = 1004,
                                            c_v = 717,
                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)

###############################################################################
# Build DG solver.

polydeg = 2

dg = DGMulti(element_type = Hex(),
             approximation_type = SBP(),
             surface_flux = flux_central,
             polydeg = polydeg; Nplot = 3)

###############################################################################
# Build mesh.

horizontal_initial_refinement = 4
vertical_layers = 3

mesh = DGMultiMeshCubedShell3D(dg, EARTH_RADIUS, EARTH_RADIUS + 30000;
                                     horizontal_initial_refinement = horizontal_initial_refinement,
                                     vertical_layers = vertical_layers)

initial_condition_transformed = transform_initial_condition(initial_condition, equations)

# define boundary conditions
boundary_conditions = (; top = boundary_condition_slip_wall,
                       bottom = boundary_condition_slip_wall)

# A semidiscretization collects data structures and functions for the spatial discretization
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_transformed, dg,
                                    source_terms = source_terms_geometric_coriolis,
                                    boundary_conditions = boundary_conditions,
                                    auxiliary_field = geopotential_earth)


###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0 to T
tspan = (0.0, 12.0 * SECONDS_PER_DAY)
ode = semidiscretize(semi, tspan)

# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation 
# setup and resets the timers
summary_callback = SummaryCallback()

# The AnalysisCallback allows to analyse the solution in regular intervals and prints the 
# results
analysis_callback = AnalysisCallback(semi, interval = 1,
                                     save_analysis = true,
                                     uEltype = real(dg))

# The SaveSolutionCallback allows to save the solution to a file in regular intervals
save_solution = SaveSolutionCallback(interval = 1,
                                     solution_variables = contravariant_cons2global_prim)

# The StepsizeCallback handles the re-calculation of the maximum Δt after each time step
stepsize_callback = StepsizeCallback(cfl = 0.007)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE 
# solver
callbacks = CallbackSet(summary_callback, save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

# OrdinaryDiffEq's `solve` method evolves the solution in time and executes the passed 
# callbacks
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 1.0, save_everystep = false, callback = callbacks)