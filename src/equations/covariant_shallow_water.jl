@muladd begin
#! format: noindent

struct CovariantShallowWaterEquations2D{RealT <: Real} <: AbstractCovariantEquations{2,
                                  3, 3}
    gravitational_acceleration::RealT
    rotation_rate::RealT
    splitting_coefficient::RealT # α in the splitting α∇ⱼτⁱʲ + (1-α)Gⁱᵏ∇ⱼτₖʲ
    function CovariantShallowWaterEquations2D(gravitational_acceleration::RealT,
                                              rotation_rate::RealT,
                                              alpha = convert(RealT, 1.0f0)) where {
                                                                                    RealT <:
                                                                                    Real
                                                                                    }
        return new{RealT}(gravitational_acceleration, rotation_rate, alpha)
    end
end

Trixi.have_nonconservative_terms(::CovariantShallowWaterEquations2D) = True()

function Trixi.varnames(::typeof(cons2cons), ::CovariantShallowWaterEquations2D)
    return ("h", "hv_con_1", "hv_con_2")
end

# Compute the entropy variables (requires element container and indices)
@inline function Trixi.cons2entropy(u, equations::CovariantShallowWaterEquations2D,
                                    elements, i, j, element)
    h, hv_con_1, hv_con_2 = u
    v_con_1, v_con_2 = hv_con_1 / h, hv_con_2 / h
    v_cov_1 = elements.covariant_metric[1, 1, i, j, element] * v_con_1 +
              elements.covariant_metric[1, 2, i, j, element] * v_con_2
    v_cov_2 = elements.covariant_metric[2, 1, i, j, element] * v_con_1 +
              elements.covariant_metric[2, 2, i, j, element] * v_con_2

    w1 = equations.gravitational_acceleration * h -
         0.5f0 * (v_cov_1 * v_con_1 + v_cov_2 * v_con_2)

    return SVector{3}(w1, v_cov_1, v_cov_2)
end

# Convert contravariant velocity/momentum components to zonal and meridional components
@inline function contravariant2spherical(u::SVector{3},
                                         ::CovariantShallowWaterEquations2D,
                                         elements, i, j, element)
    hv_lon, hv_lat = contravariant2spherical(u[2], u[3], elements, i, j, element)
    return SVector(u[1], hv_lon, hv_lat)
end

# Convert zonal and meridional velocity/momentum components to contravariant components
@inline function spherical2contravariant(u::SVector{3},
                                         ::CovariantShallowWaterEquations2D,
                                         elements, i, j, element)
    hv_con_1, hv_con_2 = spherical2contravariant(u[2], u[3], elements, i, j, element)
    return SVector(u[1], hv_con_1, hv_con_2)
end

# The flux for the covariant form takes in the element container and node/element indices
# in order to give the flux access to the geometric information
@inline function Trixi.flux(u, orientation::Integer,
                            equations::CovariantShallowWaterEquations2D,
                            elements, i, j, element)
    h, hv_con_1, hv_con_2 = u

    gravitational_term = 0.5f0 * equations.gravitational_acceleration * h^2

    G_con_1 = elements.contravariant_metric[1, orientation, i, j, element]
    G_con_2 = elements.contravariant_metric[2, orientation, i, j, element]

    hv = u[1 + orientation]
    v = hv / h
    T_con_1 = hv_con_1 * v + G_con_1 * gravitational_term
    T_con_2 = hv_con_2 * v + G_con_2 * gravitational_term

    return volume_element(elements, i, j, element) * SVector(hv, T_con_1, T_con_2)
end

# Directional flux that takes in the normal components in reference space
@inline function Trixi.flux(u, normal_direction::AbstractVector,
                            equations::CovariantShallowWaterEquations2D,
                            elements, i, j, element)
    F_con_1 = Trixi.flux(u, 1, equations, elements, i, j, element)
    F_con_2 = Trixi.flux(u, 2, equations, elements, i, j, element)
    return F_con_1 * normal_direction[1] + F_con_2 * normal_direction[2]
end

# Return a vector containing components (T_1^orientation and T_2^orientation)
# In other words, orientation is the contravariant index
@inline function momentum_flux_mixed(u, orientation::Integer,
                                     equations::CovariantShallowWaterEquations2D,
                                     elements, i, j, element)
    h, hv_con_1, hv_con_2 = u

    # Lower indices of momentum flux to obtain covariant components
    G_cov = SMatrix{2, 2}(view(elements.covariant_metric, :, :, i, j, element))
    hv_cov = G_cov * SVector(hv_con_1, hv_con_2)

    # Get contravariant velocity in correct orientation
    v_con = u[1 + orientation] / h

    # Add gravitational term in the correct orientation
    gravitational_term = 0.5f0 * equations.gravitational_acceleration * h^2
    return SVector(hv_cov[1] * v_con + I[orientation, 1] * gravitational_term,
                   hv_cov[2] * v_con + I[orientation, 2] * gravitational_term)
end

# Simple average but multiplied by splitting coefficient
@inline function flux_split_covariant(u_ll, u_rr,
                                      orientation_or_normal_direction,
                                      equations::CovariantShallowWaterEquations2D,
                                      elements, i_ll, j_ll, i_rr, j_rr, element)
    return equations.splitting_coefficient * Trixi.flux_central(u_ll, u_rr,
                              orientation_or_normal_direction, equations, elements,
                              i_ll, j_ll, i_rr, j_rr, element)
end

# Index raised outside the differential operator
@inline function flux_nonconservative_split_covariant(u_ll, u_rr,
                                                      orientation::Integer,
                                                      equations::CovariantShallowWaterEquations2D,
                                                      elements, i_ll, j_ll, i_rr, j_rr,
                                                      element)

    # Evaluate Jacobian at left and right
    J_ll = volume_element(elements, i_ll, j_ll, element)
    J_rr = volume_element(elements, i_rr, j_rr, element)

    # Evaluate symmetric contravariant mass flux
    Jhv_con_avg = J_ll * u_ll[1 + orientation] +
                  J_rr * u_rr[1 + orientation]

    # Evaluate symmetric mixed momentum flux (orientation is the contravariant index)
    JT_mixed_avg = J_ll * momentum_flux_mixed(u_ll, orientation, equations,
                                       elements, i_ll, j_ll, element) +
                   J_rr * momentum_flux_mixed(u_rr, orientation, equations,
                                       elements, i_rr, j_rr, element)

    # Contravariant metric tensor defined outside the differential operator (local part)
    G_ll_con = SMatrix{2, 2}(view(elements.contravariant_metric, :, :, i_ll, j_ll,
                                  element))

    # Non-conservative term is a product of local and symmetric parts
    return (1 - equations.splitting_coefficient) * SVector(Jhv_con_avg,
                   G_ll_con[1, 1] * JT_mixed_avg[1] + G_ll_con[1, 2] * JT_mixed_avg[2],
                   G_ll_con[2, 1] * JT_mixed_avg[1] + G_ll_con[2, 2] * JT_mixed_avg[2])
end

@inline function (dissipation::DissipationLocalLaxFriedrichs)(u_ll, u_rr,
                                                              normal_direction::AbstractVector,
                                                              equations::CovariantShallowWaterEquations2D,
                                                              elements, i, j, element)
    λ = dissipation.max_abs_speed(u_ll, u_rr, normal_direction, equations,
                                  elements, i, j, element)
    return -0.5f0 * volume_element(elements, i, j, element) * λ * (u_rr - u_ll)
end

# Geometric and Coriolis source terms for a rotating sphere
@inline function source_terms_split_covariant(u, x, t,
                                              equations::CovariantShallowWaterEquations2D,
                                              elements, i, j, element)
    h, hv_con_1, hv_con_2 = u

    # Geometric variables
    G_cov = SMatrix{2, 2}(view(elements.covariant_metric, :, :, i, j, element))
    G_con = SMatrix{2, 2}(view(elements.contravariant_metric, :, :, i, j, element))
    Gamma1 = SMatrix{2, 2}(view(elements.christoffel_symbols, :, :, 1, i, j, element))
    Gamma2 = SMatrix{2, 2}(view(elements.christoffel_symbols, :, :, 2, i, j, element))
    J = volume_element(elements, i, j, element)

    # Physical variables
    hv_con = SVector{2}(hv_con_1, hv_con_2)
    v_con = SVector{2}(hv_con_1 / h, hv_con_2 / h)

    T_con = hv_con * v_con' + 0.5f0 * equations.gravitational_acceleration * h^2 * G_con
    T_mix = G_cov * T_con

    # Coriolis parameter
    f = 2 * equations.rotation_rate * x[3] / norm(x)  # 2Ωsinθ

    # Geometric source term
    positive_geometric_source = (equations.splitting_coefficient) *
                                SVector(sum(Gamma1 .* T_con), sum(Gamma2 .* T_con))
    negative_geometric_source = (1 - equations.splitting_coefficient) *
                                G_con * (Gamma1 * T_mix[1, :] + Gamma2 * T_mix[2, :])
    s_geo = positive_geometric_source - negative_geometric_source

    # Combined source terms
    source_1 = s_geo[1] + f * J * (G_con[1, 2] * hv_con_1 - G_con[1, 1] * hv_con_2)
    source_2 = s_geo[2] + f * J * (G_con[2, 2] * hv_con_1 - G_con[2, 1] * hv_con_2)

    # Do not scale by Jacobian since apply_jacobian! is called before this
    return SVector(zero(eltype(u)), -source_1, -source_2)
end

# Maximum wave speed along the normal direction in reference space
@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, normal_direction,
                                           equations::CovariantShallowWaterEquations2D,
                                           elements, i, j, element)
    h_ll, hv_con_1_ll, hv_con_2_ll = u_ll
    h_rr, hv_con_1_rr, hv_con_2_rr = u_rr

    G_con = SMatrix{2, 2}(view(elements.contravariant_metric, :, :, i, j, element))
    G = normal_direction' * G_con * normal_direction

    v_ll = (hv_con_1_ll * normal_direction[1] + hv_con_2_ll * normal_direction[2]) /
           h_ll
    v_rr = (hv_con_1_rr * normal_direction[1] + hv_con_2_rr * normal_direction[2]) /
           h_rr

    phi_ll = max(h_ll * equations.gravitational_acceleration, 0)
    phi_rr = max(h_rr * equations.gravitational_acceleration, 0)

    return max(abs(v_ll) + sqrt(G * phi_ll), abs(v_rr) + sqrt(G * phi_rr))
end

# Maximum wave speeds with respect to the covariant basis
@inline function Trixi.max_abs_speeds(u, equations::CovariantShallowWaterEquations2D,
                                      elements, i, j, element)
    h, hv_con_1, hv_con_2 = u

    G_con_11 = elements.contravariant_metric[1, 1, i, j, element]
    G_con_22 = elements.contravariant_metric[2, 2, i, j, element]
    v_con_1 = hv_con_1 / h
    v_con_2 = hv_con_2 / h
    phi = max(h * equations.gravitational_acceleration, 0)
    return abs(v_con_1) + sqrt(G_con_11 * phi), abs(v_con_2) + sqrt(G_con_22 * phi)
end

# Steady geostrophically balanzed zonal flow
function Trixi.initial_condition_convergence_test(x, t,
                                                  equations::CovariantShallowWaterEquations2D)
    (; gravitational_acceleration, rotation_rate) = equations

    radius = sqrt(x[1]^2 + x[2]^2 + x[3]^2)
    lat = asin(x[3] / radius)

    # compute zonal and meridional components of the velocity
    V = convert(eltype(x), 2π) * radius / (12 * SECONDS_PER_DAY)
    v_lon, v_lat = V * cos(lat), zero(eltype(x))

    # compute geopotential height 
    h = 1 / gravitational_acceleration *
        (2.94f4 - (radius * rotation_rate * V + 0.5f0 * V^2) * (sin(lat))^2)

    # convert to conservative variables
    return SVector(h, h * v_lon, h * v_lat)
end

# Rossby-Haurwitz wave
function initial_condition_rossby_haurwitz(x, t,
                                           equations::CovariantShallowWaterEquations2D)
    (; gravitational_acceleration, rotation_rate) = equations

    radius = sqrt(x[1]^2 + x[2]^2 + x[3]^2)
    lon = atan(x[2], x[1])
    lat = asin(x[3] / radius)

    h_0 = 8.0f3
    K = 7.848f-6
    R = 4.0f0

    A = 0.5f0 * K * (2 * rotation_rate + K) * (cos(lat))^2 +
        0.25f0 * K^2 * (cos(lat))^(2 * R) *
        ((R + 1) * (cos(lat))^2 +
         (2 * R^2 - R - 2) - 2 * R^2 / ((cos(lat))^2))
    B = 2 * (rotation_rate + K) * K / ((R + 1) * (R + 2)) * (cos(lat))^R *
        ((R^2 + 2R + 2) - (R + 1)^2 * (cos(lat))^2)
    C = 0.25f0 * K^2 * (cos(lat))^(2 * R) * ((R + 1) * (cos(lat))^2 - (R + 2))

    h = h_0 +
        (1 / gravitational_acceleration) *
        (radius^2 * A + radius^2 * B * cos(R * lon) + radius^2 * C * cos(2 * R * lon))

    v_lon = radius * K * cos(lat) +
            radius * K * (cos(lat))^(R - 1) * (R * (sin(lat))^2 - (cos(lat))^2) *
            cos(R * lon)
    v_lat = -radius * K * R * (cos(lat))^(R - 1) * sin(lat) * sin(R * lon)

    # convert to conservative variables
    return SVector(h, h * v_lon, h * v_lat)
end

@inline function galewsky_velocity(θ, u_0, θ_0, θ_1)
    if (θ_0 < θ) && (θ < θ_1)
        u = u_0 / exp(-4 / (θ_1 - θ_0)^2) * exp(1 / (θ - θ_0) * 1 / (θ - θ_1))
    else
        u = zero(θ)
    end
    return u
end

@inline function galewsky_integrand(θ, u_0, θ_0, θ_1, a,
                                    equations::CovariantShallowWaterEquations2D)
    (; rotation_rate) = equations
    u = galewsky_velocity(θ, u_0, θ_0, θ_1)
    return u * (2 * rotation_rate * sin(θ) + u * tan(θ) / a)
end

function initial_condition_barotropic_instability(x, t,
                                                  equations::CovariantShallowWaterEquations2D)
    (; gravitational_acceleration, rotation_rate) = equations
    realT = eltype(x)
    radius = sqrt(x[1]^2 + x[2]^2 + x[3]^2)
    lon = atan(x[2], x[1])
    lat = asin(x[3] / radius)

    # compute zonal and meridional velocity components
    u_0 = 80.0f0
    lat_0 = convert(realT, π / 7)
    lat_1 = convert(realT, π / 2) - lat_0
    v_lon = galewsky_velocity(lat, u_0, lat_0, lat_1)
    v_lat = zero(eltype(x))

    # numerically integrate (here we use the QuadGK package) to get height
    galewsky_integral, _ = quadgk(latp -> galewsky_integrand(latp, u_0, lat_0, lat_1,
                                                             radius,
                                                             equations), -π / 2, lat)
    h = 10158.0f0 - radius / gravitational_acceleration * galewsky_integral

    # add perturbation to initiate instability
    α = convert(realT, 1 / 3)
    β = convert(realT, 1 / 15)
    lat_2 = convert(realT, π / 4)
    if (-π < lon) && (lon < π)
        h = h + 120.0f0 * cos(lat) * exp(-((lon / α)^2)) * exp(-((lat_2 - lat) / β)^2)
    end
    # convert to conservative variables
    return SVector(h, h * v_lon, h * v_lat)
end
end # @muladd
