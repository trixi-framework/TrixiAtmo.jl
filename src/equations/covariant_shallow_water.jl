@muladd begin
#! format: noindent

struct CovariantShallowWaterEquations2D{RealT <: Real} <: AbstractCovariantEquations{2,
                                  3, 3}
    gravitational_acceleration::RealT
    rotation_rate::RealT
    function CovariantShallowWaterEquations2D(gravitational_acceleration::RealT,
                                              rotation_rate::RealT) where {RealT <:
                                                                           Real}
        return new{RealT}(gravitational_acceleration, rotation_rate)
    end
end

Trixi.have_nonconservative_terms(::CovariantShallowWaterEquations2D) = True()

function Trixi.varnames(::typeof(cons2cons), ::CovariantShallowWaterEquations2D)
    return ("h", "h_vcon1", "h_vcon2")
end

# Compute the entropy variables (requires element container and indices)
@inline function Trixi.cons2entropy(u, equations::CovariantShallowWaterEquations2D,
                                    elements, i, j, element)
    h, h_vcon1, h_vcon2 = u
    Gcov = SMatrix{2, 2}(view(elements.covariant_metric, :, :, i, j, element))
    vcon = SVector(h_vcon1 / h, h_vcon2 / h)
    vcov = Gcov * vcon
    w1 = equations.gravitational_acceleration * h - 0.5f0 * dot(vcov, vcon)
    return SVector{3}(w1, vcov[1], vcov[2])
end

# Convert contravariant velocity/momentum components to zonal and meridional components
@inline function contravariant2spherical(u::SVector{3},
                                         ::CovariantShallowWaterEquations2D,
                                         elements, i, j, element)
    h_vlon, h_vlat = contravariant2spherical(u[2], u[3], elements, i, j, element)
    return SVector(u[1], h_vlon, h_vlat)
end

# Convert zonal and meridional velocity/momentum components to contravariant components
@inline function spherical2contravariant(u::SVector{3},
                                         ::CovariantShallowWaterEquations2D,
                                         elements, i, j, element)
    h_vcon1, h_vcon2 = spherical2contravariant(u[2], u[3], elements, i, j, element)
    return SVector(u[1], h_vcon1, h_vcon2)
end

# The flux for the covariant form takes in the element container and node/element indices
# in order to give the flux access to the geometric information
@inline function Trixi.flux(u, orientation::Integer,
                            equations::CovariantShallowWaterEquations2D,
                            elements, i, j, element)
    h, h_vcon1, h_vcon2 = u
    (; contravariant_metric) = elements

    J = volume_element(elements, i, j, element)
    gravitational_term = 0.5f0 * equations.gravitational_acceleration * h^2
    h_vcon = SVector(h_vcon1, h_vcon2)
    vcon_orientation = h_vcon[orientation] / h

    momentum_flux_1 = h_vcon1 * vcon_orientation +
                      contravariant_metric[1, orientation, i, j, element] *
                      gravitational_term
    momentum_flux_2 = h_vcon2 * vcon_orientation +
                      contravariant_metric[2, orientation, i, j, element] *
                      gravitational_term

    return SVector(J * h_vcon[orientation], J * momentum_flux_1, J * momentum_flux_2)
end

# Directional flux that takes in the normal components in reference space
@inline function Trixi.flux(u, normal_direction::AbstractVector,
                            equations::CovariantShallowWaterEquations2D,
                            elements, i, j, element)
    fcon1 = Trixi.flux(u, 1, equations, elements, i, j, element)
    fcon2 = Trixi.flux(u, 2, equations, elements, i, j, element)
    return fcon1 * normal_direction[1] + fcon2 * normal_direction[2]
end

# Symmetric part of entropy-conservative flux
@inline function flux_split_covariant(u_ll, u_rr, orientation::Integer,
                                      equations::CovariantShallowWaterEquations2D,
                                      elements, i_ll, j_ll, i_rr, j_rr, element)
    h_ll, h_vcon1_ll, h_vcon2_ll = u_ll
    h_rr, h_vcon1_rr, h_vcon2_rr = u_rr

    J_ll = volume_element(elements, i_ll, j_ll, element)
    J_rr = volume_element(elements, i_rr, j_rr, element)

    h_vcon_ll = SVector(h_vcon1_ll, h_vcon2_ll)
    h_vcon_rr = SVector(h_vcon1_rr, h_vcon2_rr)

    # Scaled mass flux in conservative form
    mass_flux = 0.5f0 * (J_ll * h_vcon_ll[orientation] + J_rr * h_vcon_rr[orientation])

    # Half of scaled inertial flux in conservative form
    momentum_flux = 0.25f0 * (J_ll * h_vcon_ll * h_vcon_ll[orientation] / h_ll +
                     J_rr * h_vcon_rr * h_vcon_rr[orientation] / h_rr)

    return SVector(mass_flux, momentum_flux[1], momentum_flux[2])
end

# Split-covariant flux with local Lax-Friedrichs dissipation
const flux_split_covariant_lax_friedrichs = FluxPlusDissipation(flux_split_covariant,
                                                                DissipationLocalLaxFriedrichs(max_abs_speed_naive))

# Non-symmetric part of entropy-conservative flux
@inline function flux_nonconservative_split_covariant(u_ll, u_rr,
                                                      orientation::Integer,
                                                      equations::CovariantShallowWaterEquations2D,
                                                      elements, i_ll, j_ll, i_rr, j_rr,
                                                      element)
    h_ll, h_vcon1_ll, h_vcon2_ll = u_ll
    h_rr, h_vcon1_rr, h_vcon2_rr = u_rr

    # Geometric variables
    Gcov_ll = SMatrix{2, 2}(view(elements.covariant_metric, :, :, i_ll, j_ll, element))
    Gcov_rr = SMatrix{2, 2}(view(elements.covariant_metric, :, :, i_rr, j_rr, element))
    Gcon_ll = SMatrix{2, 2}(view(elements.contravariant_metric, :, :, i_ll, j_ll,
                                 element))
    J_ll = volume_element(elements, i_ll, j_ll, element)
    J_rr = volume_element(elements, i_rr, j_rr, element)

    # Contravariant and covariant momentum and velocity components
    h_vcon_ll = SVector(h_vcon1_ll, h_vcon2_ll)
    h_vcon_rr = SVector(h_vcon1_rr, h_vcon2_rr)
    vcov_ll = Gcov_ll * h_vcon_ll / h_ll
    vcov_rr = Gcov_rr * h_vcon_rr / h_rr

    # Half of inertial term in non-conservative form
    nonconservative_term_inertial = 0.5f0 * Gcon_ll *
                                    (J_ll * h_vcon_ll[orientation] * vcov_rr +
                                     J_rr * h_vcon_rr[orientation] * vcov_ll)

    # Gravity term in non-conservative form
    nonconservative_term_gravitational = equations.gravitational_acceleration *
                                         J_ll * Gcon_ll[:, orientation] *
                                         h_ll * h_rr # the same as h_ll * (h_ll + h_rr)

    nonconservative_term = nonconservative_term_inertial +
                           nonconservative_term_gravitational

    return SVector(zero(eltype(u_ll)), nonconservative_term[1], nonconservative_term[2])
end

@inline function source_terms_weak_form(u, x, t,
                                        equations::CovariantShallowWaterEquations2D,
                                        elements, i, j, element)
    h, h_vcon1, h_vcon2 = u

    # Geometric variables
    Gcon = SMatrix{2, 2}(view(elements.contravariant_metric, :, :, i, j, element))
    Gamma1 = SMatrix{2, 2}(view(elements.christoffel_symbols, :, :, 1, i, j, element))
    Gamma2 = SMatrix{2, 2}(view(elements.christoffel_symbols, :, :, 2, i, j, element))
    J = volume_element(elements, i, j, element)

    # Physical variables
    h_vcon = SVector{2}(h_vcon1, h_vcon2)
    v_con = SVector{2}(h_vcon1 / h, h_vcon2 / h)

    # Doubly-contravariant and mixed inertial flux tensors
    T = h_vcon * v_con' + 0.5f0 * equations.gravitational_acceleration * h^2 * Gcon

    # Coriolis parameter
    f = 2 * equations.rotation_rate * x[3] / norm(x)  # 2Ωsinθ

    # Geometric source term
    s_geo = SVector(sum(Gamma1 .* T), sum(Gamma2 .* T))

    # Combined source terms
    source_1 = s_geo[1] + f * J * (Gcon[1, 2] * h_vcon1 - Gcon[1, 1] * h_vcon2)
    source_2 = s_geo[2] + f * J * (Gcon[2, 2] * h_vcon1 - Gcon[2, 1] * h_vcon2)

    # Do not scale by Jacobian since apply_jacobian! is called before this
    return SVector(zero(eltype(u)), -source_1, -source_2)
end

# Geometric and Coriolis source terms for a rotating sphere
@inline function source_terms_split_covariant(u, x, t,
                                              equations::CovariantShallowWaterEquations2D,
                                              elements, i, j, element)
    h, h_vcon1, h_vcon2 = u

    # Geometric variables
    Gcov = SMatrix{2, 2}(view(elements.covariant_metric, :, :, i, j, element))
    Gcon = SMatrix{2, 2}(view(elements.contravariant_metric, :, :, i, j, element))
    Gamma1 = SMatrix{2, 2}(view(elements.christoffel_symbols, :, :, 1, i, j, element))
    Gamma2 = SMatrix{2, 2}(view(elements.christoffel_symbols, :, :, 2, i, j, element))
    J = volume_element(elements, i, j, element)

    # Physical variables
    h_vcon = SVector{2}(h_vcon1, h_vcon2)
    v_con = SVector{2}(h_vcon1 / h, h_vcon2 / h)

    # Doubly-contravariant and mixed inertial flux tensors
    h_vcon_vcon = h_vcon * v_con'
    h_vcov_vcon = Gcov * h_vcon_vcon

    # Coriolis parameter
    f = 2 * equations.rotation_rate * x[3] / norm(x)  # 2Ωsinθ

    # Geometric source term
    s_geo = 0.5f0 * (SVector(sum(Gamma1 .* h_vcon_vcon), sum(Gamma2 .* h_vcon_vcon)) -
             Gcon * (Gamma1 * h_vcov_vcon[1, :] + Gamma2 * h_vcov_vcon[2, :]))

    # Combined source terms
    source_1 = s_geo[1] + f * J * (Gcon[1, 2] * h_vcon1 - Gcon[1, 1] * h_vcon2)
    source_2 = s_geo[2] + f * J * (Gcon[2, 2] * h_vcon1 - Gcon[2, 1] * h_vcon2)

    # Do not scale by Jacobian since apply_jacobian! is called before this
    return SVector(zero(eltype(u)), -source_1, -source_2)
end

# Maximum wave speed along the normal direction in reference space
@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, orientation,
                                           equations::CovariantShallowWaterEquations2D,
                                           elements, i_ll, j_ll, i_rr, j_rr, element)
    h_ll, h_vcon1_ll, h_vcon2_ll = u_ll
    h_rr, h_vcon1_rr, h_vcon2_rr = u_rr

    h_vcon_ll = SVector(h_vcon1_ll, h_vcon2_ll)
    h_vcon_rr = SVector(h_vcon1_rr, h_vcon2_rr)

    Gcon = elements.contravariant_metric[orientation, orientation, i_ll, j_ll, element]

    phi_ll = max(h_ll * equations.gravitational_acceleration, 0)
    phi_rr = max(h_rr * equations.gravitational_acceleration, 0)

    return max(abs(h_vcon_ll[orientation] / h_ll) + sqrt(Gcon * phi_ll),
               abs(h_vcon_rr[orientation] / h_rr) + sqrt(Gcon * phi_rr))
end

# Maximum wave speeds with respect to the covariant basis
@inline function Trixi.max_abs_speeds(u, equations::CovariantShallowWaterEquations2D,
                                      elements, i, j, element)
    h, h_vcon1, h_vcon2 = u
    Gcon_11 = elements.contravariant_metric[1, 1, i, j, element]
    Gcon_22 = elements.contravariant_metric[2, 2, i, j, element]
    phi = max(h * equations.gravitational_acceleration, 0)
    return abs(h_vcon1 / h) + sqrt(Gcon_11 * phi),
           abs(h_vcon2 / h) + sqrt(Gcon_22 * phi)
end

# Steady geostrophically balanzed zonal flow
function Trixi.initial_condition_convergence_test(x, t,
                                                  equations::CovariantShallowWaterEquations2D)
    (; gravitational_acceleration, rotation_rate) = equations

    radius = sqrt(x[1]^2 + x[2]^2 + x[3]^2)
    lat = asin(x[3] / radius)

    # compute zonal and meridional components of the velocity
    V = convert(eltype(x), 2π) * radius / (12 * SECONDS_PER_DAY)
    vlon, vlat = V * cos(lat), zero(eltype(x))

    # compute geopotential height 
    h = 1 / gravitational_acceleration *
        (2.94f4 - (radius * rotation_rate * V + 0.5f0 * V^2) * (sin(lat))^2)

    # convert to conservative variables
    return SVector(h, h * vlon, h * vlat)
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

    vlon = radius * K * cos(lat) +
           radius * K * (cos(lat))^(R - 1) * (R * (sin(lat))^2 - (cos(lat))^2) *
           cos(R * lon)
    vlat = -radius * K * R * (cos(lat))^(R - 1) * sin(lat) * sin(R * lon)

    # convert to conservative variables
    return SVector(h, h * vlon, h * vlat)
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
    lat_0, lat_1 = convert(realT, π / 7), convert(realT, π / 2) - lat_0
    vlon = galewsky_velocity(lat, u_0, lat_0, lat_1)
    vlat = zero(eltype(x))

    # numerically integrate (here we use the QuadGK package) to get height
    galewsky_integral, _ = quadgk(latp -> galewsky_integrand(latp, u_0, lat_0, lat_1,
                                                             radius,
                                                             equations), -π / 2, lat)
    h = 10158.0f0 - radius / gravitational_acceleration * galewsky_integral

    # add perturbation to initiate instability
    α, β = convert(realT, 1 / 3), convert(realT, 1 / 15)
    lat_2 = convert(realT, π / 4)
    if (-π < lon) && (lon < π)
        h = h + 120.0f0 * cos(lat) * exp(-((lon / α)^2)) * exp(-((lat_2 - lat) / β)^2)
    end
    # convert to conservative variables
    return SVector(h, h * vlon, h * vlat)
end
end # @muladd
