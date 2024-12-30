@muladd begin
#! format: noindent

struct CovariantShallowWaterEquations2D{GlobalCoordinateSystem, RealT <: Real} <:
       AbstractCovariantEquations{2, 3, GlobalCoordinateSystem, 3}
    gravity::RealT
    rotation_rate::RealT
    global_coordinate_system::GlobalCoordinateSystem
    function CovariantShallowWaterEquations2D(gravity::RealT,
                                              rotation_rate::RealT;
                                              global_coordinate_system = GlobalCartesianCoordinates()) where {RealT <:
                                                                                                              Real}
        return new{typeof(global_coordinate_system), RealT}(gravity, rotation_rate)
    end
end

# The conservative variables are the height and contravariant momentum components
function Trixi.varnames(::typeof(cons2cons), ::CovariantShallowWaterEquations2D)
    return ("h", "h_vcon1", "h_vcon2")
end

# Convenience function to extract the velocity
function velocity(u, ::CovariantShallowWaterEquations2D)
    return SVector(u[2] / u[1], u[3] / u[1])
end

# Convenience function to extract the momentum
function momentum(u, ::CovariantShallowWaterEquations2D)
    return SVector(u[2], u[3])
end

# Our implementation of flux-differencing formulation uses nonconservative terms, but the 
# standard weak form does not. To handle both options, we have defined a dummy kernel for 
# the nonconservative terms that does nothing when VolumeIntegralWeakForm is used with a 
# nonconservative system.
Trixi.have_nonconservative_terms(::CovariantShallowWaterEquations2D) = True()

# Entropy function (total energy per unit volume)
@inline function Trixi.entropy(u, aux_vars, equations::CovariantShallowWaterEquations2D)
    h, h_vcon1, h_vcon2 = u
    Gcov = metric_covariant(aux_vars, equations)
    vcon = SVector(h_vcon1 / h, h_vcon2 / h)
    vcov = Gcov * vcon
    return 0.5f0 * (dot(vcov, vcon) + equations.gravity * h^2)
end

@inline function Trixi.cons2prim(u, ::CovariantShallowWaterEquations2D)
    h, h_vcon1, h_vcon2 = u
    return SVector(h, h_vcon1 / h, h_vcon2 / h)
end

@inline function Trixi.prim2cons(u, ::CovariantShallowWaterEquations2D)
    h, vcon1, vcon2 = u
    return SVector(h, h * vcon1, h * vcon2)
end

# Entropy variables (partial derivatives of entropy with respect to conservative variables)
@inline function Trixi.cons2entropy(u, aux_vars,
                                    equations::CovariantShallowWaterEquations2D)
    h, h_vcon1, h_vcon2 = u
    Gcov = metric_covariant(aux_vars, equations)
    vcon = SVector(h_vcon1 / h, h_vcon2 / h)
    vcov = Gcov * vcon
    w1 = equations.gravity * h - 0.5f0 * dot(vcov, vcon)
    return SVector{3}(w1, vcov[1], vcov[2])
end


# Height and three global momentum components
function Trixi.varnames(::typeof(contravariant2global),
                        ::CovariantShallowWaterEquations2D)
    return ("h", "h_vglo1", "h_vglo2", "h_vglo3")
end

# Convert contravariant momentum components to the global coordinate system
@inline function contravariant2global(u, aux_vars,
                                      equations::CovariantShallowWaterEquations2D)
    vglo1, vglo2, vglo3 = basis_covariant(aux_vars, equations) * momentum(u, equations)
    return SVector(u[1], vglo1, vglo2, vglo3)
end

# Convert momentum components in the global coordinate system to contravariant components
@inline function global2contravariant(u, aux_vars,
                                      equations::CovariantShallowWaterEquations2D)
    vcon1, vcon2 = basis_contravariant(aux_vars, equations) * SVector(u[2], u[3], u[4])
    return SVector(u[1], vcon1, vcon2)
end
# The flux for the covariant form takes in the element container and node/element indices
# in order to give the flux access to the geometric information
@inline function Trixi.flux(u, aux_vars, orientation::Integer,
                            equations::CovariantShallowWaterEquations2D)
    h, h_vcon1, h_vcon2 = u

    J = area_element(aux_vars, equations)
    Gcon = metric_contravariant(aux_vars, equations)
    gravitational_term = 0.5f0 * equations.gravity * h^2
    h_vcon = SVector(h_vcon1, h_vcon2)
    vcon_orientation = h_vcon[orientation] / h

    momentum_flux_1 = h_vcon1 * vcon_orientation +
                      Gcon[1, orientation] * gravitational_term
    momentum_flux_2 = h_vcon2 * vcon_orientation +
                      Gcon[2, orientation] * gravitational_term

    return SVector(J * h_vcon[orientation], J * momentum_flux_1, J * momentum_flux_2)
end

# Symmetric part of entropy-conservative flux
@inline function flux_split_covariant(u_ll, u_rr, aux_vars_ll, aux_vars_rr,
                                      orientation::Integer,
                                      equations::CovariantShallowWaterEquations2D)
    h_ll, h_vcon1_ll, h_vcon2_ll = u_ll
    h_rr, h_vcon1_rr, h_vcon2_rr = u_rr

    J_ll = area_element(aux_vars_ll, equations)
    J_rr = area_element(aux_vars_rr, equations)

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
@inline function flux_nonconservative_split_covariant(u_ll, u_rr, aux_vars_ll,
                                                      aux_vars_rr,
                                                      orientation::Integer,
                                                      equations::CovariantShallowWaterEquations2D)
    h_ll, h_vcon1_ll, h_vcon2_ll = u_ll
    h_rr, h_vcon1_rr, h_vcon2_rr = u_rr

    # Geometric variables
    Gcov_ll = metric_covariant(aux_vars_ll, equations)
    Gcov_rr = metric_covariant(aux_vars_rr, equations)
    Gcon_ll = metric_contravariant(aux_vars_ll, equations)
    J_ll = area_element(aux_vars_ll, equations)
    J_rr = area_element(aux_vars_rr, equations)

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
    nonconservative_term_gravitational = equations.gravity *
                                         J_ll * Gcon_ll[:, orientation] *
                                         h_ll * h_rr # the same as h_ll * (h_ll + h_rr)

    nonconservative_term = nonconservative_term_inertial +
                           nonconservative_term_gravitational

    return SVector(zero(eltype(u_ll)), nonconservative_term[1], nonconservative_term[2])
end

@inline function source_terms_weak_form(u, x, t, aux_vars,
                                        equations::CovariantShallowWaterEquations2D)
    h, h_vcon1, h_vcon2 = u

    # Geometric variables
    Gcon = metric_contravariant(aux_vars, equations)
    (Gamma1, Gamma2) = christoffel_symbols(aux_vars, equations)
    J = area_element(aux_vars, equations)

    # Physical variables
    h_vcon = SVector{2}(h_vcon1, h_vcon2)
    v_con = SVector{2}(h_vcon1 / h, h_vcon2 / h)

    # Doubly-contravariant flux tensor
    T = h_vcon * v_con' + 0.5f0 * equations.gravity * h^2 * Gcon

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
@inline function source_terms_split_covariant(u, x, t, aux_vars,
                                              equations::CovariantShallowWaterEquations2D)
    h, h_vcon1, h_vcon2 = u

    # Geometric variables
    Gcov = metric_covariant(aux_vars, equations)
    Gcon = metric_contravariant(aux_vars, equations)
    (Gamma1, Gamma2) = christoffel_symbols(aux_vars, equations)
    J = area_element(aux_vars, equations)

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
@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, aux_vars_ll, aux_vars_rr,
                                           orientation,
                                           equations::CovariantShallowWaterEquations2D)
    h_ll, h_vcon1_ll, h_vcon2_ll = u_ll
    h_rr, h_vcon1_rr, h_vcon2_rr = u_rr

    h_vcon_ll = SVector(h_vcon1_ll, h_vcon2_ll)
    h_vcon_rr = SVector(h_vcon1_rr, h_vcon2_rr)

    Gcon_ll = metric_contravariant(aux_vars_ll, equations)
    Gcon_rr = metric_contravariant(aux_vars_rr, equations)

    phi_ll = max(h_ll * equations.gravity, 0)
    phi_rr = max(h_rr * equations.gravity, 0)

    return max(abs(h_vcon_ll[orientation] / h_ll) +
               sqrt(Gcon_ll[orientation, orientation] * phi_ll),
               abs(h_vcon_rr[orientation] / h_rr) +
               sqrt(Gcon_rr[orientation, orientation] * phi_rr))
end

# Maximum wave speeds with respect to the covariant basis
@inline function Trixi.max_abs_speeds(u, aux_vars,
                                      equations::CovariantShallowWaterEquations2D)
    h, h_vcon1, h_vcon2 = u
    Gcon = metric_contravariant(aux_vars, equations)
    phi = max(h * equations.gravity, 0)
    return abs(h_vcon1 / h) + sqrt(Gcon[1, 1] * phi),
           abs(h_vcon2 / h) + sqrt(Gcon[2, 2] * phi)
end

# Rossby-Haurwitz wave
function initial_condition_rossby_haurwitz(x, t,
                                           equations::CovariantShallowWaterEquations2D{<:GlobalSphericalCoordinates})
    (; gravity, rotation_rate) = equations

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
        (1 / gravity) *
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
                                                  equations::CovariantShallowWaterEquations2D{GlobalSphericalCoordinates})
    (; gravity) = equations
    realT = eltype(x)
    radius = sqrt(x[1]^2 + x[2]^2 + x[3]^2)
    lon = atan(x[2], x[1])
    lat = asin(x[3] / radius)

    # compute zonal and meridional velocity components
    u_0 = 80.0f0
    lat_0 = convert(realT, π / 7)
    lat_1 = convert(realT, π / 2) - lat_0
    vlon = galewsky_velocity(lat, u_0, lat_0, lat_1)
    vlat = zero(eltype(x))

    # numerically integrate (here we use the QuadGK package) to get height
    galewsky_integral, _ = quadgk(latp -> galewsky_integrand(latp, u_0, lat_0, lat_1,
                                                             radius,
                                                             equations), -π / 2, lat)
    h = 10158.0f0 - radius / gravity * galewsky_integral

    # add perturbation to initiate instability
    α, β = convert(realT, 1 / 3), convert(realT, 1 / 15)
    lat_2 = convert(realT, π / 4)
    if (-π < lon) && (lon < π)
        h = h + 120.0f0 * cos(lat) * exp(-((lon / α)^2)) * exp(-((lat_2 - lat) / β)^2)
    end
    # convert to conservative variables
    return SVector(h, h * vlon, h * vlat)
end

# If the initial velocity field is defined in Cartesian coordinates and the chosen global 
# coordinate system is spherical, perform the appropriate conversion
@inline function cartesian2global(u, x,
    ::CovariantShallowWaterEquations2D{GlobalSphericalCoordinates})
    h_vlon, h_vlat, h_vrad = cartesian2spherical(u[2], u[3], u[4], x)
return SVector(u[1], h_vlon, h_vlat, h_vrad)
end

# If the initial velocity field is defined in spherical coordinates and the chosen global 
# coordinate system is Cartesian, perform the appropriate conversion
@inline function spherical2global(u, x,
    ::CovariantShallowWaterEquations2D{GlobalCartesianCoordinates})
h_vx, h_vy, h_vz = spherical2cartesian(u[2], u[3], u[4], x)
return SVector(u[1], h_vx, h_vy, h_vz)
end

# If the initial velocity field is defined in spherical coordinates and the chosen global 
# coordinate system is spherical, do not convert
@inline function spherical2global(u, x,
    ::CovariantShallowWaterEquations2D{GlobalSphericalCoordinates})
return u
end
end # @muladd
