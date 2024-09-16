@muladd begin
#! format: noindent

# Variable-coefficient linear advection equation in covariant form on a two-dimensional
# surface in three-dimensional space
struct CovariantShallowWaterEquations2D{RealT <: Real} <: AbstractCovariantEquations{2,
                                  3, 3}
    gravitational_acceleration::RealT
    rotation_rate::RealT
end

function Trixi.varnames(::typeof(cons2cons), ::CovariantShallowWaterEquations2D)
    return ("h", "hv_con_1", "hv_con_2")
end

# TODO: actual entropy variables
Trixi.cons2entropy(u, ::CovariantShallowWaterEquations2D) = SVector{3}(u[1],
                                                                       zero(eltype(u)),
                                                                       zero(eltype(u)))

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

    v = u[1 + orientation] / h
    T_con_1 = hv_con_1 * v + G_con_1 * gravitational_term
    T_con_2 = hv_con_2 * v + G_con_2 * gravitational_term

    return volume_element(elements, i, j, element) * SVector(hv_con_1, T_con_1, T_con_2)
end

# Directional flux that takes in the normal components in reference space
@inline function Trixi.flux(u, normal_direction::AbstractVector,
                            ::CovariantShallowWaterEquations2D,
                            elements, i, j, element)
    h, hv_con_1, hv_con_2 = u

    gravitational_term = 0.5f0 * equations.gravitational_acceleration * h^2

    G_con_1 = elements.contravariant_metric[1, 1, i, j, element] * normal_direction[1] +
              elements.contravariant_metric[1, 2, i, j, element] * normal_direction[2]
    G_con_2 = elements.contravariant_metric[2, 1, i, j, element] * normal_direction[1] +
              elements.contravariant_metric[2, 2, i, j, element] * normal_direction[2]

    v = (hv_con_1 * normal_direction[1] + hv_con_2 * normal_direction[2]) / h
    T_con_1 = hv_con_1 * v + G_con_1 * gravitational_term
    T_con_2 = hv_con_2 * v + G_con_2 * gravitational_term

    return volume_element(elements, i, j, element) * SVector(hv_con_1, T_con_1, T_con_2)
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
@inline function source_terms_spherical(u, x, t,
                                        equations::CovariantShallowWaterEquations2D,
                                        elements, i, j, element)
    h, hv_con_1, hv_con_2 = u

    # Geometric variables
    G_con = SMatrix{2, 2}(view(elements.contravariant_metric, :, :, i, j, element))
    Gamma1 = SMatrix{2, 2}(view(elements.christoffel_symbols, :, :, 1, i, j, element))
    Gamma2 = SMatrix{2, 2}(view(elements.christoffel_symbols, :, :, 2, i, j, element))
    J = volume_element(elements, i, j, element)

    # Physical variables
    hv_con = SVector{2}(hv_con_1, hv_con_2)
    v_con = SVector{2}(hv_con_1 / h, hv_con_2 / h)
    T = hv_con * v_con' + 0.5f0 * equations.gravitational_acceleration * h^2 * G_con

    # Coriolis parameter
    f = 2 * equations.rotation_rate * x[3] / norm(x)  # 2Ωsinθ

    # Combined source terms
    source_1 = sum(Gamma1 .* T) +
               f * J * (G_con[1, 2] * hv_con_1 - G_con[1, 1] * hv_con_2)
    source_2 = sum(Gamma2 .* T) +
               f * J * (G_con[2, 2] * hv_con_1 - G_con[2, 1] * hv_con_2)

    # Scale by negative Jacobian
    return SVector(zero(eltype(u)), -J * source_1, -J * source_2)
end

# Maximum wave speed along the normal direction in reference space
@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, normal_direction,
                                           ::CovariantShallowWaterEquations2D,
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
@inline function Trixi.max_abs_speeds(u, ::CovariantShallowWaterEquations2D,
                                      elements, i, j, element)
    h, hv_con_1, hv_con_2 = u

    G_con_11 = elements.contravariant_metric[1, 1, i, j, element]
    G_con_22 = elements.contravariant_metric[2, 2, i, j, element]
    v_con_1 = hv_con_1 / h
    v_con_2 = hv_con_2 / h
    phi = max(h_ll * equations.gravitational_acceleration, 0)
    return abs(v_con_1) + sqrt(G_con_11 * phi), abs(v_con_2) + sqrt(G_con_22 * phi)
end
end # @muladd
