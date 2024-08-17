###############################################################################
# Shallow water equations in covariant form
###############################################################################

@muladd begin
#! format: noindent

struct CovariantShallowWaterEquations2D{RealT <: Real} <:
       AbstractCovariantEquations2D{6}
    gravitational_acceleration::RealT
    rotation_rate::RealT
end

function Trixi.varnames(::typeof(cons2cons), ::CovariantShallowWaterEquations2D)
    return ("h", "M_con_1", "M_con_2", "G_cov_11", "G_cov_12", "G_cov_22")
end

function Trixi.cons2entropy(u, equations::CovariantShallowWaterEquations2D)
    z = zero(eltype(u))
    h, M_con_1, M_con_2, G_cov_11, G_cov_12, G_cov_22 = u
    (; gravitational_acceleration) = equations

    v_con_1 = M_con_1 / h
    v_con_2 = M_con_2 / h

    v_cov_1 = G_cov_11 * v_con_1 + G_cov_12 * v_con_2
    v_cov_2 = G_cov_12 * v_con_1 + G_cov_22 * v_con_2

    w_1 = gravitational_acceleration * h -
          0.5 * (v_cov_1 * v_con_1 + v_cov_2 * v_con_2)

    return SVector(w_1, v_cov_1, v_cov_2, z, z, z)
end

"""
    Custom dissipation to ensure no flux is applied to metric coefficients
"""
@inline function (dissipation::DissipationLocalLaxFriedrichs)(u_ll, u_rr,
                                                              normal_direction::AbstractVector,
                                                              equations::CovariantShallowWaterEquations2D)
    z = zero(eltype(u_ll))
    λ = dissipation.max_abs_speed(u_ll, u_rr, normal_direction, equations)
    return -0.5f0 * λ *
           SVector(u_rr[1] - u_ll[1],
                   u_rr[2] - u_ll[2],
                   u_rr[3] - u_ll[3],
                   z, z, z)
end

"""
    Compute a given contravariant flux component
"""
@inline function Trixi.flux(u, orientation::Integer,
                            equations::CovariantShallowWaterEquations2D)
    z = zero(eltype(u))

    h, M_con_1, M_con_2, G_cov_11, G_cov_12, G_cov_22 = u

    half_g_h_squared = 0.5f0 * equations.gravitational_acceleration * h^2

    G = G_cov_11 * G_cov_22 - G_cov_12^2
    G_con_11 = G_cov_22 / G
    G_con_12 = -G_cov_12 / G
    G_con_22 = G_cov_11 / G

    if orientation == 1
        M_con_j = M_con_1
        T_con_1j = (M_con_1 * M_con_1 / h) + G_con_11 * half_g_h_squared
        T_con_2j = (M_con_2 * M_con_1 / h) + G_con_12 * half_g_h_squared
    else # orientation == 2
        M_con_j = M_con_2
        T_con_1j = (M_con_2 * M_con_1 / h) + G_con_12 * half_g_h_squared
        T_con_2j = (M_con_2 * M_con_2 / h) + G_con_22 * half_g_h_squared
    end

    return SVector(M_con_j, T_con_1j, T_con_2j, z, z, z)
end

"""
    Compute the flux component in a given normal direction
"""
@inline function Trixi.flux(u, normal_direction::AbstractVector,
                            equations::CovariantShallowWaterEquations2D)
    return flux(u, 1, equations) * normal_direction[1] +
           flux(u, 2, equations) * normal_direction[2]
end

@inline function source_terms_coriolis_sphere(u, x, t,
                                              equations::CovariantShallowWaterEquations2D)
    _, M_con_1, M_con_2, G_cov_11, G_cov_12, G_cov_22 = u

    z = zero(eltype(u))

    G = G_cov_11 * G_cov_22 - G_cov_12^2
    G_con_11 = G_cov_22 / G
    G_con_12 = -G_cov_12 / G
    G_con_22 = G_cov_11 / G

    fG = 2 * equations.rotation_rate * x[3] / norm(x) * G

    return SVector(z, fG * (G_con_12 * M_con_1 - G_con_11 * M_con_2),
                   fG * (G_con_22 * M_con_1 - G_con_12 * M_con_2),
                   z, z, z)
end

@inline function flux_nonconservative_naive(u_ll, u_rr, orientation::Integer,
                                            equations::CovariantShallowWaterEquations2D)
    _, _, _, G_cov_11_rr, G_cov_12_rr, G_cov_22_rr = u_rr

    z = zero(eltype(u_ll))

    _, T_con_11_ll, _ = flux(u_ll, orientation, equations)
    _, T_con_12_ll, T_con_22_ll = flux(u_ll, orientation, equations)

    G_ll = G_cov_11_ll * G_cov_22_ll - G_cov_12_ll^2

    if orientation == 1
        a11_ll = T_con_11_ll * G_cov_22_ll / G_ll
        a12_ll = 2 * T_con_11_ll * G_cov_12_ll / G_ll
        a13_ll = (T_con_11_ll * G_cov_11 - 2 * T_con_12_ll * G_cov_12_ll -
                  T_con_22_ll * G_cov_22) / (2 * G_ll)

        a21_ll = (T_con_12_ll * G_cov_22_ll - T_con_11_ll * G_cov_12_ll) / (2 * G_ll)
        a22_ll = (T_con_11_ll * G_cov_11_ll - T_con_12_ll * G_cov_12_ll) / G_ll
        a23_ll = (T_con_22_ll * G_cov_12 + 3 * T_con_12 * G_cov_11) / (2 * G_ll)

    else # orientation == 2
        a11_ll = (T_con_11_ll * G_cov_12_ll + 3 * T_con_12_ll * G_cov_22_ll) /
                 (2 * G_ll)
        a12_ll = (T_con_22_ll * G_cov_22_ll - T_con_12_ll * G_cov_12_ll) / G_ll
        a13_ll = (T_con_12_ll * G_cov_11_ll - T_con_22_ll * G_cov_12_ll) / (2 * G_ll)

        a21_ll = (T_con_22_ll * G_cov_22_ll - T_con_11 * G_cov_11 -
                  2 * T_con_12_ll * G_cov_12_ll) / (2 * G_ll)
        a22_ll = -2 * T_con_22_ll * G_cov_12_ll / G_ll
        a23_ll = T_con_22_ll * G_cov_11_ll / G_ll
    end

    return SVector(z,
                   a11_ll * G_cov_11_rr + a12_ll * G_cov_12_rr + a13_ll * G_cov_22_rr,
                   a21_ll * G_cov_21_rr + a22_ll * G_cov_12_rr + a23_ll * G_cov_22_rr,
                   z, z, z)
end

"""
    Maximum directional wave speed for Lax-Friedrichs dissipation
"""
@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, normal_direction,
                                           ::CovariantShallowWaterEquations2D)
    h_ll, M_con_1_ll, M_con_2_ll, G_cov_11, G_cov_12, G_cov_22 = u_ll
    h_rr, M_con_1_rr, M_con_2_rr = u_rr

    (n_con_1, n_con_2) = normal_direction

    G = G_cov_11 * G_cov_22 - G_cov_12^2
    G_con_11 = G_cov_22 / G
    G_con_12 = -G_cov_12 / G
    G_con_22 = G_cov_11 / G

    G_nn = n_con_1 * G_con_11 * n_con_1 + n_con_1 * G_con_12 * n_con_2 +
           n_con_2 * G_con_12 * n_con_1 + n_con_2 * G_con_22 * n_con_2

    v_n_ll = (M_con_1_ll * n_con_1 + M_con_2_ll * n_con_2) / h_ll
    v_n_rr = (M_con_1_rr * n_con_1 + M_con_2_rr * n_con_2) / h_rr

    gh_ll = max(h_ll * equations.gravitational_acceleration, 0)
    gh_rr = max(h_rr * equations.gravitational_acceleration, 0)

    return max(abs(v_n_ll) + sqrt(G_nn * gh_ll), abs(v_n_rr) + sqrt(G_nn * gh_rr))
end

"""
    Maximum wave speeds with respect to the contravariant basis
"""
@inline function Trixi.max_abs_speeds(u, ::CovariantShallowWaterEquations2D)
    h, M_con_1, M_con_2, G_cov_11, G_cov_12, G_cov_22 = u

    G = G_cov_11 * G_cov_22 - G_cov_12^2
    G_con_11 = G_cov_22 / G
    G_con_22 = G_cov_11 / G

    gh = max(h * equations.gravitational_acceleration, 0)

    v_con_1 = M_con_1 / h
    v_con_2 = M_con_2 / h

    return abs(v_con_1) + sqrt(G_con_11 * gh), abs(v_con_2) + sqrt(G_con_22 * gh)
end

@inline function cartesian2contravariant(u_cartesian,
                                         ::CovariantShallowWaterEquations2D,
                                         i, j, element, cache)
    (; contravariant_vectors, inverse_jacobian) = cache.elements

    A11, A21, A31, A12, A22, A32 = view(cache.elements.jacobian_matrix, :, :, i, j,
                                        element)

    G_cov_11 = A11 * A11 + A21 * A21 + A31 * A31
    G_cov_12 = A11 * A12 + A21 * A22 + A31 * A32
    G_cov_22 = A12 * A12 + A22 * A22 + A32 * A32

    Ja11, Ja12, Ja13 = Trixi.get_contravariant_vector(1, contravariant_vectors, i, j,
                                                      element)
    Ja21, Ja22, Ja23 = Trixi.get_contravariant_vector(2, contravariant_vectors, i, j,
                                                      element)
    return SVector(u_cartesian[1],
                   inverse_jacobian[i, j, element] * (Ja11 * u_cartesian[2] +
                    Ja12 * u_cartesian[3] +
                    Ja13 * u_cartesian[4]),
                   inverse_jacobian[i, j, element] * (Ja21 * u_cartesian[2] +
                    Ja22 * u_cartesian[3] +
                    Ja23 * u_cartesian[4]),
                   G_cov_11, G_cov_12, G_cov_22)
end

@inline function contravariant2cartesian(u_node, ::CovariantShallowWaterEquations2D,
                                         i, j, element, cache)
    A11, A21, A31, A12, A22, A32 = view(cache.elements.jacobian_matrix, :, :, i, j,
                                        element)
    return SVector(u_node[1],
                   A11 * u_node[2] + A12 * u_node[3],
                   A21 * u_node[2] + A22 * u_node[3],
                   A31 * u_node[2] + A32 * u_node[3])
end
end # @muladd
