# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Coefficients of a Rosenbrock-W method in the transformed form of Hairer and Wanner (1996),
# Section IV.7, eq. (7.25): `a = alpha * inv(Gamma)`, `c = inv(diag(Gamma)) - inv(Gamma)`,
# `m = transpose(b) * inv(Gamma)`, the stage times `c_stage` (the row sums of `alpha`), and
# the constant diagonal `gamma` of `Gamma`. Both `a` and `c` are strictly lower triangular.
struct RosenbrockTableau{N, RealT <: Real, N2}
    a::SMatrix{N, N, RealT, N2}
    c::SMatrix{N, N, RealT, N2}
    m::SVector{N, RealT}
    c_stage::SVector{N, RealT}
    gamma::RealT
end

# Transformation of the coefficients `alpha`, `b` (explicit part), and `Gamma` (implicit
# part, lower triangular with constant diagonal) of a Rosenbrock-W method
function RosenbrockTableau(alpha::SMatrix{N, N, RealT}, b::SVector{N, RealT},
                           Gamma::SMatrix{N, N, RealT}) where {N, RealT}
    gamma = Gamma[1, 1]
    inv_Gamma = inv(Gamma)
    a = alpha * inv_Gamma
    c = inv(gamma) * I - inv_Gamma
    m = transpose(inv_Gamma) * b
    c_stage = alpha * ones(SVector{N, RealT})

    # The entries on and above the diagonal of `a` and `c` are set to zero exactly (instead
    # of round-off errors), since the stage combination relies on them to skip the stages
    # not computed yet
    a = SMatrix{N, N, RealT}(j < i ? a[i, j] : 0 for i in 1:N, j in 1:N)
    c = SMatrix{N, N, RealT}(j < i ? c[i, j] : 0 for i in 1:N, j in 1:N)

    return RosenbrockTableau{N, RealT, N * N}(a, c, m, c_stage, gamma)
end

"""
    SSPRW32()

Three-stage, second-order, L-stable Rosenbrock-W method for horizontally explicit,
vertically implicit (HEVI) time integration, whose explicit part is the SSPRK(3,3) method of
Shu and Osher. The approximate Jacobian (e.g., [`VerticalJacobianTwoSided`](@ref)) is passed
to `Trixi.solve` as keyword argument `jacobian`.

## References

- Marco Artiano, Andrés Rueda-Ramírez, Oswald Knoth, Hendrik Ranocha, Gregor Gassner (2026)
  GPU-Accelerated Nested Schur-Complement Solvers for Atmospheric Flow Simulations using
  Horizontally-Explicit, Vertically-Implicit Schemes and Discontinuous Galerkin Methods
  [arXiv:2610.15775](https://arxiv.org/abs/2610.15775)
- Chi-Wang Shu, Stanley Osher (1988)
  Efficient Implementation of Essentially Non-oscillatory Shock-Capturing Schemes (Eq. 2.18)
  [DOI: 10.1016/0021-9991(88)90177-5](https://doi.org/10.1016/0021-9991(88)90177-5)
"""
struct SSPRW32 <: AbstractRosenbrockAlgorithm{3} end

function tableau(::SSPRW32, RealT)
    # Explicit part: SSPRK(3,3)
    alpha = SMatrix{3, 3, RealT}(0, 1, 1 // 4,
                                 0, 0, 1 // 4,
                                 0, 0, 0)
    b = SVector{3, RealT}(1 // 6, 1 // 6, 2 // 3)

    # Implicit part. `gamma21 = -3 gamma` cancels the third-order W-residual
    # `e1 = sum b_i alpha_ij gamma_j`.
    sqrt2 = sqrt(convert(RealT, 2))
    gamma = 1 - 1 / sqrt2
    gamma21 = -3 * gamma
    gamma31 = (3 * sqrt2 - 4) / 8
    gamma32 = convert(RealT, -1 // 4)
    Gamma = SMatrix{3, 3, RealT}(gamma, gamma21, gamma31,
                                 0, gamma, gamma32,
                                 0, 0, gamma)

    return RosenbrockTableau(alpha, b, Gamma)
end
end # @muladd
