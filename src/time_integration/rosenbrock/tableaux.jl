# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Coefficients of a Rosenbrock-W method in the transformed form of Hairer and Wanner (1996),
# Section IV.7, Eq. (7.25): `a = alpha * inv(Gamma)`, `c = inv(diag(Gamma)) - inv(Gamma)`,
# `m = transpose(b) * inv(Gamma)`, the stage times `c_stage` (the row sums of `alpha`), and
# the constant diagonal `gamma` of `Gamma`. Both `a` and `c` are strictly lower triangular.
struct RosenbrockTableau{RealT <: Real, MatrixT, VectorT}
    a::MatrixT
    c::MatrixT
    m::VectorT
    c_stage::VectorT
    gamma::RealT
end

# Transformation of the coefficients `alpha`, `b` (explicit part), and `Gamma` (implicit
# part, lower triangular with constant diagonal) of a Rosenbrock-W method
function RosenbrockTableau(alpha, b, Gamma)
    gamma = Gamma[1, 1]
    inv_Gamma = inv(Gamma)
    a = alpha * inv_Gamma
    c = inv(gamma) * I - inv_Gamma
    m = transpose(inv_Gamma) * b
    c_stage = vec(sum(alpha, dims = 2))

    # The entries on and above the diagonal of `a` and `c` are set to zero exactly (instead
    # of round-off errors), since the stage combination relies on them to skip the stages
    # not computed yet
    a = a - UpperTriangular(a)
    c = c - UpperTriangular(c)

    return RosenbrockTableau{typeof(gamma), typeof(a), typeof(m)}(a, c, m, c_stage,
                                                                  gamma)
end

"""
    SSPRW32()

Three-stage, second-order, L-stable Rosenbrock-W method, whose explicit part
is the SSPRK(3,3) method of Shu and Osher.

The approximate Jacobian `J` of the right-hand side is passed as keyword argument
`jacobian` to `Trixi.solve`, e.g.,
`Trixi.solve(ode, SSPRW32(); dt, jacobian = VerticalJacobianTwoSided())`, and its arrays
are created by `Trixi.init` for the semidiscretization of the ODE problem. The default
`jacobian = nothing` corresponds to `J = 0`, i.e., the explicit SSPRK(3,3) method.

## References

- Marco Artiano, Andrés Rueda-Ramírez, Oswald Knoth, Hendrik Ranocha, Gregor Gassner (2026)
  GPU-Accelerated Nested Schur-Complement Solvers for Atmospheric Flow Simulations using
  Horizontally-Explicit, Vertically-Implicit Schemes and Discontinuous Galerkin Methods
  [arXiv: 2610.15775](https://arxiv.org/abs/2610.15775)
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
