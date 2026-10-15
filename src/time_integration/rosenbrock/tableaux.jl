# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

"""
    SSPRW32(; jacobian = nothing)

Three-stage, second-order, L-stable Rosenbrock-W method ([`AbstractRosenbrockAlgorithm`](@ref))
for horizontally explicit, vertically implicit (HEVI) time integration, whose explicit
part is the SSPRK(3,3) method of Shu and Osher, with the approximate `jacobian` of the
vertical dynamics, e.g., [`VerticalJacobianTwoSided`](@ref). With `jacobian = nothing`,
the method reduces to the explicit SSPRK(3,3) method. The time integrator is used like
the ones of Trixi.jl, e.g.,
`Trixi.solve(ode, SSPRW32(jacobian = jacobian); dt = 85.0, callback = callbacks)`.

## References

- Marco Artiano, Andrés Rueda-Ramírez, Oswald Knoth, Hendrik Ranocha, Gregor Gassner (2026)
  GPU-Accelerated Nested Schur-Complement Solvers for Atmospheric Flow Simulations using
  Horizontally-Explicit, Vertically-Implicit Schemes and Discontinuous Galerkin Methods
  [arXiv:2610.15775](https://arxiv.org/abs/2610.15775)
- Chi-Wang Shu, Stanley Osher (1988)
  Efficient Implementation of Essentially Non-oscillatory Shock-Capturing Schemes (Eq. 2.18)
  [DOI: 10.1016/0021-9991(88)90177-5](https://doi.org/10.1016/0021-9991(88)90177-5)
"""
struct SSPRW32{Jacobian} <: AbstractRosenbrockAlgorithm{3}
    # Coefficients of the transformed form: `a = alpha * inv(Gamma)`,
    # `c = inv(diag(Gamma)) - inv(Gamma)`, `m = transpose(b) * inv(Gamma)`. Both `a` and `c`
    # are strictly lower triangular.
    a::SMatrix{3, 3, Float64, 9}
    c::SMatrix{3, 3, Float64, 9}
    m::SVector{3, Float64}
    # Stage times `t + c_stage * dt`, i.e., the row sums of `alpha`
    c_stage::SVector{3, Float64}
    # Constant diagonal of `Gamma`
    gamma::Float64
    jacobian::Jacobian

    function SSPRW32(; jacobian = nothing)
        # Explicit part: SSPRK(3,3)
        alpha = @SMatrix [0.0 0.0 0.0;
                          1.0 0.0 0.0;
                          0.25 0.25 0.0]
        b = SVector(1 / 6, 1 / 6, 2 / 3)

        # Implicit part. `gamma21 = -3 gamma` cancels the third-order W-residual
        # `e1 = sum b_i alpha_ij gamma_j`.
        gamma = 1 - 1 / sqrt(2)
        gamma21 = -3 * gamma
        gamma31 = (3 * sqrt(2) - 4) / 8
        gamma32 = -1 / 4
        Gamma = @SMatrix [gamma 0.0 0.0;
                          gamma21 gamma 0.0;
                          gamma31 gamma32 gamma]

        inv_Gamma = inv(Gamma)
        a = alpha * inv_Gamma
        c = inv(gamma) * I - inv_Gamma
        m = transpose(inv_Gamma) * b
        c_stage = alpha * SVector(1.0, 1.0, 1.0)

        # Both `a` and `c` are strictly lower triangular. The entries on and above the
        # diagonal are set to zero exactly (instead of round-off errors), since the stage
        # combination relies on them to skip the stages not computed yet.
        a = SMatrix{3, 3}(j < i ? a[i, j] : 0.0 for i in 1:3, j in 1:3)
        c = SMatrix{3, 3}(j < i ? c[i, j] : 0.0 for i in 1:3, j in 1:3)

        return new{typeof(jacobian)}(a, c, m, c_stage, gamma, jacobian)
    end
end
end # @muladd
