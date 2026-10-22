# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Coefficients of a linear IMEX Runge-Kutta method: the explicit coefficients `a` and `b`,
# the stage times `c` (the row sums of `a`), `d = (a_im - a) / gamma` of the auxiliary
# variables of the stages, and the constant diagonal `gamma` of the implicit coefficients
# `a_im`. Both `a` and `d` are strictly lower triangular.
struct LinearIMEXTableau{RealT <: Real, MatrixT, VectorT}
    a::MatrixT
    d::MatrixT
    b::VectorT
    c::VectorT
    gamma::RealT
end

# Coefficients of the explicit part `a`, `b` and of the implicit part `a_im` (lower
# triangular with an explicit first stage and a constant diagonal of the other stages) of
# an IMEX Runge-Kutta method. The final update uses only the right-hand sides of the stages,
# so the weights of both parts must agree.
function LinearIMEXTableau(a, b, a_im)
    gamma = a_im[2, 2]
    d = (a_im - a) / gamma
    c = vec(sum(a, dims = 2))

    # The entries on and above the diagonal of `d` are set to zero exactly (instead of
    # round-off errors), since the stage combination relies on them to skip the stages not
    # computed yet
    d = d - UpperTriangular(d)

    return LinearIMEXTableau{typeof(gamma), typeof(a), typeof(b)}(a, d, b, c, gamma)
end

"""
    ARK2()

A second-order, three-stage type II IMEX method developed by Giraldo, Kelly, and
Constantinescu (2013). The implicit part is L-stable and stiffly accurate.

## References

- Francis X. Giraldo, James F. Kelly, Emil M. Constantinescu (2013)
  Implicit-Explicit Formulations of a Three-Dimensional Nonhydrostatic Unified Model of
  the Atmosphere (NUMA)
  [DOI: 10.1137/120876034](https://doi.org/10.1137/120876034)
"""
struct ARK2 <: AbstractLinearIMEXAlgorithm{3} end

function tableau(::ARK2, RealT)
    # ARK2 with the free parameter `a32` of Giraldo, Kelly, and Constantinescu (2013)
    sqrt2 = sqrt(convert(RealT, 2))
    gamma = 1 - 1 / sqrt2
    a32 = (3 + 2 * sqrt2) / 6
    a = SMatrix{3, 3, RealT}(0, 2 * gamma, 1 - a32,
                             0, 0, a32,
                             0, 0, 0)
    b = SVector{3, RealT}(1 / (2 * sqrt2), 1 / (2 * sqrt2), gamma)

    a_im = SMatrix{3, 3, RealT}(0, gamma, 1 / (2 * sqrt2),
                                0, gamma, 1 / (2 * sqrt2),
                                0, 0, gamma)

    return LinearIMEXTableau(a, b, a_im)
end

"""
    ARS232()

A second-order, three-stage type II IMEX method developed by Ascher, Ruuth, and Spiteri
(1997). The implicit part is L-stable and stiffly accurate.

## References

- Uri M. Ascher, Steven J. Ruuth, Raymond J. Spiteri (1997)
  Implicit-explicit Runge-Kutta methods for time-dependent partial differential equations
  [DOI: 10.1016/S0168-9274(97)00056-1](https://doi.org/10.1016/S0168-9274(97)00056-1)
"""
struct ARS232 <: AbstractLinearIMEXAlgorithm{3} end

function tableau(::ARS232, RealT)
    # ARS(2,3,2) IMEX Runge-Kutta - Second order
    sqrt2 = sqrt(convert(RealT, 2))
    gamma = 1 - 1 / sqrt2
    delta = -2 * sqrt2 / 3
    a = SMatrix{3, 3, RealT}(0, gamma, delta,
                             0, 0, 1 - delta,
                             0, 0, 0)
    b = SVector{3, RealT}(0, 1 - gamma, gamma)

    a_im = SMatrix{3, 3, RealT}(0, 0, 0,
                                0, gamma, 1 - gamma,
                                0, 0, gamma)

    return LinearIMEXTableau(a, b, a_im)
end
"""
    ARS343()

A third-order, four-stage type II IMEX method developed by Ascher, Ruuth, and Spiteri
(1997). The implicit part is L-stable and stiffly accurate.

## References

- Uri M. Ascher, Steven J. Ruuth, Raymond J. Spiteri (1997)
  Implicit-explicit Runge-Kutta methods for time-dependent partial differential equations
  [DOI: 10.1016/S0168-9274(97)00056-1](https://doi.org/10.1016/S0168-9274(97)00056-1)
"""
struct ARS343 <: AbstractLinearIMEXAlgorithm{4} end

function tableau(::ARS343, RealT)
    # ARS(3,4,3) IMEX Runge-Kutta - Third order
    # `gamma` is the middle root of `x^3 - 3 x^2 + 3/2 x - 1/6 = 0`, computed in `RealT`
    # from the trigonometric solution of the cubic equation
    sqrt2 = sqrt(convert(RealT, 2))
    gamma = 1 + sqrt2 * cos((acos(2 * sqrt2 / 3) - 2 * convert(RealT, pi)) / 3)
    c3 = (1 + gamma) / 2
    b2 = -3 * gamma^2 / 2 + 4 * gamma - 1 / convert(RealT, 4)
    b3 = 3 * gamma^2 / 2 - 5 * gamma + 5 / convert(RealT, 4)

    # The explicit coefficients are given to ten digits by Ascher, Ruuth, and Spiteri
    # (1997), Section 2.7. Only `a32` is taken from there: with `a42 = a43`, the remaining
    # ones follow from the row sums `c` and the third-order condition `b' * a * c = 1/6`,
    # which are thus satisfied to the precision of `RealT`.
    a32 = convert(RealT, 0.3966543747)
    a31 = c3 - a32
    a42 = (1 / convert(RealT, 6) - b3 * a32 * gamma) / (gamma * (gamma + c3))
    a43 = a42
    a41 = 1 - a42 - a43
    a = SMatrix{4, 4, RealT}(0, gamma, a31, a41,
                             0, 0, a32, a42,
                             0, 0, 0, a43,
                             0, 0, 0, 0)
    b = SVector{4, RealT}(0, b2, b3, gamma)

    a_im = SMatrix{4, 4, RealT}(0, 0, 0, 0,
                                0, gamma, (1 - gamma) / 2, b2,
                                0, 0, gamma, b3,
                                0, 0, 0, gamma)

    return LinearIMEXTableau(a, b, a_im)
end
end # @muladd
