@muladd begin
#! format: noindent

@doc raw"""
    DeardorffTKE(; c_k = 0.125)

Vertical model of the turbulent kinetic energy (TKE) of Deardorff (1980) in the form of
CGDycore.jl, used as vertical model of a [`HorizontalVerticalTurbulence`](@ref) whose
horizontal model has no prognostic variables. The prognostic variable is the density of the
TKE ``\varrho e``. The eddy viscosity and diffusivity are
```math
\nu_v = \kappa_v = c_k \Delta \sqrt{e},
```
with the vertical grid spacing ``\Delta``, i.e., the height of the element divided by the
number of nodes in the vertical, also for the TKE itself. The TKE equation is
```math
\partial_t (\varrho e) + \nabla \cdot (\varrho e \vec{v})
  = \nabla \cdot (\varrho D \nabla e)
  + \varrho \left(2 \nu_v \lVert S - P S P \rVert^2 - \kappa_v N^2 - \varepsilon\right),
\qquad
\varepsilon = C_\varepsilon \frac{e^{3/2}}{l},
```
with the shear production, which is the dissipation of the vertical part of the turbulent
stress of [`TurbulentDiffusion3D`](@ref), the buoyancy production with the squared
Brunt-Väisälä frequency ``N^2 = \lVert \nabla \phi \rVert \hat{k} \cdot \nabla \theta / \theta``,
the length scale ``l = \Delta`` for ``N^2 \leq 0`` and
``l = \min(\Delta, 0.76 \sqrt{e} / N)`` for ``N^2 > 0``, and
``C_\varepsilon = 1.9 c_k + (0.93 - 1.9 c_k) l / \Delta``. The dissipation heats the air,
``\partial_t (\varrho \theta) = \varrho \varepsilon / (c_p \pi)`` with the Exner pressure
``\pi = (p / p_0)^{R / c_p}``, such that the total energy is conserved.

The source terms are not part of the right-hand side. They are integrated exactly by
[`turbulence_sources_exact!`](@ref), passed as `source_splitting!` to `Trixi.solve`, with the
gradients and all other variables frozen: then, ``s = \sqrt{e}`` satisfies
```math
\frac{\mathrm{d} s}{\mathrm{d} t} = \frac{A - a s - b s^2}{2}, \qquad
A = c_k \Delta \left(2 \lVert S - P S P \rVert^2 - N^2\right),
```
with ``a = 0`` and ``b = 0.93 / \Delta`` for ``l = \Delta``, and
``a = 1.9 c_k N / 0.76`` and ``b = (0.93 - 1.9 c_k) / \Delta`` for ``l < \Delta``, whose
solution is known in closed form. Since ``s = 0`` is absorbing for ``A \leq 0``, the TKE stays
nonnegative. The density of the TKE ``\varrho e``, see
[`density_turbulent_kinetic_energy`](@ref), is kept nonnegative during the transport by the
`PositivityPreservingLimiterZhangShu` of Trixi.jl as `stage_limiter!` and `step_limiter!`.
The limiter requires nonnegative mean values in the elements, which the turbulent diffusion
of ``\varrho e`` with the parabolic scheme of Bassi and Rebay does not guarantee, so small
negative values of the TKE can remain, for which the eddy viscosity and the source terms use
``\max(e, 0)``. The model is available for the
[`CompressibleEulerPotentialTemperatureEquationsWithGravity3D`](@ref).

## References

- James W. Deardorff (1980)
  Stratocumulus-capped mixed layers derived from a three-dimensional model
  [DOI: 10.1007/BF00119502](https://doi.org/10.1007/BF00119502)
"""
struct DeardorffTKE{RealT <: Real}
    c_k::RealT # coefficient of the eddy viscosity
end

DeardorffTKE(; c_k = 0.125) = DeardorffTKE(c_k)

function Base.similar(turbulence::DeardorffTKE, ::Type{NewRealT}) where {NewRealT}
    return DeardorffTKE(convert(NewRealT, turbulence.c_k))
end

@inline nturbulencevariables(::DeardorffTKE) = 1

varnames_turbulence(::typeof(cons2cons), ::DeardorffTKE) = ("rho_e",)
varnames_turbulence(::typeof(cons2prim), ::DeardorffTKE) = ("e",)

@inline have_turbulence_sources(::DeardorffTKE) = true

"""
    density_turbulent_kinetic_energy(u, equations)

Density of the turbulent kinetic energy of the [`DeardorffTKE`](@ref), i.e., the last
conservative variable, e.g., as variable of the `PositivityPreservingLimiterZhangShu` of
Trixi.jl.
"""
@inline density_turbulent_kinetic_energy(u, equations) = last(u)

# Squared Brunt-Väisälä frequency of the specific variables `w` with the potential
# temperature `w[5]`
@inline function brunt_vaisala_frequency_squared(w, gradients, k, aux)
    gradients_x, gradients_y, gradients_z = gradients
    k_dot_gradient_theta = k[1] * gradients_x[5] + k[2] * gradients_y[5] +
                           k[3] * gradients_z[5]
    return aux[4] * k_dot_gradient_theta / w[5]
end

# Squared norm `||S - P S P||^2` of the vertical part of the rate-of-strain tensor `S`, with
# `S - P S P = k (S k)' + (S k) k' - (k' S k) k k'`
@inline function strain_rate_vertical_squared(gradients, k)
    gradients_x, gradients_y, gradients_z = gradients
    velocity_gradient = SMatrix{3, 3}(gradients_x[2], gradients_x[3], gradients_x[4],
                                      gradients_y[2], gradients_y[3], gradients_y[4],
                                      gradients_z[2], gradients_z[3], gradients_z[4])
    S = 0.5f0 * (velocity_gradient + velocity_gradient')
    S_k = S * k
    k_S_k = dot(k, S_k)
    return sum(abs2, k * S_k' + S_k * k' - k_S_k * (k * k'))
end

# Dissipation `C_epsilon e^(3/2) / l` of the TKE `e`, written with
# `l / Delta = min(1, 0.76 sqrt(e) / (N Delta))` and
# `e^(3/2) / l = max(e^(3/2) / Delta, e N / 0.76)` for `N^2 > 0`, such that it vanishes for
# `e = 0` without dividing by `l = 0`
@inline function dissipation_tke(e, N2, grid_spacing, turbulence::DeardorffTKE)
    RealT = typeof(e)
    N = sqrt(max(N2, zero(N2)))
    sqrt_e = sqrt(e)
    stable = N2 > 0
    length_ratio = ifelse(stable,
                          min(one(RealT),
                              convert(RealT, 0.76) * sqrt_e / (N * grid_spacing)),
                          one(RealT))
    e_to_3_2_over_length = ifelse(stable,
                                  max(e * sqrt_e / grid_spacing,
                                      e * N / convert(RealT, 0.76)),
                                  e * sqrt_e / grid_spacing)
    (; c_k) = turbulence
    c_epsilon = convert(RealT, 1.9) * c_k +
                (convert(RealT, 0.93) - convert(RealT, 1.9) * c_k) * length_ratio
    return c_epsilon * e_to_3_2_over_length
end

@inline function eddy_diffusivities(w, gradients, k, aux, turbulence::DeardorffTKE,
                                    equations)
    e = max(last(w), zero(eltype(w)))
    nu = turbulence.c_k * aux[3] * sqrt(e)
    return nu, nu
end

# Right-hand sides of the source terms of the conservative variables, i.e., the derivative
# of `integrate_turbulence_sources` at `dt = 0`
@inline function turbulence_sources(w, gradients, k, aux, turbulence::DeardorffTKE,
                                    equations)
    rho = w[1]
    e = max(last(w), zero(eltype(w)))
    nu, kappa = eddy_diffusivities(w, gradients, k, aux, turbulence, equations)
    shear_production = 2 * nu * strain_rate_vertical_squared(gradients, k)
    N2 = brunt_vaisala_frequency_squared(w, gradients, k, aux)
    buoyancy_production = -kappa * N2
    dissipation = dissipation_tke(e, N2, aux[3], turbulence)

    # Heating by the dissipation at constant pressure
    p = pressure_specific(w, equations)
    exner = (p / equations.p_0)^(equations.R / equations.c_p)
    heating = rho * dissipation / (equations.c_p * exner)

    zero_ = zero(rho)
    return SVector(zero_, zero_, zero_, zero_, heating,
                   rho * (shear_production + buoyancy_production - dissipation))
end

# Increments of the conservative variables by the exact solution of the source terms over
# the time `dt` with frozen gradients, density, and pressure: the TKE of the solution `s^2`
# of `ds/dt = (A - a s - b s^2) / 2` for `s = sqrt(e)`, and the heating by the dissipation
# `int epsilon dt = A int s dt - (e(dt) - e(0))`, see `DeardorffTKE`
@inline function integrate_turbulence_sources(w, gradients, k, aux, dt,
                                              turbulence::DeardorffTKE, equations)
    RealT = eltype(w)
    rho = w[1]
    e = max(last(w), zero(RealT))
    grid_spacing = aux[3]
    (; c_k) = turbulence
    N2 = brunt_vaisala_frequency_squared(w, gradients, k, aux)
    A = c_k * grid_spacing * (2 * strain_rate_vertical_squared(gradients, k) - N2)

    # Coefficients of `l = Delta` and, for `N^2 > 0`, of `l = 0.76 s / N < Delta`, which
    # holds for `s < N Delta / 0.76`
    b_large = convert(RealT, 0.93) / grid_spacing
    N = sqrt(max(N2, zero(N2)))
    a_small = convert(RealT, 1.9) * c_k * N / convert(RealT, 0.76)
    b_small = (convert(RealT, 0.93) - convert(RealT, 1.9) * c_k) / grid_spacing
    s_switch = ifelse(N2 > 0, N * grid_spacing / convert(RealT, 0.76),
                      zero(RealT))

    s, integral_s = integrate_two_regimes(sqrt(e), dt, A, a_small, b_small,
                                          zero(RealT), b_large, s_switch)
    e_new = s^2
    integral_dissipation = max(A * integral_s - (e_new - e), zero(RealT))

    p = pressure_specific(w, equations)
    exner = (p / equations.p_0)^(equations.R / equations.c_p)
    zero_ = zero(rho)
    return SVector(zero_, zero_, zero_, zero_,
                   rho * integral_dissipation / (equations.c_p * exner),
                   rho * (e_new - e))
end

# Solution `s(dt)` and `int_0^dt s dt` of `ds/dt = (A - a s - b s^2) / 2`, where the
# coefficients are `(a_small, b_small)` for `s < s_switch` and `(a_large, b_large)` otherwise.
# The solution is monotone, so it switches the regime at most once.
@inline function integrate_two_regimes(s0, dt, A, a_small, b_small, a_large, b_large,
                                       s_switch)
    if s0 < s_switch
        s, integral = integrate_riccati(s0, dt, A, a_small, b_small)
        if s > s_switch
            t_switch = time_riccati(s0, s_switch, A, a_small, b_small)
            _, integral = state_riccati(s0, t_switch, A, a_small, b_small)
            s, integral_large = integrate_riccati(s_switch, dt - t_switch, A, a_large,
                                                  b_large)
            integral += integral_large
        end
    else
        s, integral = integrate_riccati(s0, dt, A, a_large, b_large)
        if s < s_switch
            t_switch = time_riccati(s0, s_switch, A, a_large, b_large)
            _, integral = state_riccati(s0, t_switch, A, a_large, b_large)
            s, integral_small = integrate_riccati(s_switch, dt - t_switch, A, a_small,
                                                  b_small)
            integral += integral_small
        end
    end
    return s, integral
end

# Solution and its time integral of `ds/dt = (A - a s - b s^2) / 2` with `a >= 0`, `b > 0`,
# and `s(0) = s0 >= 0`, where `s = 0` is absorbing for `A <= 0`
@inline function integrate_riccati(s0, dt, A, a, b)
    if A < 0
        # Both roots of `A - a s - b s^2` are negative, so `s` reaches zero
        t_zero = time_riccati(s0, zero(s0), A, a, b)
        if dt >= t_zero
            _, integral = state_riccati(s0, t_zero, A, a, b)
            return zero(s0), integral
        end
    end
    return state_riccati(s0, dt, A, a, b)
end

# Closed-form solution `s(t)` and `int_0^t s dt` with the discriminant `D = a^2 + 4 b A`
@inline function state_riccati(s0, t, A, a, b)
    D = a^2 + 4 * b * A
    if D > 0
        # `s(t) = s_minus + (s_plus - s_minus) / (1 - K exp(-lambda t))` with the roots
        # `s_minus < 0` and `s_plus`
        sqrt_D = sqrt(D)
        s_plus = (sqrt_D - a) / (2 * b)
        s_minus = -(sqrt_D + a) / (2 * b)
        lambda = sqrt_D / 2
        K = (s0 - s_plus) / (s0 - s_minus)
        E = exp(-lambda * t)
        s = s_minus + (s_plus - s_minus) / (1 - K * E)
        integral = s_minus * t +
                   (s_plus - s_minus) * (t + log((1 - K * E) / (1 - K)) / lambda)
    elseif D < 0
        # `s(t) + c = omega tan(phi0 - b omega t / 2)`
        c = a / (2 * b)
        omega = sqrt(-D) / (2 * b)
        phi0 = atan((s0 + c) / omega)
        phi = phi0 - b * omega * t / 2
        s = omega * tan(phi) - c
        integral = -c * t + 2 * log(cos(phi) / cos(phi0)) / b
    else
        # Double root `s_root = -a / (2 b)`
        s_root = -a / (2 * b)
        s = s_root + 1 / (1 / (s0 - s_root) + b * t / 2)
        integral = s_root * t + 2 * log(1 + b * (s0 - s_root) * t / 2) / b
    end
    return s, integral
end

# Time at which the solution of `ds/dt = (A - a s - b s^2) / 2` from `s0` reaches `s_target`,
# which must lie between `s0` and the limit of the solution
@inline function time_riccati(s0, s_target, A, a, b)
    D = a^2 + 4 * b * A
    if D > 0
        sqrt_D = sqrt(D)
        s_plus = (sqrt_D - a) / (2 * b)
        s_minus = -(sqrt_D + a) / (2 * b)
        K = (s0 - s_plus) / (s0 - s_minus)
        K_target = (s_target - s_plus) / (s_target - s_minus)
        return 2 * log(K / K_target) / sqrt_D
    elseif D < 0
        c = a / (2 * b)
        omega = sqrt(-D) / (2 * b)
        return 2 * (atan((s0 + c) / omega) - atan((s_target + c) / omega)) / (b * omega)
    else
        s_root = -a / (2 * b)
        return 2 * (1 / (s_target - s_root) - 1 / (s0 - s_root)) / b
    end
end
end # @muladd
