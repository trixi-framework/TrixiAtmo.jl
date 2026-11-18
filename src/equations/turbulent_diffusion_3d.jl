@muladd begin
#! format: noindent

@doc raw"""
    TurbulentDiffusion3D(equations; lowest_level_height = 65)

Turbulent fluxes of the turbulence model `equations.turbulence`, e.g., a
[`HorizontalVerticalTurbulence`](@ref), of the
[`CompressibleEulerPotentialTemperatureEquationsWithGravity3D`](@ref) or the
[`CompressibleEulerMoistPotentialTemperatureEquations3D`](@ref), used as parabolic terms of a
`SemidiscretizationHyperbolicParabolic`,
```math
\begin{aligned}
\partial_t (\varrho \vec{v}) &= \nabla \cdot \tau, &
\tau &= 2 \varrho \left(\nu_h P S P + \nu_v (S - P S P)\right), \\
\partial_t (\varrho \psi) &= \nabla \cdot (\varrho D \nabla \psi), &
D &= \kappa_h P + \kappa_v \hat{k} \hat{k}^T,
\end{aligned}
```
with the rate-of-strain tensor ``S = (\nabla \vec{v} + \nabla \vec{v}^T) / 2``, the
projection ``P = I - \hat{k} \hat{k}^T`` onto the plane orthogonal to the upward unit
vector ``\hat{k}``, and the eddy viscosities ``\nu_h, \nu_v`` and diffusivities
``\kappa_h, \kappa_v`` of the horizontal and the vertical model. Here, ``\psi`` runs over
the specific variables of [`cons2specific`](@ref) except for the velocity, i.e., the
potential temperature ``\theta`` or the virtual potential temperature ``\theta_v`` and the
mass fractions ``\varrho_x / \varrho`` of the moisture variables, whose turbulent fluxes
sum to zero, such that the density has no turbulent flux. The stress ``\tau`` is
symmetric, it reduces to ``2 \varrho \nu S`` for ``\nu_h = \nu_v = \nu``, the vertical
shear of the horizontal velocity is diffused with ``\nu_v``, and its dissipation is
```math
\tau : S = 2 \varrho \left(\nu_h \lVert P S P \rVert^2
                           + \nu_v \lVert S - P S P \rVert^2\right) \geq 0.
```

The upward unit vector is ``\hat{k} = \nabla \phi / \lVert \nabla \phi \rVert`` of the
geopotential ``\phi`` given by the `aux_field` passed to
`SemidiscretizationHyperbolicParabolic`. For the
[`CompressibleEulerMoistPotentialTemperatureEquations3D`](@ref), it is the ``\hat{k}`` of
their auxiliary variables, i.e., the direction of the falling water. The pointwise
functions take ``\hat{k}`` and the variables `column` of the lowest level of the column
of the node as additional arguments after the gradients, e.g.,
`flux(w, gradients, k, column, orientation, equations)`. They are the speed of the horizontal
velocity `column[1]` and the height `column[2]` of the lowest level above the surface, used
by column models such as the [`ReedJablonowskiBoundaryLayer`](@ref). The lowest level of a
column is the point of the column at the height `lowest_level_height` in m above the surface,
which must lie in the lowest element, and its state is the value of the polynomial of the
conservative variables there.
"""
struct TurbulentDiffusion3D{E, NVARS, RealT <: Real} <:
       Trixi.AbstractEquationsParabolic{3, NVARS, GradientVariablesPrimitive}
    equations_hyperbolic::E
    lowest_level_height::RealT # height of the lowest level of the columns in m
end

function TurbulentDiffusion3D(equations::Union{CompressibleEulerPotentialTemperatureEquationsWithGravity3D,
                                               CompressibleEulerMoistPotentialTemperatureEquations3D};
                              lowest_level_height = 65)
    if !(equations.turbulence isa HorizontalVerticalTurbulence)
        throw(ArgumentError("`TurbulentDiffusion3D` requires equations with a " *
                            "`HorizontalVerticalTurbulence`, but the turbulence of the " *
                            "equations is `$(nameof(typeof(equations.turbulence)))`"))
    end
    lowest_level_height = convert(typeof(equations.g), lowest_level_height)
    return TurbulentDiffusion3D{typeof(equations), nvariables(equations),
                                typeof(lowest_level_height)}(equations,
                                                             lowest_level_height)
end

function Base.similar(equations_parabolic::TurbulentDiffusion3D,
                      ::Type{NewRealT}) where {NewRealT}
    return TurbulentDiffusion3D(similar(equations_parabolic.equations_hyperbolic,
                                        NewRealT);
                                lowest_level_height = equations_parabolic.lowest_level_height)
end

function varnames(variable_mapping, equations_parabolic::TurbulentDiffusion3D)
    return varnames(variable_mapping, equations_parabolic.equations_hyperbolic)
end

Trixi.gradient_variable_transformation(::TurbulentDiffusion3D) = cons2specific

@doc raw"""
    cons2specific(u, equations)

Specific variables, i.e., the density ``\varrho`` followed by the other conservative
variables divided by ``\varrho``: the velocity, the potential temperature ``\theta`` or the
virtual potential temperature ``\theta_v``, and the mass fractions ``\varrho_x / \varrho`` of
the moisture variables. These are the gradient variables of
[`TurbulentDiffusion3D`](@ref).
"""
@inline function cons2specific(u,
                               equations::Union{CompressibleEulerPotentialTemperatureEquationsWithGravity3D,
                                                CompressibleEulerMoistPotentialTemperatureEquations3D})
    rho = u[1]
    return SVector(rho, ntuple(i -> u[i + 1] / rho, Val(nvariables(equations) - 1))...)
end

@inline function cons2specific(u, equations_parabolic::TurbulentDiffusion3D)
    return cons2specific(u, equations_parabolic.equations_hyperbolic)
end

# Conservative variables of the specific variables `w`, i.e., the inverse of `cons2specific`
@inline function specific2cons(w)
    rho = w[1]
    return SVector(rho, ntuple(i -> rho * w[i + 1], Val(length(w) - 1))...)
end

# Pressure `p = K (rho theta)^gamma` of the specific variables `w = cons2specific(u, equations)`,
# with the virtual potential temperature for the moist equations
@inline function pressure_specific(w,
                                   equations::Union{CompressibleEulerPotentialTemperatureEquationsWithGravity3D,
                                                    CompressibleEulerMoistPotentialTemperatureEquations3D})
    return equations.K * exp(equations.gamma * log(w[1] * w[5]))
end

# Upward unit vector `k = ∇φ / |∇φ|` of the geopotential `φ = aux_field(x, equations)[1]`
@inline function upward_unit_vector(x, aux_field, equations)
    gradient_phi = ForwardDiff.gradient(y -> first(aux_field(y, equations)), x)
    return gradient_phi / norm(gradient_phi)
end

# The auxiliary variables `aux[2:4]` are the upward unit vector, see `vertical_component`
@inline function upward_unit_vector(x, aux_field,
                                    equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    aux = aux_field(x, equations)
    return SVector(aux[2], aux[3], aux[4])
end

# Turbulent flux in the direction `orientation` of the specific variables `w` with the
# gradients `(dw/dx, dw/dy, dw/dz)`, the upward unit vector `k`, and the variables `column` of
# the lowest level of the column, see `eddy_diffusivities`
@inline function flux(w, gradients, k, column, orientation::Integer,
                      equations_parabolic::TurbulentDiffusion3D)
    gradients_x, gradients_y, gradients_z = gradients
    rho = w[1]
    (; equations_hyperbolic) = equations_parabolic
    nu_h, kappa_h, nu_v, kappa_v = eddy_diffusivities(w, gradients, k, column,
                                                      equations_hyperbolic.turbulence,
                                                      equations_hyperbolic)

    # Rate-of-strain tensor `S` of the velocity gradient `dv_i/dx_j = gradients_j[1 + i]`
    velocity_gradient = SMatrix{3, 3}(gradients_x[2], gradients_x[3], gradients_x[4],
                                      gradients_y[2], gradients_y[3], gradients_y[4],
                                      gradients_z[2], gradients_z[3], gradients_z[4])
    S = 0.5f0 * (velocity_gradient + velocity_gradient')
    S_k = S * k
    k_S_k = dot(k, S_k)

    # Row `orientation` of the symmetric stress
    # `tau = 2 rho (nu_h S + (nu_v - nu_h) (S - P S P))` with
    # `S - P S P = k (S k)' + (S k) k' - (k' S k) k k'`
    k_d = k[orientation]
    stress = 2 * rho *
             (nu_h * S[:, orientation] +
              (nu_v - nu_h) * (k_d * S_k + S_k[orientation] * k - k_S_k * k_d * k))

    # Fluxes `rho (kappa_h g_d + (kappa_v - kappa_h) (k . g) k_d)` of the other specific
    # variables with the gradient `g`
    nscalars = Val(nvariables(equations_parabolic) - 4)
    gradients_d = gradients[orientation]
    scalars = ntuple(nscalars) do i
        v = i + 4
        k_dot_g = k[1] * gradients_x[v] + k[2] * gradients_y[v] + k[3] * gradients_z[v]
        return rho * (kappa_h * gradients_d[v] + (kappa_v - kappa_h) * k_dot_g * k_d)
    end

    return SVector(zero(rho), stress[1], stress[2], stress[3], scalars...)
end

# No turbulent flux through the boundary for a vanishing `boundary_normal_flux_function`,
# as for the `LaplaceDiffusion3D` of Trixi.jl. The `flux_inner` of the `Gradient` is the
# boundary value of the specific variables.
@inline function (boundary_condition::BoundaryConditionNeumann)(flux_inner, u_inner,
                                                                normal::AbstractVector,
                                                                x, t,
                                                                operator_type::Trixi.Gradient,
                                                                equations_parabolic::TurbulentDiffusion3D)
    return flux_inner
end

@inline function (boundary_condition::BoundaryConditionNeumann)(flux_inner, u_inner,
                                                                normal::AbstractVector,
                                                                x, t,
                                                                operator_type::Trixi.Divergence,
                                                                equations_parabolic::TurbulentDiffusion3D)
    # `flux_inner` is the turbulent flux dotted with `normal`, which need not be a unit
    # vector, so the prescribed normal flux is scaled by `norm(normal)` as well
    return (boundary_condition.boundary_normal_flux_function(x, t,
                                                             equations_parabolic) *
            norm(normal))
end

@doc raw"""
    BoundaryConditionSurfaceFluxes(sea_surface_temperature)

Bulk fluxes through the surface of the momentum and the potential temperature of the
[`ReedJablonowskiBoundaryLayer`](@ref), used as parabolic boundary condition of
[`TurbulentDiffusion3D`](@ref) for the
[`CompressibleEulerPotentialTemperatureEquationsWithGravity3D`](@ref) with a
`HorizontalVerticalTurbulence` whose vertical model is a `ReedJablonowskiBoundaryLayer`.
With the density ``\varrho_a``, the horizontal velocity ``\vec{v}_a``, the potential
temperature ``\theta_a``, and the pressure ``p_a`` of the lowest level of the column (see
[`TurbulentDiffusion3D`](@ref)), the turbulent fluxes into the domain are
```math
-\varrho_a C_d \lVert \vec{v}_a \rVert \vec{v}_a, \qquad
\varrho_a C_E \lVert \vec{v}_a \rVert
\left(T_s \left(\frac{p_0}{p_a}\right)^{R / c_p} - \theta_a\right)
```
of the momentum and of ``\varrho \theta`` (Ullrich et al., 2016, Appendix C), with the
drag coefficient ``C_d`` of [`drag_coefficient`](@ref), the exchange coefficient ``C_E``,
and the sea surface temperature ``T_s`` `= sea_surface_temperature(x, equations)` at the
surface node `x`. The density has no flux through the surface.
"""
struct BoundaryConditionSurfaceFluxes{SeaSurfaceTemperature}
    sea_surface_temperature::SeaSurfaceTemperature
end

# The `flux_inner` of the `Gradient` is the boundary value of the specific variables
@inline function (boundary_condition::BoundaryConditionSurfaceFluxes)(flux_inner,
                                                                      u_inner,
                                                                      normal::AbstractVector,
                                                                      x, t,
                                                                      operator_type::Trixi.Gradient,
                                                                      equations_parabolic::TurbulentDiffusion3D)
    return flux_inner
end

# Turbulent flux through the surface in the direction of the outward `normal`, which need not
# be a unit vector, from the specific variables `w_a`, the upward unit vector `k_a`, and the
# height above the surface `z_a` of the lowest level of the column
@inline function (boundary_condition::BoundaryConditionSurfaceFluxes)(w_a, k_a, z_a,
                                                                      normal::AbstractVector,
                                                                      x,
                                                                      equations_parabolic::TurbulentDiffusion3D{<:CompressibleEulerPotentialTemperatureEquationsWithGravity3D{<:Any,
                                                                                                                                                                              <:HorizontalVerticalTurbulence{<:Any,
                                                                                                                                                                                                             <:ReedJablonowskiBoundaryLayer}}})
    (; equations_hyperbolic) = equations_parabolic
    turbulence = equations_hyperbolic.turbulence.vertical
    rho_a = w_a[1]
    v_a = SVector(w_a[2], w_a[3], w_a[4])
    v_horizontal = v_a - dot(v_a, k_a) * k_a
    wind_speed = norm(v_horizontal)
    p_a = pressure_specific(w_a, equations_hyperbolic)
    T_s = boundary_condition.sea_surface_temperature(x, equations_hyperbolic)
    theta_s = T_s *
              (equations_hyperbolic.p_0 / p_a)^(equations_hyperbolic.R /
                                                equations_hyperbolic.c_p)

    norm_ = norm(normal)
    momentum = -norm_ * rho_a * drag_coefficient(wind_speed, turbulence) * wind_speed *
               v_horizontal
    heat = norm_ * rho_a * turbulence.exchange_coefficient * wind_speed *
           (theta_s - w_a[5])
    return SVector(zero(heat), momentum[1], momentum[2], momentum[3], heat)
end
end # @muladd
