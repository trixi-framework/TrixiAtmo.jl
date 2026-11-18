# Turbulence models of the equations with gravity on the sphere, e.g.,
# `CompressibleEulerPotentialTemperatureEquationsWithGravity3D` and
# `CompressibleEulerMoistPotentialTemperatureEquations3D`, whose turbulent fluxes are the
# parabolic terms of `TurbulentDiffusion3D`. A model only sees the specific variables `w` of
# `cons2specific`, their gradients, the upward unit vector `k`, and the variables `column`
# of the lowest level of the column of the node, see `TurbulentDiffusion3D`: the speed of the
# horizontal velocity `column[1]` and the height above the surface `column[2]`. A model
# defines
# - `nturbulencevariables(turbulence)`: the number of prognostic variables of the model,
# - `eddy_diffusivities(w, gradients, k, column, turbulence, equations)`: the eddy viscosity
#   of the momentum and the eddy diffusivity of the other specific variables in m²/s, where
#   `equations` are the hyperbolic equations, e.g., to compute the pressure of `w`.
# A `HorizontalVerticalTurbulence` combines a horizontal and a vertical model, whose
# diffusivities act in the plane orthogonal to `k` and in the direction of `k`, respectively,
# see `flux(w, gradients, k, orientation, equations_parabolic::TurbulentDiffusion3D)`.

"""
    NoTurbulence()

No turbulence model, the default of the equations, which then have no parabolic terms.
"""
struct NoTurbulence end

Base.similar(turbulence::NoTurbulence, ::Type) = turbulence

@inline nturbulencevariables(::NoTurbulence) = 0

"""
    HorizontalVerticalTurbulence(horizontal, vertical)

Turbulence model whose eddy viscosity and diffusivity are given by the model `horizontal`
in the plane orthogonal to the upward unit vector ``\\hat{k}`` and by the model `vertical`
in the direction of ``\\hat{k}``, e.g., `ConstantEddyViscosity`s. See
[`TurbulentDiffusion3D`](@ref) for the resulting turbulent fluxes.
"""
struct HorizontalVerticalTurbulence{Horizontal, Vertical}
    horizontal::Horizontal
    vertical::Vertical
end

function Base.similar(turbulence::HorizontalVerticalTurbulence,
                      ::Type{NewRealT}) where {NewRealT}
    return HorizontalVerticalTurbulence(similar(turbulence.horizontal, NewRealT),
                                        similar(turbulence.vertical, NewRealT))
end

@inline function nturbulencevariables(turbulence::HorizontalVerticalTurbulence)
    return (nturbulencevariables(turbulence.horizontal) +
            nturbulencevariables(turbulence.vertical))
end

# Eddy viscosities and diffusivities `(nu_h, kappa_h, nu_v, kappa_v)` of the horizontal and
# the vertical model
@inline function eddy_diffusivities(w, gradients, k, column,
                                    turbulence::HorizontalVerticalTurbulence, equations)
    nu_h, kappa_h = eddy_diffusivities(w, gradients, k, column, turbulence.horizontal,
                                       equations)
    nu_v, kappa_v = eddy_diffusivities(w, gradients, k, column, turbulence.vertical,
                                       equations)
    return nu_h, kappa_h, nu_v, kappa_v
end

include("constant_eddy_viscosity.jl")
include("reed_jablonowski.jl")
