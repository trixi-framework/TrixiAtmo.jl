# Turbulence models of the equations with gravity on the sphere, e.g.,
# `CompressibleEulerPotentialTemperatureEquationsWithGravity3D` and
# `CompressibleEulerMoistPotentialTemperatureEquations3D`, whose turbulent fluxes are the
# parabolic terms of `TurbulentDiffusion3D`. A model only sees the specific variables `w` of
# `cons2specific`, their gradients, the upward unit vector `k`, and the auxiliary variables
# `aux` of the node: the speed of the horizontal velocity `aux[1]` and the height above the
# surface `aux[2]` of the lowest level of the column, see `TurbulentDiffusion3D`, the vertical
# grid spacing `aux[3]`, i.e., the height of the element divided by the number of nodes in the
# vertical, and the gravitational acceleration `aux[4] = |∇φ|`. A model defines
# - `nturbulencevariables(turbulence)`: the number of prognostic variables of the model,
# - `varnames_turbulence(variables, turbulence)`: their names for `cons2cons` and
#   `cons2prim`, where the primitive variables are the specific ones,
# - `eddy_diffusivities(w, gradients, k, aux, turbulence, equations)`: the eddy viscosity of
#   the momentum and the eddy diffusivity of the other specific variables in m²/s, where
#   `equations` are the hyperbolic equations, e.g., to compute the pressure of `w`,
# - `have_turbulence_sources(turbulence)`, `turbulence_sources(w, gradients, k, aux,
#   turbulence, equations)`, and `integrate_turbulence_sources(w, gradients, k, aux, dt,
#   turbulence, equations)`: whether the model has source terms, e.g., of its prognostic
#   variables, their right-hand sides, and the increments of the conservative variables by
#   their exact solution over the time `dt` with frozen gradients, which
#   `turbulence_sources_exact!` applies by splitting.
# A `HorizontalVerticalTurbulence` combines a horizontal and a vertical model, whose
# diffusivities act in the plane orthogonal to `k` and in the direction of `k`, respectively,
# see `flux(w, gradients, k, aux, orientation, equations_parabolic::TurbulentDiffusion3D)`.

varnames_turbulence(variables, turbulence) = ()

@inline have_turbulence_sources(turbulence) = false

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
@inline function eddy_diffusivities(w, gradients, k, aux,
                                    turbulence::HorizontalVerticalTurbulence, equations)
    nu_h, kappa_h = eddy_diffusivities(w, gradients, k, aux, turbulence.horizontal,
                                       equations)
    nu_v, kappa_v = eddy_diffusivities(w, gradients, k, aux, turbulence.vertical,
                                       equations)
    return nu_h, kappa_h, nu_v, kappa_v
end

function varnames_turbulence(variables, turbulence::HorizontalVerticalTurbulence)
    return (varnames_turbulence(variables, turbulence.horizontal)...,
            varnames_turbulence(variables, turbulence.vertical)...)
end

@inline function have_turbulence_sources(turbulence::HorizontalVerticalTurbulence)
    return (have_turbulence_sources(turbulence.horizontal) ||
            have_turbulence_sources(turbulence.vertical))
end

# Only the vertical model has source terms and prognostic variables, see `DeardorffTKE`
@inline function turbulence_sources(w, gradients, k, aux,
                                    turbulence::HorizontalVerticalTurbulence, equations)
    return turbulence_sources(w, gradients, k, aux, turbulence.vertical, equations)
end

@inline function integrate_turbulence_sources(w, gradients, k, aux, dt,
                                              turbulence::HorizontalVerticalTurbulence,
                                              equations)
    return integrate_turbulence_sources(w, gradients, k, aux, dt, turbulence.vertical,
                                        equations)
end

include("constant_eddy_viscosity.jl")
include("reed_jablonowski.jl")
include("deardorff_tke.jl")
