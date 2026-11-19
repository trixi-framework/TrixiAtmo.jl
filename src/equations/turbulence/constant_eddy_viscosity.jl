@muladd begin
#! format: noindent

@doc raw"""
    ConstantEddyViscosity(viscosity; prandtl_number = 1)

Constant eddy viscosity ``\nu`` `= viscosity` of the momentum in m²/s and the eddy
diffusivity ``\kappa = \nu / \mathrm{Pr}_t`` of the other specific variables with the
turbulent Prandtl number ``\mathrm{Pr}_t`` `= prandtl_number`, used as the horizontal or the
vertical model of a [`HorizontalVerticalTurbulence`](@ref).
"""
struct ConstantEddyViscosity{RealT <: Real}
    viscosity::RealT # eddy viscosity of the momentum in m²/s
    diffusivity::RealT # eddy diffusivity of the other specific variables in m²/s
end

function ConstantEddyViscosity(viscosity; prandtl_number = 1)
    viscosity, diffusivity = promote(viscosity, viscosity / prandtl_number)
    return ConstantEddyViscosity(viscosity, diffusivity)
end

function Base.similar(turbulence::ConstantEddyViscosity,
                      ::Type{NewRealT}) where {NewRealT}
    return ConstantEddyViscosity(convert(NewRealT, turbulence.viscosity),
                                 convert(NewRealT, turbulence.diffusivity))
end

@inline nturbulencevariables(::ConstantEddyViscosity) = 0

@inline function eddy_diffusivities(w, gradients, k, aux,
                                    turbulence::ConstantEddyViscosity, equations)
    return turbulence.viscosity, turbulence.diffusivity
end
end # @muladd
