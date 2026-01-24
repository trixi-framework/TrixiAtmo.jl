@muladd begin

struct VariableCoefficientAdvectionEquation2D{} <:
       AbstractVariableCoefficientEquations{2, 1} end 

varnames(::typeof(cons2cons), ::VariableCoefficientAdvectionEquation2D) = ("scalar",)
varnames(::typeof(cons2prim), ::VariableCoefficientAdvectionEquation2D) = ("scalar",)
varnames(::typeof(cons2aux), ::VariableCoefficientAdvectionEquation2D) = ("v1", "v2")
varnames(::typeof(cons2prim_and_aux), ::VariableCoefficientAdvectionEquation2D) = ("scalar", "v1", "v2")

@inline function Trixi.flux(u, aux_vars, orientation::Integer, equations::VariableCoefficientAdvectionEquation2D)
       a = aux_vars[orientation]
       return a * u
end

@inline function Trixi.flux(u, aux_vars, normal_direction::AbstractVector,
       equation::VariableCoefficientAdvectionEquation2D)
       a = dot(aux_vars, normal_direction) # velocity in normal direction
return a * u
end

@inline function (numflux::FluxPlusDissipation)(u_ll, u_rr, aux_ll, aux_rr,
                                                orientation_or_normal_direction,
                                                equations)
    @unpack numerical_flux, dissipation = numflux

    return (numerical_flux(u_ll, u_rr, aux_ll, aux_rr,
                           orientation_or_normal_direction, equations)
            +
            dissipation(u_ll, u_rr, aux_ll, aux_rr,
                           orientation_or_normal_direction, equations))
end

const flux_lax_friedrichs = FluxLaxFriedrichs()

#function FluxLaxFriedrichs(max_abs_speed = max_abs_speed)
#    FluxPlusDissipation(flux_central, DissipationLocalLaxFriedrichs(max_abs_speed))
#end

# same as above for equations with auxiliary variables
@inline function (dissipation::DissipationLocalLaxFriedrichs)(u_ll, u_rr, aux_ll,
                                                              aux_rr,
                                                              orientation_or_normal_direction,
                                                              equations)
    λ = dissipation.max_abs_speed(u_ll, u_rr, aux_ll, aux_rr,
                                  orientation_or_normal_direction, equations)
    return -0.5f0 * λ * (u_rr - u_ll)
end

# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation
@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, aux_ll, aux_rr, normal_direction::AbstractVector,
       equation::VariableCoefficientAdvectionEquation2D)
       λ_ll = dot(aux_ll, normal_direction)
       λ_rr = dot(aux_rr, normal_direction)
       return max(abs(λ_ll), abs(λ_rr)) 
end

# Maximum wave speeds in each direction for CFL calculation
@inline function Trixi.max_abs_speeds(u, aux_vars, equations::VariableCoefficientAdvectionEquation2D)
       return abs.(aux_vars)
end



@inline Trixi.cons2entropy(u, equations::VariableCoefficientAdvectionEquation2D) = u 
@inline Trixi.cons2prim(u, equations::VariableCoefficientAdvectionEquation2D) = u 
#@inline Trixi.cons2prim(u, aux, equations::LinearVariableScalarAdvectionEquation2D) = u


end #muladd
