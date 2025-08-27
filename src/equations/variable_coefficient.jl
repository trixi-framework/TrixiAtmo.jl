@muladd begin

struct VariableCoefficientAdvectionEquation2D{} <:
       AbstractVariableCoefficientEquations{2, 1} end 

varnames(::typeof(cons2cons), ::VariableCoefficientAdvectionEquation2D) = ("scalar",)
varnames(::typeof(cons2prim), ::VariableCoefficientAdvectionEquation2D) = ("scalar",)

end

@inline function Trixi.flux(u, aux_vars, orientation::Integer, equations::VariableCoefficientAdvectionEquation2D)
       a = aux_vars[orientation]
       return a * u
end

@inline function Trixi.flux(u, aux_vars, normal_direction::AbstractVector,
       equation::VariableCoefficientAdvectionEquation2D)
       a = dot(aux_vars, normal_direction) # velocity in normal direction
return a * u
end


# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation
#@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, aux_vars, orientation::Integer,
#       equation::VariableCoefficientAdvectionEquation2D)
#       λ_max = abs(aux_vars[orientation])
#       return λ_max 
#end

# Maximum wave speeds in each direction for CFL calculation
#@inline function Trixi.max_abs_speeds(u, aux_vars,
#       equations::VariableCoefficientAdvectionEquation2D)
#return abs.(aux_vars)
#end

@inline Trixi.cons2entropy(u, equations::VariableCoefficientAdvectionEquation2D) = u 
@inline Trixi.cons2prim(u, equations::VariableCoefficientAdvectionEquation2D) = u 
