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


# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation
@inline function max_abs_speed_naive(u_ll, u_rr, aux_vars, orientation::Integer,
       equation::LinearScalarAdvectionEquation2D)
λ_max = abs(aux_vars[orientation])
end

end