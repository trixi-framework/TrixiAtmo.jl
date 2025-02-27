@muladd begin

struct VariableCoefficientAdvectionEquation2D{} <:
       AbstractVariableCoefficientEquations{2, 1} end 

varnames(::typeof(cons2cons), ::VariableCoefficientAdvectionEquation2D) = ("scalar",)
varnames(::typeof(cons2prim), ::VariableCoefficientAdvectionEquation2D) = ("scalar",)

end

@inline function Trixi.flux(u, aux_vars, orientation::Integer, equations::VariableCoefficientAdvectionEquation2D)
       
end
