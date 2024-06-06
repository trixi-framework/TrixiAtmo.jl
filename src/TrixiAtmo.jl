module TrixiAtmo

using Trixi
# Import additional symbols that are not exported by Trixi.jl
# using Trixi:
using MuladdMacro: @muladd
using StaticArrays: SVector
using Static: True, False
using LinearAlgebra: norm

include("equations/equations.jl")

end # module TrixiAtmo
