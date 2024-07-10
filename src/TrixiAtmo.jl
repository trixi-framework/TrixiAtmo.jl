module TrixiAtmo

using Trixi
# Import additional symbols that are not exported by Trixi.jl
# using Trixi:
using MuladdMacro: @muladd
using StaticArrays: SVector
using Static: True, False
using LinearAlgebra: norm

foo() = true
bar() = false
baz() = Trixi.examples_dir()

include("equations/equations.jl")

end # module TrixiAtmo
