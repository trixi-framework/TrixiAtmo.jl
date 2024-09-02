using Trixi: AbstractEquations

abstract type AbstractCompressibleMoistEulerEquations{NDIMS, NVARS} <:
              AbstractEquations{NDIMS, NVARS} end

abstract type AbstractCompressibleRainyEulerEquations{NDIMS, NVARS} <:
              AbstractEquations{NDIMS, NVARS} end

include("compressible_rainy_euler_2d.jl")
include("compressible_moist_euler_2d_lucas.jl")