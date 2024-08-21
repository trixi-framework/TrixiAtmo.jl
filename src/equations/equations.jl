using Trixi: AbstractEquations

# Physical constants
const EARTH_RADIUS = 6.37122f6  # m
const EARTH_GRAVITATIONAL_ACCELERATION = 9.80616f0  # m/s²
const EARTH_ROTATION_RATE = 7.292f-5  # rad/s
const SECONDS_PER_DAY = 8.64f4

abstract type AbstractCovariantEquations2D{NVARS} <: AbstractEquations{2, NVARS} end
include("covariant_advection.jl")

abstract type AbstractCompressibleMoistEulerEquations{NDIMS, NVARS} <:
              AbstractEquations{NDIMS, NVARS} end
include("compressible_moist_euler_2d_lucas.jl")
