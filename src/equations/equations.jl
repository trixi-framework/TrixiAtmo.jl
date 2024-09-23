using Trixi: AbstractEquations

abstract type AbstractCompressibleMoistEulerEquations{NDIMS, NVARS} <:
              AbstractEquations{NDIMS, NVARS} end
abstract type AbstractCompressibleEulerPotentialTemperatureEquations{NDIMS, NVARS} <:
              AbstractEquations{NDIMS, NVARS} end


              struct flux_LMARS{SpeedOfSound}
                # Estimate for the speed of sound
                speed_of_sound::SpeedOfSound
              end
                              
include("compressible_moist_euler_2d_lucas.jl")
include("compressible_euler_potential_temperature_2d.jl")
