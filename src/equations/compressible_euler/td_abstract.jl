# Abstract type for thermodynamic parts of CompressibleEulerAtmo equations
#
# Generic functions go here
# Files including concrete types are included at the end

abstract type AbstractThermodynamicEquation end

include("td_total_energy.jl")
include("td_potential_temperature.jl")
