# Abstract type for thermodynamic parts of CompressibleEulerAtmo equations
#
# Generic functions go here
# Files including concrete types are included at the end

abstract type AbstractThermodynamicEquation end

include("td_energy_total.jl")
include("td_energy_internal.jl")
#include("td_energy_total_potential.jl")
include("td_potential_temperature.jl")
