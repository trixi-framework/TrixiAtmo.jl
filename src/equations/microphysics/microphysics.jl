# Microphysics schemes of the moist equations, e.g.,
# `CompressibleEulerMoistPotentialTemperatureEquations3D`.
# A scheme only sees the thermodynamic state, i.e., the densities of the moisture variables,
# the density of dry air, the temperature, and the pressure, such that it does not depend on
# the spatial dimension or on the prognostic thermodynamic variable of the equations.
# The first moisture variable is the water vapor. A scheme defines
# - `nmoisturevariables(microphysics)`: the number of moisture variables,
# - `varnames_moisture(variables, microphysics)`: their names for `cons2cons` and
#   `cons2prim`,
# - `density_water(moisture, microphysics)`: the sum of the moisture variables that are
#   densities of water, which are part of the total density,
# - `sedimentation(moisture, rho_dry, microphysics)`: the downward fluxes of the falling
#   water through a horizontal surface, i.e., of the total density and of each moisture
#   variable,
# - `max_terminal_velocity(moisture, rho_dry, microphysics)`: the largest terminal velocity
#   of the falling water,
# - `microphysics_tendencies(moisture, rho_dry, T, p, microphysics)`: the isobaric heating
#   rate in K/s and the right-hand sides of the moisture variables.

include("kessler.jl")
