# # Rising bubble (dry Euler + gravity)
#
# This tutorial demonstrates how to set up a dry rising-bubble test case in
# `TrixiAtmo.jl`.
#
# It is based on the existing buoyancy example distributed with the package.
# The goal here is to explain the physical setup, numerical ingredients, and
# expected solution behavior step by step.
#
# ## Load packages

using Trixi
using TrixiAtmo

# ## Overview
#
# The rising-bubble problem starts from a background atmospheric state with a
# localized warm perturbation. Because the perturbation is buoyant, it rises as
# the simulation evolves.
#
# In the following sections, we adapt the existing example into a more explicit
# tutorial.

# TODO: copy and explain the example blocks here