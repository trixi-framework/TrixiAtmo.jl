# Changelog

TrixiAtmo.jl follows the interpretation of
[semantic versioning (semver)](https://julialang.github.io/Pkg.jl/dev/compatibility/#Version-specifier-format-1)
used in the Julia ecosystem. Notable changes will be documented in this file
for human readability.

## Changes in the v0.1 lifecycle

#### Added
- The `CompressibleEulerVectorInvariantEquationsWithGravity3D` are the compressible Euler
  equations with gravity in vector-invariant form with the potential temperature as
  prognostic variable. They come with the entropy and total energy conservative volume flux
  `flux_etec`, which supports `FluxTurbo`, and the entropy-stable and well-balanced surface
  flux `flux_lmars_es` ([#253]).
- The new nonconservative flux `flux_nonconservative_isothermal` is a well-balanced gravity
  term for isothermal background states of the
  `CompressibleEulerEnergyEquationsWithGravity2D` and the
  `CompressibleEulerEnergyEquationsWithGravity3D` ([#248]).
- GPU support including adaptation to `Float32` via `trixi_adapt` ([#242]).

#### Changed
- The `CompressibleEulerPotentialTemperatureEquationsWithGravity2D`, the
  `CompressibleEulerEnergyEquationsWithGravity2D`, and the
  `CompressibleEulerInternalEnergyEquationsWithGravity2D` store the geopotential as
  auxiliary variable of Trixi.jl (see `Trixi.nauxvariables`) instead of as last variable
  of the state. Thus, they have one variable less, their pointwise functions take the
  auxiliary variables as additional arguments after the states, e.g.,
  `flux_ec(u_ll, u_rr, aux_ll, aux_rr, normal_direction, equations)`, and
  `SemidiscretizationHyperbolic` requires the keyword argument `aux_field`, e.g., the new
  `aux_field_geopotential`. Since auxiliary variables are only available for `P4estMesh`es,
  the elixirs `elixir_potential_temperature_inertia_gravity_waves.jl` and
  `elixir_energy_inertia_gravity_waves.jl` use a `P4estMesh` instead of a `StructuredMesh`.
- Fixes a bug in `flux_nonconservative_waruszewski_etal` for the
  `CompressibleEulerEnergyEquationsWithGravity2D` and the
  `CompressibleEulerEnergyEquationsWithGravity3D`: The gravity term of the momentum did
  not use the density of Eq. (73d) of Waruszewski et al. (2022). Thus, results of
  simulations with `flux_nonconservative_waruszewski_etal` change except for isothermal
  states. The previous gravity term is available as `flux_nonconservative_isothermal`
  ([#248]).
- The gravity term `flux_nonconservative_waruszewski_etal` of the
  `CompressibleEulerPotentialTemperatureEquationsWithGravity1D`, the
  `CompressibleEulerPotentialTemperatureEquationsWithGravity2D`, and the
  `CompressibleEulerPotentialTemperatureEquationsWithGravity3D` is renamed to
  `flux_nonconservative_isothermal`, since it does not correspond to the one of
  Waruszewski et al. (2022) ([#248]).
