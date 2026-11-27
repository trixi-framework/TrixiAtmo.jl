# Changelog

TrixiAtmo.jl follows the interpretation of
[semantic versioning (semver)](https://julialang.github.io/Pkg.jl/dev/compatibility/#Version-specifier-format-1)
used in the Julia ecosystem. Notable changes will be documented in this file
for human readability.

## Changes when updating to v0.2 from v0.1.x

#### Changed
- The `CompressibleEulerPotentialTemperatureEquationsWithGravity3D`, the
  `CompressibleEulerEnergyEquationsWithGravity3D`, and the
  `CompressibleEulerVectorInvariantEquationsWithGravity3D` store the geopotential as
  auxiliary variable of Trixi.jl (see `Trixi.nauxvariables`) instead of as last variable
  of the state. Thus, they have one variable less, their pointwise functions take the
  auxiliary variables as additional arguments after the states, e.g.,
  `flux_ec(u_ll, u_rr, aux_ll, aux_rr, normal_direction, equations)`,
  `source_terms_coriolis(u, aux, x, t, equations)`, and `cons2prim(u, aux, equations)`
  (also `cons2entropy` and `entropy`, as required by Trixi.jl v0.18), and
  `SemidiscretizationHyperbolic` requires the keyword argument `aux_field`, e.g.,
  `aux_field_gravity_earth`.
  With adaptive time stepping, results change slightly since the geopotential no longer
  enters the error estimate of the time integrator.
- `energy_total` of the `CompressibleEulerEnergyEquationsWithGravity{2,3}D` and the
  `CompressibleEulerPotentialTemperatureEquationsWithGravity{2,3}D` is renamed to
  `energy_internal_kinetic`, since it is the sum of the internal and the kinetic energy
  without the potential energy of gravity.

## Changes in the v0.1 lifecycle

#### Added
- The `CompressibleEulerMoistPotentialTemperatureEquations3D` are the compressible Euler
  equations with gravity with the virtual potential temperature as prognostic variable and
  the moisture variables of a microphysics scheme, whose right-hand sides are given by
  `source_terms_microphysics`. The `KesslerMicrophysics` is the warm-rain Kessler
  microphysics of DCMIP2016 with water vapor, cloud water, and falling rain water.
  The equations store the geopotential and the upward unit vector,
  which gives the direction of the falling water,
  as auxiliary variables of Trixi.jl, e.g., given by `aux_field_gravity_earth`.
  Their volume flux `flux_kennedy_gruber_souza_etal` supports `FluxTurbo` with the
  `KesslerMicrophysics`.
  The moist baroclinic wave of DCMIP2016 is available as
  `initial_condition_moist_baroclinic_instability` ([#271]).
- `JacobianVerticalTwoSided` enables horizontally explicit, vertically implicit (HEVI) time
  integration of the `CompressibleEulerPotentialTemperatureEquationsWithGravity3D` on a
  `P4estMeshCubedSphereTopography` as linear operator of linearly implicit methods, i.e.,
  Rosenbrock-W methods such as `SSPRW32` and linear IMEX methods such as `ARS343`, also on
  GPUs ([#243]).
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
