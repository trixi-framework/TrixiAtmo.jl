# Changelog

TrixiAtmo.jl follows the interpretation of
[semantic versioning (semver)](https://julialang.github.io/Pkg.jl/dev/compatibility/#Version-specifier-format-1)
used in the Julia ecosystem. Notable changes will be documented in this file
for human readability.

## Changes when updating to v0.2 from v0.1.x

#### Removed
- TrixiAtmo.jl now consistently uses `P4estMesh` for the Cartesian equations in 2D and 3D.
  Thus, the methods based on `orientation` (for `TreeMesh` and `StructuredMesh`) of
  `max_abs_speed_naive`, `max_abs_speed`, `flux_ec_rain`, `boundary_condition_slip_wall`,
  and `boundary_condition_simple_slip_wall` are removed for the
  `CompressibleEulerPotentialTemperatureEquations2D`, the
  `CompressibleEulerPotentialTemperatureEquations3D`, the
  `CompressibleEulerPotentialTemperatureEquationsWithGravity2D`, the
  `CompressibleEulerPotentialTemperatureEquationsWithGravity3D`, the
  `CompressibleEulerEnergyEquationsWithGravity2D`, the
  `CompressibleEulerEnergyEquationsWithGravity3D`, the
  `CompressibleMoistEulerEquations2D`, and the `CompressibleRainyEulerEquations2D`. The
  methods of `boundary_condition_slip_wall` with the additional argument `direction`,
  which fixed the sign convention of `StructuredMesh`, are removed for the
  `CompressibleMoistEulerEquations2D` and the `CompressibleRainyEulerEquations2D`. Use the
  methods based on `normal_direction` together with `P4estMesh` instead. The physical flux
  `flux(u, orientation, equations)` is kept since `VolumeIntegralWeakForm` uses it also on
  `P4estMesh`. The 1D equations and the equations in covariant form are not affected.

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
