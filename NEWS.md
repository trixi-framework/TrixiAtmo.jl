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
- The `p4est` connectivities of `P4estMeshCubedSphere2D`, `P4estMeshQuadIcosahedron2D`, and
  `P4estMeshCubedSphereTopography` now contain all corners (and edges in 3D) connecting
  trees, as required by `p4est` for 2:1 balancing. Previously, they only contained the face
  connections, so the 3D meshes could violate the 2:1 balance after local refinement
  (e.g., with AMR). In 2D, this was hidden by Trixi.jl calling `p4est_balance` twice,
  which is no longer necessary with complete connectivities. Simulations with uniform
  refinement are not affected.
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
