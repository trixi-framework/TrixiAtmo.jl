# 🌎 TrixiAtmo.jl 🌍

[comment]: <> ([![Docs-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://trixi-framework.github.io/TrixiAtmo.jl/stable))  
[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://trixi-framework.github.io/TrixiAtmo.jl/dev)
[![Slack](https://img.shields.io/badge/chat-slack-e01e5a)](https://join.slack.com/t/trixi-framework/shared_invite/zt-sgkc6ppw-6OXJqZAD5SPjBYqLd8MU~g)
[![Build Status](https://github.com/trixi-framework/TrixiAtmo.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/trixi-framework/TrixiAtmo.jl/actions?query=workflow%3ACI)
[![Codecov](https://codecov.io/gh/trixi-framework/TrixiAtmo.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/trixi-framework/TrixiAtmo.jl)
[![Coveralls](https://coveralls.io/repos/github/trixi-framework/TrixiAtmo.jl/badge.svg?branch=main)](https://coveralls.io/github/trixi-framework/TrixiAtmo.jl?branch=main)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[comment]: <> ([![DOI](https://zenodo.org/badge/DOI/TODO.svg)](https://doi.org/TODO))  

<p align="center">
  <img width="60%" src="https://trixi-framework.github.io/assets/logo_atmo.png">
</p>

**Note: This repository is still in its alpha stage and anything might change at
any time and without warning.**

**TrixiAtmo.jl** is a numerical simulation package focused on atmospheric flows. It builds
upon [Trixi.jl](https://github.com/trixi-framework/Trixi.jl), a generic flow solver for
conservation laws, implementing discontinuous Galerkin methods and written in Julia.

Currently available features include:

* Compressible Euler and shallow water models on cubed sphere meshes, discretizing the
  atmosphere or its two-dimensional surface
* Moist compressible Euler equations, including cloud and rain microphysics
* Flux-differencing formulations, including entropy-stable schemes

## Installation

If you have not yet installed Julia, please [follow the instructions for your
operating system](https://julialang.org/downloads/platform/). TrixiAtmo.jl works
with Julia v1.10 and newer. We recommend using the latest stable release of Julia.

TrixiAtmo.jl is **not** a registered Julia package yet, and therefore needs to be
downloaded manually and then run from within the cloned directory:
```bash
git clone https://github.com/trixi-framework/TrixiAtmo.jl.git
julia --project=@.
```
In addition TrixiAtmo.jl requires the numerical solver framework
[Trixi.jl](https://github.com/trixi-framework/Trixi.jl), relevant sub-packages of
[OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) for time integration, and
[Plots.jl](https://github.com/JuliaPlots/Plots.jl) for visualization, which can be
installed by executing the following in the Julia REPL:
```julia
julia> using Pkg

julia> Pkg.add(["Trixi", "Trixi2Vtk", "OrdinaryDiffEqLowStorageRK", "OrdinaryDiffEqSSPRK", "Plots"])
```


## Usage

In the Julia REPL, first load the package Trixi.jl
```julia
julia> using Trixi
```
Then start a simulation by executing
```julia
julia> trixi_include("examples/elixir_euler_warmbubble.jl")
```
Please see our documentation for more advanced setups.


## Authors
TrixiAtmo.jl is maintained by the
[Trixi authors](https://github.com/trixi-framework/Trixi.jl/blob/main/AUTHORS.md).
It was initiated by
[Andrés Rueda-Ramírez](https://andres.rueda-ramirez.com)
(Polytechnic University of Madrid (UPM), Spain),
[Benedict Geihe](https://www.mi.uni-koeln.de/NumSim/), and
[Tristan Montoya](https://tjbmontoya.com/)
(University of Cologne, Germany).

## License and contributing
TrixiAtmo.jl is licensed under the MIT license (see [LICENSE.md](LICENSE.md)).
Since TrixiAtmo.jl is an open-source project, we are very happy to accept contributions from the
community. To get in touch with the developers,
[join us on Slack](https://join.slack.com/t/trixi-framework/shared_invite/zt-sgkc6ppw-6OXJqZAD5SPjBYqLd8MU~g)
or [create an issue](https://github.com/trixi-framework/TrixiAtmo.jl/issues/new).

## Acknowledgments
<p align="center">
  <!-- BMBF -->
  <img width="200px" src="https://user-images.githubusercontent.com/3637659/231436391-b28a76a4-f027-40f9-bd28-14e3a2f3e16a.png"/>
</p>

This project has benefited from funding from the 
[Federal Ministry of Education and Research](https://www.bmbf.de) (BMBF) 
through the following grants:

* Project grant "Adaptive earth system modeling with significantly reduced computation time for
  exascale supercomputers (ADAPTEX)" (funding id: 16ME0668K)
* Project grant "ICON-DG" of the WarmWorld initiative (funding id: 01LK2315B)
