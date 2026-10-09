# AGENTS.md

Canonical instructions for AI coding agents (Claude Code, Codex, Copilot,
Cursor, Gemini CLI, ...) working on TrixiAtmo.jl. Most tools read this file
natively; Claude Code reads it via the `@AGENTS.md` import in `CLAUDE.md`. Keep
this file the single source of truth; put only tool-specific notes into
tool-specific files. TrixiAtmo.jl follows the conventions of
[Trixi.jl](https://github.com/trixi-framework/Trixi.jl), whose developer
documentation (in particular the pages on conventions, style guide, and testing)
also applies here; this file summarizes what agents need most.

When editing this file, only add instructions and commands that have actually
been tested in this repository. Instructions that do not work are worse than no
instructions.

TrixiAtmo.jl is a Julia package for numerical simulations of atmospheric flows
built on top of Trixi.jl, which provides adaptive high-order discontinuous
Galerkin spectral element methods (DGSEM) for conservation laws. TrixiAtmo.jl
adds atmospheric models (compressible Euler equations with gravity in several
formulations, moist and rainy Euler equations, shallow water equations), meshes
of the sphere, and discretizations on two-dimensional manifolds embedded in 3D,
including a covariant formulation. Correctness of the numerics matters more than
anything else: many bugs are silent (e.g., a wrong index or sign that only
reduces accuracy), so be careful and verify mathematical changes.

Before editing, inspect `git status` and preserve all existing work. Do not
discard, reset, or stash unrelated changes while testing a fix.

## Repository layout

- `src/`: package code (see [Architecture](#architecture)).
- `examples/`: elixirs, i.e., scripts setting up and running a simulation,
  grouped by physics, e.g., `advection/{cartesian,covariant}`,
  `euler/{dry_air,moist_air,precipitation}/{buoyancy,global_circulation,mountain_flow,tests}`,
  and `shallow_water/{cartesian,covariant}`. Here, `cartesian` and `covariant`
  refer to the formulation on the sphere (see below). `examples_dir()` returns
  the path to this directory.
- `test/`: test suite (entry point `test/runtests.jl`, test items in
  `test/test_*.jl`, shared test helpers in `test/test_trixiatmo.jl`).
- `docs/`: documentation source (Documenter.jl). The API reference
  (`docs/src/reference.md`) collects all docstrings automatically; tutorials are
  Literate.jl scripts in `docs/src/tutorials/`.
- `utils/`: developer scripts (setup of the agent environments, see below).
- `run/` (if present, ignored by git): the human developer's own Julia project.
  Do not modify it; agents use `run_agents/` instead (see below).

## Architecture

TrixiAtmo.jl extends Trixi.jl via multiple dispatch. As in Trixi.jl, a
simulation combines **equations × mesh × solver** into a **semidiscretization**,
which is integrated in time. TrixiAtmo.jl does not define its own abstract type
hierarchy but subtypes the abstract types of Trixi.jl (e.g.,
`Trixi.AbstractEquations`) and adds methods to (often internal) functions of
Trixi.jl such as `Trixi.rhs_hyperbolic!`, `Trixi.create_cache`,
`Trixi.calc_volume_integral!`, or `Trixi.max_dt`. Thus, changes of Trixi.jl
internals can break TrixiAtmo.jl. When extending such a function, look up its
current signature and implementation in the version of Trixi.jl that is actually
used (e.g., `julia --project=run_agents -e 'using Trixi; print(pkgdir(Trixi))'`);
never edit the installed Trixi.jl sources. Trixi.jl runs the TrixiAtmo.jl test
items tagged `upstream` in its own downstream CI.

- The module file `src/TrixiAtmo.jl` contains all imports and exports. Its
  include order is deliberate (`auxiliary` → `equations` → `meshes` →
  `semidiscretization` → `solvers` → `callbacks_step` → `callbacks_stage`);
  respect it when adding cross-cutting types.
- **Equations** (`src/equations/`) are concrete types defined mainly by
  *pointwise* functions dispatched on them: `flux`, numerical fluxes
  (`flux_ec`, `flux_tec`, `flux_etec`, `FluxLMARS`, nonconservative fluxes such
  as `flux_nonconservative_souza_etal`, ...), wave speed estimates,
  `initial_condition_*`, `boundary_condition_*`, `source_terms_*`, and variable
  conversions (`cons2prim`, `cons2entropy`, ...). There are
  - equations in Cartesian coordinates, e.g.,
    `CompressibleEulerPotentialTemperatureEquations{1,2,3}D` (with
    `...WithGravity{1,2,3}D` variants including the geopotential as an
    auxiliary variable), `CompressibleEulerEnergyEquationsWithGravity{2,3}D`,
    `CompressibleMoistEulerEquations2D`, `CompressibleRainyEulerEquations2D`,
    and `ShallowWaterEquations3D` (shallow water equations on the sphere in 3D
    Cartesian coordinates);
  - equations in covariant form on manifolds, subtypes of
    `AbstractCovariantEquations{NDIMS, NDIMS_AMBIENT, GlobalCoordinateSystem, NVARS}`,
    e.g., `CovariantLinearAdvectionEquation2D`,
    `CovariantShallowWaterEquations2D`, `SplitCovariantShallowWaterEquations2D`,
    and `CovariantEulerEnergyEquations2D`. Their prognostic variables and fluxes
    use contravariant vector components, and their pointwise functions take an
    additional argument `aux_vars` with the geometric information (metric
    terms) at each node, e.g., `flux(u, aux_vars, orientation, equations)`.
    `GlobalCoordinateSystem` is `GlobalCartesianCoordinates()` or
    `GlobalSphericalCoordinates()`; initial conditions are given in this system
    and converted via `transform_initial_condition(initial_condition, equations)`.
- **Meshes** (`src/meshes/`) of the sphere based on Trixi.jl meshes:
  `P4estMeshCubedSphere2D`, `P4estMeshQuadIcosahedron2D`,
  `DGMultiMeshTriIcosahedron2D`, and `P4estMeshCubedSphereTopography` (3D shell
  with terrain-following coordinates, e.g., `GalChen` and `Sleve`).
- **Semidiscretization** (`src/semidiscretization/`): metric terms
  (`MetricTermsCrossProduct`, `MetricTermsInvariantCurl`, `MetricTermsCovariant`
  with Christoffel symbols) and specializations of
  `SemidiscretizationHyperbolic` for 2D manifolds in 3D.
- **Solvers** (`src/solvers/`): specializations of the Trixi.jl solvers for
  `dgsem_p4est` (2D manifolds in 3D, in Cartesian and covariant form) and
  `dgmulti` (covariant form on triangles, including flux differencing).
- **Callbacks**: `src/callbacks_step/` contains analysis, output, and time step
  specializations (e.g., for the covariant form); `src/callbacks_stage/`
  contains `NonlinearSolveDG`, which solves the nonlinear system of the rainy
  Euler equations within Runge-Kutta stages.
- Physical constants such as `EARTH_RADIUS`, `EARTH_GRAVITATIONAL_ACCELERATION`,
  `EARTH_ROTATION_RATE`, and `SECONDS_PER_DAY` are exported; use them instead of
  literal values.
- **Elixirs** (`examples/`) are the main user-facing entry point and are run via
  `trixi_include`. Every elixir must be used in at least one test.

## Testing: never run the full test suite

**Never run the complete test suite of TrixiAtmo.jl.** It takes more than 20
minutes and is run in CI (GitHub Actions) for every pull request. In
particular, do **not** run

- `Pkg.test("TrixiAtmo")` or `julia test/runtests.jl` without a restrictive
  setup (without `TRIXI_TEST`, they run **all** tests, including the MPI and
  threaded tests),
- `TRIXI_TEST=all` or `TRIXI_TEST=upstream`,
- whole test files such as `include("test/test_3d_potential_temperature.jl")`,
- the special suites `mpi` and `threaded`, which relaunch Julia via `mpiexec`
  or with multiple threads.

Instead, run only the individual test items related to the code you changed.
Tests are isolated `@testitem`s (TestItems.jl) run with TestItemRunner.jl. Test
files are named like `test_<dim>_<equation>.jl` (e.g.,
`test_2d_moist_euler.jl`) plus files such as `test_spherical_advection.jl`,
`test_unit.jl` (consistency checks of fluxes), `test_type.jl` (type
stability), and `test_trixi_consistency.jl` (comparison with Trixi.jl). Test
item names start with a prefix for the file (e.g., `"Moist Euler 2D: ..."`).
Each test item lists its setup snippets via `setup = [...]` (e.g., `Setup`,
defined in `test/runtests.jl` and providing the helpers from
`test/test_trixiatmo.jl`, and file-specific snippets such as `MoistEuler2D`,
defined at the top of the corresponding test file and setting `EXAMPLES_DIR`)
and its test group via `tags = [...]`. The value of the environment variable
`TRIXI_TEST` selects the test items with the corresponding tag (see
`test/runtests.jl`); CI runs all tests (`.github/workflows/ci.yml`).

### The `run_agents` Julia project

Use a local Julia project in the directory `run_agents/` at the root of this
repository for running tests, elixirs, and scratch scripts. It develops the
local version of TrixiAtmo.jl (`Pkg.develop`) and contains all dependencies of
the tests (copied from `test/Project.toml`, including compat bounds). Recycle it
from one run (and session) to the next instead of creating new environments.
The directory is ignored by git (`run*/` in `.gitignore`) and must never be
committed.

The formatter JuliaFormatter.jl v1.0.60 (see [Code style](#code-style)) cannot
be installed in `run_agents/` since its dependencies are incompatible with those
of Trixi.jl. Thus, it lives in the separate project `run_agents/formatter/`,
where it is pinned to this version.

Create or update both projects by running the script
`utils/setup_run_agents.jl` **before running tests or the formatter**, e.g.,
from the root of the repository,

```bash
julia utils/setup_run_agents.jl
```

The script works from any directory and is idempotent: if `run_agents/` is up to
date with `test/Project.toml` and still points to this checkout, it reuses the
environment; otherwise, it recreates the project. Similarly, it only
(re-)creates `run_agents/formatter/` if JuliaFormatter.jl v1.0.60 is not pinned
there. The first setup may take several minutes; afterward, it takes only a few
seconds unless packages need to be precompiled again.

The setup script does not check the dependencies of TrixiAtmo.jl itself (in the
root `Project.toml`). If they changed (e.g., after you added a dependency to
TrixiAtmo.jl or pulled such a change from `main`), loading TrixiAtmo.jl fails
with an error such as
`ArgumentError: Package TrixiAtmo does not have X in its dependencies`.
Only in this case, run
`julia --project=run_agents -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'`.
Do not run `Pkg.resolve()` routinely, since it takes some time.

If you need additional packages for your own experiments (e.g., BenchmarkTools),
add them with `julia --project=run_agents -e 'using Pkg; Pkg.add("PackageName")'`.
Note that they are removed when the project is re-created after a change of
`test/Project.toml`; never add them to `test/Project.toml` or `Project.toml`
unless this is part of the actual task. Temporary scripts can be stored in
`run_agents/` as well.

`run_agents/Manifest.toml` pins resolved package versions (in particular of
Trixi.jl), and the setup script does not upgrade them. CI uses the latest
compatible versions, so the local environment can fall behind. When
investigating CI failures or numerical differences, compare the package
versions in the CI logs with
`julia --project=run_agents -e 'using Pkg; Pkg.status()'` before attributing a
change to TrixiAtmo.jl. To move to the latest compatible versions, run
`julia --project=run_agents -e 'using Pkg; Pkg.Registry.update(); Pkg.update()'`.
Use a separate ignored or temporary environment to test other dependency
versions so that the usual `run_agents/` setup remains reproducible.

### Running selected test items

From the root of the repository, run, e.g., the following command (or, if
available, the corresponding Julia code in a persistent session via MCP, see
[below](#optional-persistent-julia-sessions-via-mcp)):

```bash
julia --project=run_agents -e '
using TestItemRunner
cd("test") do
    @run_package_tests filter = ti -> ti.name in (
        "Moist Euler 2D: elixir_gemein_bubble dry",
        "Unit: Consistency check for 3D shallow water fluxes: SWE",
    )
end'
```

Useful filters (`ti` has the fields `name`, `filename`, and `tags`):

- exact names: `ti -> ti.name in ("name 1", "name 2")`
- substring: `ti -> occursin("elixir_hoeck_bubble", ti.name)`
- one file (use sparingly, some files take minutes):
  `ti -> endswith(ti.filename, "test_2d_rainy_euler.jl")`
- one tag: `ti -> :unit_fluxes in ti.tags`

Notes:

- **Check the test summary.** A filter that matches no test item does not fail;
  TestItemRunner just prints a summary without any tests (`None`). Make sure
  that the number of passed tests is plausible, and copy test item names
  exactly from the `@testitem "..."` lines.
- If a persistent Julia session is available, prefer it (see
  [the next section](#optional-persistent-julia-sessions-via-mcp)). Otherwise,
  batch related test items into one call to avoid repeated startup costs.
- Always call `@run_package_tests` inside `cd(".../test") do ... end` as
  above. The macro searches for test items in the parent directory of the file
  it is called from; outside of a file (e.g., in `julia -e` or an MCP
  session), this is the parent of the current working directory. From the
  repository root, it would search the parent directory of the repository
  (e.g., failing with `UndefVarError: TrixiAtmo not defined`). The `do` block
  restores the previous working directory afterward. Output files of test
  items are written to `test/out/`, as in `Pkg.test`.
- Do not use `julia --project=.` (the package environment lacks the test
  dependencies) or `julia --project=test` (the test environment lacks
  TrixiAtmo.jl), and do not stack environments via `JULIA_LOAD_PATH` (e.g.,
  `"@:test:@stdlib"`); the latter can mix incompatible package versions and fail
  during precompilation. Use the `run_agents` project instead.
- The first run after setting up or updating `run_agents` may take several
  minutes because of precompilation. Later runs are much faster.
- Find the relevant test items by searching `test/` for the elixir name,
  function name, or equation type you changed, e.g.,
  `rg -n "elixir_hoeck_bubble" test/` (or `grep -rn`). Remember that elixirs
  can be included by other elixirs and tests via `trixi_include`.
- Some reference values differ slightly on macOS. The affected test items use
  a larger tolerance only on GitHub Actions (`CI_ON_MACOS` in
  `test/test_trixiatmo.jl`), so they can fail locally on macOS (e.g.,
  `"Euler potential temperature 3D: elixir_potential_temperature_held_suarez"`).
  Compare with the behavior on the unchanged code before attributing such a
  failure to your change.
- To quickly try an elixir without the test harness, use the `run_agents`
  project (elixirs need packages such as `OrdinaryDiffEqLowStorageRK`, which are
  test-only dependencies) and `trixi_include`, which can override assignments in
  the elixir, e.g., in `julia --project=run_agents`:
  ```julia
  using Trixi, TrixiAtmo
  trixi_include(joinpath(TrixiAtmo.examples_dir(), "euler", "dry_air", "buoyancy",
                         "elixir_gemein_bubble.jl"),
                tspan = (0.0, 1.0))
  ```
  Use the qualified name `TrixiAtmo.examples_dir()` since Trixi.jl exports a
  function `examples_dir` as well. Only variables that are assigned in the
  elixir can be overridden (e.g., not every elixir defines `tspan`). Output
  files are written to `out/`, which must not be committed.
- When developing new kernels, run Julia with `--check-bounds=yes` so that
  out-of-bounds accesses are caught despite `@inbounds`.

### Optional: persistent Julia sessions via MCP

Some developers configure MCP servers providing persistent Julia sessions, in
which packages are loaded and compiled only once. **Only if the corresponding
tools are available to you**, use them to run Julia code such as tests and
elixirs, in the following order of preference:

1. a [Kaimon.jl](https://github.com/kahliburke/Kaimon.jl) session started by
   the developer that loads TrixiAtmo.jl from this checkout (see
   [below](#kaimonjl)), otherwise
2. a session of [julia-mcp](https://github.com/aplavin/julia-mcp) (see
   [below](#julia-mcp)), otherwise
3. the shell commands described above.

Everything else in this file applies unchanged. Run the `run_agents` setup
script in the shell in any case (e.g., for the formatter), even if the
description of an MCP tool tells you to never run Julia via the command line.
Sessions may run with non-default Julia flags chosen by the developer (e.g.,
`--threads=1 --check-bounds=yes`). Check them with `Threads.nthreads()` and
`Base.JLOptions().check_bounds` (`1` means bounds checking is enabled) if they
matter, e.g., for timings or allocation tests.

#### Kaimon.jl

Kaimon.jl connects agents to Julia sessions (REPLs) started by the developer,
who sees the code you run in their REPL. Its tools include `ping` and `ex`
(in Claude Code, e.g., `mcp__kaimon__ex`). They are only available while the
Kaimon.jl server is running.

- Call `ping` to list the connected sessions with their 8-character keys. Use
  a session only if
  `ex(e = "using TrixiAtmo; pkgdir(TrixiAtmo)", q = false, ses = "<key>")`
  returns the absolute path of the root of this repository. Otherwise (e.g.,
  for sessions of other packages such as Trixi.jl, clones, or worktrees, or if
  no session is connected), use julia-mcp or the shell instead.
- Always pass `q = false` (the default `q = true` does not return the result)
  and `ses = "<key>"`.
- `ex` returns the value of the last expression or the error, together with
  the output printed to `stdout` and `stderr` (e.g., test summaries and
  details of failures). However, output capture may be bypassed for code that
  loads packages via `using` or `import`. Thus, load packages in a separate
  call first, e.g., `ex(e = "using TestItemRunner", q = false, ses = "<key>")`,
  and then run test items in another call such as
  ```julia
  cd("/path/to/TrixiAtmo.jl/test") do
      @run_package_tests filter = ti -> ti.name in (
          "Moist Euler 2D: elixir_gemein_bubble dry",
          "Unit: Consistency check for 3D shallow water fluxes: SWE",
      )
  end;
  ```
  ending with `;` to avoid printing the returned test set.
- Pass `max_output = 25000` (the maximum; the default is 6000 characters) when
  running test items. Longer output is truncated in the middle, keeping its
  beginning and end (including the test summary). Since a single test item
  running an elixir can print more than 10000 characters, run only a few test
  items per call and rerun failed test items individually if their details
  were truncated.
- Kaimon.jl may remove calls such as `println` from your code; use the value
  of the last expression to return results instead.
- The session belongs to the developer. Do not restart or shut it down, do not
  add or remove packages, do not start new sessions, and avoid persistent
  changes to its global state (e.g., use `cd(...) do ... end` instead of
  `cd(...)`). If the session cannot be used (e.g., missing packages, Revise.jl
  not loaded so that changes in `src/` are not picked up, or a restart is
  required after changing type definitions), ask the developer or use
  julia-mcp or the shell instead.
- Do not use the Kaimon.jl tools `run_tests` (not designed for the test items
  of TrixiAtmo.jl, may run the complete test suite) and `format_code` (see
  [Code style](#code-style) for the formatter).
- If an evaluation takes longer than about 30 seconds (e.g., the first run of
  test items, which includes compilation), `ex` returns a job ID (`eval_id`)
  before the evaluation is finished. In this case, wait at least 30 seconds and
  call `check_eval(eval_id = "<id>")` until its status is `completed` or
  `failed`; then, the output and the result are shown.
  Do not start other evaluations in the same session while a job is running.
- If an evaluation does not produce any output for 10 minutes, Kaimon.jl
  reports a timeout, but the evaluation may still be running in the session.
- `trixi_include` writes its output files to `out/` in the working directory of
  the session (check it with `pwd()`).

#### julia-mcp

julia-mcp provides the tools `julia_eval`, `julia_restart`, and
`julia_list_sessions` (in Claude Code, e.g., `mcp__julia__julia_eval`) and
starts the Julia sessions itself.

- Always pass the **absolute path** of the `run_agents` directory as
  `env_path` (e.g., `/path/to/TrixiAtmo.jl/run_agents`). Do not call
  `Pkg.activate` in the code.
- `julia_eval` returns only the output printed to `stdout` and `stderr`, not
  the value of the last expression. Use `println(...)` or `display(...)` to see
  results, e.g., `println(format([...]))` for the formatter.
- The working directory of the session is `run_agents/`, not the repository
  root. Use absolute paths, e.g., wrap `@run_package_tests` in
  `cd("/path/to/TrixiAtmo.jl/test") do ... end` to run test items.
- The default `timeout` is 60 seconds, and a call exceeding its timeout **kills
  the session**. Pass a generous `timeout` (e.g., `1800`) for the first call
  (`using Trixi, TrixiAtmo, TestItemRunner`) and for the first run of test
  items or elixirs, which include compilation.
- Afterward, repeat `@run_package_tests` calls with different filters in
  the same session. Revise.jl is loaded automatically (if it is installed in
  the global environment), so changes to files in `src/` are picked up without
  restarting; test files are re-read by TestItemRunner.
- Call `julia_restart` with the same `env_path` after changing type
  definitions (e.g., fields of a `struct`), which Revise cannot handle, after
  re-creating or updating `run_agents`, or if the session is in a broken state.
- `trixi_include` writes its output files to `out/` in the working directory,
  i.e., to `run_agents/out/`.
- Run the formatter (see [Code style](#code-style)) in a separate session with
  the absolute path of `run_agents/formatter` as `env_path`, e.g.,
  `using JuliaFormatter; format(["/path/to/TrixiAtmo.jl/src/file.jl"])`, since
  JuliaFormatter.jl is not available in the `run_agents` environment. Use
  absolute paths since the working directory is `run_agents/formatter/`.

### Adding and changing tests

- New or changed behavior must be covered by focused tests; the overall test
  coverage (reported by Codecov and Coveralls) must not decrease. For a bug
  fix, verify that the regression test detects the original behavior when
  feasible, using a separate checkout or temporarily reverting only your own
  changes. Preserve unrelated work in the shared working tree.
- Add tests as `@testitem`s to the matching `test/test_xxx.jl` file (e.g., an
  elixir for the 2D moist Euler equations goes into
  `test/test_2d_moist_euler.jl`, flux consistency checks into
  `test/test_unit.jl`). Copy the name prefix, `setup = [...]`, and
  `tags = [...]` of neighboring test items. New test files are discovered
  automatically; give their test items a new tag (and a new snippet setting
  `EXAMPLES_DIR` if needed). Add the tag `:upstream` only to a few
  representative test items that should also run in the downstream CI of
  Trixi.jl.
- Exercise every new elixir in at least one test item. For a new equation or
  numerical method, run a convergence study and report the observed order.
- Test item names must be unique within a file.
- Test items should run fast (at most a few seconds after compilation). Use
  short `tspan`s and coarse meshes (e.g., small `cells_per_dimension`,
  `trees_per_cube_face`, or `initial_refinement_level`) where possible.
- Regression tests via `@test_trixi_include` (defined in
  `test/test_trixiatmo.jl`) compare `l2`/`linf` errors with reference values.
  If a change intentionally modifies results, inspect the actual values and
  relevant numerical behavior first, and check that the change is plausible;
  for example, errors of smooth test cases with an analytical solution (e.g.,
  convergence tests with manufactured solutions) should usually not increase
  after a bug fix. Otherwise, the errors are computed with respect to the
  initial condition evaluated at the final time, which is in general not the
  exact solution, so they may change in either direction. Then run the affected
  test items once, take the actual values from the
  `Evaluated: isapprox(expected, actual; ...)` lines of the failures, replace
  the reference values **only within the affected test items**, rerun them, and
  explain the change in the PR.
- Check allocations of the right-hand side after running an elixir with
  `@test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)` as in
  the existing test items.
- Never loosen tolerances (`atol`, `rtol`) just to make a failing test pass.
  A local, commented tolerance is acceptable only for small differences whose
  cause is understood and outside TrixiAtmo.jl (e.g., platform-dependent
  floating-point differences as for `CI_ON_MACOS`); mention it to the user.
- When adding a numerical flux or wave speed estimate, add unit tests to
  `test/test_unit.jl` (see the existing `"Unit: Consistency check for ..."`
  test items). For a conservative numerical flux `f`, check consistency with
  the physical flux, i.e.,
  `f(u, u, orientation, equations) ≈ flux(u, orientation, equations)` and the
  same with `normal_direction` (with `aux_vars` arguments for covariant
  equations). Check that the 1D, 2D, and 3D versions of an equation agree for
  corresponding states and directions. For all fluxes and wave speed estimates
  with `normal_direction` methods, check agreement with the `orientation`
  methods for axis-aligned unit normals.
- Except for consistency checks, unit tests for numerical fluxes and wave speed
  estimates should use *different* left and right states: many copy-paste bugs
  are invisible for `u_ll == u_rr`.
- When adding an equation or new pointwise functions, add type stability tests
  for `Float32` and `Float64` to `test/test_type.jl` (see the existing
  `"Type stability: ..."` test items using `@inferred`).

## Code style

- Format all changed Julia files with JuliaFormatter.jl v1.0.60 (configuration
  in `.JuliaFormatter.toml`, SciML style); CI checks the formatting of the
  whole repository. Other versions of JuliaFormatter.jl format differently, so
  use the pinned version in `run_agents/formatter/` (created by the
  [setup command](#the-run_agents-julia-project)), e.g., from the root of the
  repository,
  ```bash
  julia --project=run_agents/formatter -e 'using JuliaFormatter; format(["src/file.jl", "test/test_file.jl"])'
  ```
  or via julia-mcp if available (see [above](#julia-mcp)). `format` returns
  `true` if all files were already formatted; otherwise, it formats them in
  place and returns `false`. Pass only the files you changed.
- Maximum line length 92, 4 spaces for indentation (never tabs), `snake_case`
  for functions and variables, `CamelCase` for modules and types, trailing `!`
  for mutating functions. Use ASCII for new identifiers in executable code;
  comments and docstrings may use Unicode.
- Argument order: the modified argument first (e.g., `du`, or `cache` if only
  the cache is modified), otherwise `mesh, equations, solver, cache`; dispatch
  helpers before the general argument (e.g.,
  `have_nonconservative_terms(equations), equations`). For covariant equations,
  `aux_vars` follows the state(s), e.g., `flux(u, aux_vars, orientation, equations)`
  or `flux_ec(u_ll, u_rr, aux_vars_ll, aux_vars_rr, orientation, equations)`.
- Name related functions with the general category first, e.g.,
  `flux_central` (not `central_flux`), `source_terms_coriolis`.
- Use a suffix `_` (`name_`) for local variables that would otherwise shadow
  existing names and a prefix `_` (`_name`) for fragile internal functions that
  are not part of the documented API.
- Prefer `for i in ...` over `for i = ...`. Group multiline expressions with
  explicit parentheses instead of relying on implicit line continuation.
- Most source files are wrapped in `@muladd begin ... end` (with
  `#! format: noindent`); keep this structure.
- Extend functions of Trixi.jl explicitly (e.g., `function Trixi.max_dt(...)`
  or via the `import Trixi: ...` list in `src/TrixiAtmo.jl`) instead of
  defining new functions with the same name. New exported names go into the
  export lists in `src/TrixiAtmo.jl`.
- Keep code type-stable and generic in the real type (`Float32` and `Float64`
  are supported, as well as AD types): use `0.5f0` for exactly representable
  constants and `convert(RealT, pi)` with `RealT = eltype(u)` (or `eltype(x)`)
  otherwise. Use StaticArrays.jl (`SVector`, `MVector`, ...) for small arrays of
  fixed size.
- Performance-critical code must not allocate; many regression tests check this
  with `@test_allocations`.
- Document new public functions and types with docstrings (they appear in the
  API reference automatically); for internal code, a top-level comment can
  suffice. Reference publications where appropriate. Comments should explain
  *why*, not only *what*.
- Match the style of the surrounding code, and look for existing helpers in
  TrixiAtmo.jl and Trixi.jl before adding new ones.
- In elixirs, load packages via `using OrdinaryDiffEqLowStorageRK, Trixi, TrixiAtmo`
  (or the OrdinaryDiffEq sub-package needed) and use the standard variable
  names `polydeg`, `surface_flux`, `volume_flux`, and `cells_per_dimension`,
  `trees_per_cube_face`, or `initial_refinement_level` for the resolution, so
  that `trixi_include` and `convergence_test` can modify them. With a
  `StepsizeCallback`, pass an integer dummy time step `dt = 1` (not `1.0`) to
  `solve`.

## Numerical conventions and common pitfalls

- `orientation`s: `1 => x, 2 => y, 3 => z`; `direction`s:
  `1 => -x, 2 => +x, 3 => -y, 4 => +y, 5 => -z, 6 => +z`. For covariant
  equations, orientations refer to the local reference coordinates of the
  manifold, and vector quantities are stored as contravariant components.
- A `normal_direction` argument may be non-unit and carry surface scaling,
  especially on curved meshes such as the cubed sphere. Fluxes and wave speeds
  using it must account for its magnitude exactly once. Numerical fluxes and
  wave speed estimates must be homogeneous of degree one in
  `normal_direction`, i.e., `f(u_ll, u_rr, α * n) == α * f(u_ll, u_rr, n)` for
  `α > 0` (e.g., multiply sound speeds by `norm(normal_direction)`). Check
  agreement with the corresponding `orientation` method for axis-aligned unit
  normals.
- Nonconservative surface terms depend on the argument order: the element
  receiving the flux passes its own state first,
  `nonconservative_flux(u_own, u_neighbor, ...)`.
- Some equations carry auxiliary variables as part of the state vector, e.g.,
  the geopotential in the `...WithGravity` equations or the bottom topography in
  `ShallowWaterEquations3D`. Fluxes must leave them unchanged (zero flux), and
  well-balancedness tests (`elixir_*well_balanced*`) check that steady states
  are preserved.
- Arrays visible to the time integrator have the suffix `_ode` (`u_ode`,
  `du_ode`); when passed with `semi`, normally wrap them with
  `wrap_array(u_ode, semi)` before internal use. Wrapped arrays have no suffix.
- Compare 1D/2D/3D and x/y/z implementations of the same function when
  changing one of them (e.g., the `CompressibleEulerPotentialTemperature...`
  equations exist in 1D, 2D, and 3D, with and without gravity); copy-paste
  errors between dimensions and directions are the most common source of bugs.

## Pull requests and documentation

- Develop in a branch; branch names usually follow
  `<developer initials>/<short description>` (e.g., `hr/fix_mortar_fluxes`).
- Keep changes focused on a single goal; PRs should usually change at most about
  500 lines.
- If a PR is created based on major contributions from an LLM/AI tool, disclose
  the LLM/AI assistance in the PR description.
- Do not commit generated files (e.g., `out*/` directories, `.vtu`/`.h5`
  output, `Manifest.toml` files) or temporary notes.
- Do not commit or push unless asked to do so by the user.
