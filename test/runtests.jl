using TestItemRunner

# We run subsets of the test suite in parallel CI jobs by setting the `TRIXI_TEST`
# environment variable. Its value is matched against the `tags` attached to each
# `@testitem`. By default (`TRIXI_TEST == "all"`), we run all tests.
const TRIXI_TEST = get(ENV, "TRIXI_TEST", "all")

# With `TRIXI_TEST_VERBOSE=true`, `@run_package_tests` prints every `@testitem` (together
# with its run time) in the final test summary instead of only the failing ones. It is off
# by default (to keep local runs quiet) and enabled in CI so that per-testitem timings
# show up in the job logs.
const TRIXI_TEST_VERBOSE = get(ENV, "TRIXI_TEST_VERBOSE", "false") == "true"

# The MPI tests check the partitioning of the mesh for exactly this number of ranks.
const TRIXI_MPI_NPROCS = 4
const TRIXI_NTHREADS = clamp(Sys.CPU_THREADS, 2, 3)

# A few suites cannot run in the ordinary in-process `@run_package_tests` model:
# they need a specially-launched Julia process - `mpi` (multiple ranks via
# `mpiexec`) and `threaded` (multiple threads via `--threads`). We handle them by
# *relaunching* Julia/`mpiexec` on this very file: the launched worker re-enters
# `runtests.jl` with `TRIXI_TEST_RUN_ITEMS` set and then runs the tag-filtered test
# items in-process (`TestItemRunner` evaluates items in the current process, so this
# also works for every MPI rank). The `TRIXI_TEST` value of such a suite selects
# the items via the equally-named tag.
const SPECIAL_PROCESS_SUITES = ("mpi", "threaded")
const IN_WORKER = haskey(ENV, "TRIXI_TEST_RUN_ITEMS")

# Remove the output directory `out`, where examples write solution/restart/mesh
# files. We do this once at the start of a run *and* register it to run again when
# the process exits (`atexit`, so it also fires when a test set fails and throws),
# leaving a clean working tree behind. Only the parent process handles this: workers
# share the working directory and finish before the parent exits, so a single
# cleanup here avoids races on `out`.
function clean_outdir()
    # When running tests via `Pkg.test`, the working directory is the
    # `test` directory of the package to be tested.
    outdir = joinpath(@__DIR__, "out")
    isdir(outdir) && rm(outdir, recursive = true, force = true)
end
if !IN_WORKER
    clean_outdir()
    atexit(clean_outdir)
end

# Special suites this (parent) invocation will dispatch into their own processes.
const SUITES_TO_DISPATCH = if IN_WORKER
    String[]
elseif TRIXI_TEST == "all"
    collect(SPECIAL_PROCESS_SUITES)
elseif TRIXI_TEST in SPECIAL_PROCESS_SUITES
    [TRIXI_TEST]
else
    String[]
end

# `import` is only allowed at top level, so load here what `dispatch_special_suite`
# needs: `MPI` for `mpiexec`.
if "mpi" in SUITES_TO_DISPATCH
    import MPI
end

# Relaunch Julia/`mpiexec` for a suite that needs a special process. The worker
# re-enters this file with `TRIXI_TEST_RUN_ITEMS=true` and `TRIXI_TEST=<suite>`.
function run_worker(cmd, suite)
    run(addenv(cmd, "TRIXI_TEST_RUN_ITEMS" => "true", "TRIXI_TEST" => suite))
end

function dispatch_special_suite(suite)
    project = dirname(Base.active_project())
    julia = Base.julia_cmd()

    if suite == "mpi"
        cmd = `$(MPI.mpiexec()) -n $TRIXI_MPI_NPROCS $julia --threads=1 --check-bounds=yes --project=$project $(@__FILE__)`
        run_worker(cmd, suite)
    elseif suite == "threaded"
        cmd = `$julia --threads=$TRIXI_NTHREADS --check-bounds=yes --code-coverage=none --project=$project $(@__FILE__)`
        run_worker(cmd, suite)
    else
        error("We should not reach this branch; something is wrong.")
    end
end

if !isempty(SUITES_TO_DISPATCH)
    # Dispatch each requested special suite into its own process; the worker(s)
    # run the tagged items. For a single special suite that is all we do; for
    # `all` we additionally run the remaining (in-process) items below.
    foreach(dispatch_special_suite, SUITES_TO_DISPATCH)
end

if !IN_WORKER && TRIXI_TEST in SPECIAL_PROCESS_SUITES
    # A single special suite was requested and dispatched above; nothing to run
    # in this process.
else
    # In-process run. Either we are inside a launched worker (run just that suite's
    # tagged items), or this is an ordinary partition / `all` (in which case we
    # exclude the special suites, which must run in their own processes).
    special_tags = Symbol.(SPECIAL_PROCESS_SUITES)
    tag = Symbol(TRIXI_TEST)

    testitem_filter = ti -> begin
        if TRIXI_TEST == "all"
            # The special suites run in dedicated processes (see above).
            return !any(t -> t in ti.tags, special_tags)
        else
            return tag in ti.tags
        end
    end
    @run_package_tests filter=testitem_filter verbose=TRIXI_TEST_VERBOSE
end

# Common setup shared by all `@testitem`s. Listing `setup=[Setup]` on a test item
# makes the `@test_trixi_include` helper macro (and the packages `test_trixiatmo.jl`
# pulls in) available inside it. `using TrixiAtmo`/`using Test` are already provided
# automatically by the default imports of every `@testitem`/`@testsnippet`.
@testsnippet Setup begin
    include("test_trixiatmo.jl")
end
