using Test

# Trixi.jl runs tests in parallel with CI jobs setting the `TRIXI_TEST` environment
# variable to determine the subset of tests to execute.
#
# We could do the same once we have a lot of tests
const TRIXIATMO_TEST = get(ENV, "TRIXIATMO_TEST", "all")
const TRIXIATMO_MPI_NPROCS = clamp(Sys.CPU_THREADS, 2, 3)
const TRIXIATMO_NTHREADS = clamp(Sys.CPU_THREADS, 2, 3)

@time @testset verbose=true showtiming=true "TrixiAtmo.jl tests" begin
    @time if TRIXIATMO_TEST == "all" || TRIXIATMO_TEST == "trixi_consistency"
        include("test_trixi_consistency.jl")
    end

    @time if TRIXIATMO_TEST == "all" || TRIXIATMO_TEST == "moist_euler"
        include("test_2d_moist_euler.jl")
    end

    @time if TRIXIATMO_TEST == "all" || TRIXIATMO_TEST == "spherical_advection"
        include("test_spherical_advection.jl")
    end

    @time if TRIXIATMO_TEST == "all" || TRIXIATMO_TEST == "shallow_water_3d"
        include("test_3d_shallow_water.jl")
    end

    @time if TRIXIATMO_TEST == "all" || TRIXIATMO_TEST == "threaded"
        # Do a dummy `@test true`:
        # If the process errors out the testset would error out as well,
        # cf. https://github.com/JuliaParallel/MPI.jl/pull/391
        @test true

        run(`$(Base.julia_cmd()) --threads=$TRIXIATMO_NTHREADS --check-bounds=yes --code-coverage=none $(abspath("test_threaded.jl"))`)
    end
end
