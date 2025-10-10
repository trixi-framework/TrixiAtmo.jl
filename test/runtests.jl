using Test
using MPI: mpiexec

# Trixi.jl runs tests in parallel with CI jobs setting the `TRIXI_TEST` environment
# variable to determine the subset of tests to execute.
#
# We could do the same once we have a lot of tests
const TRIXI_TEST = get(ENV, "TRIXI_TEST", "all")
const TRIXI_MPI_NPROCS = clamp(Sys.CPU_THREADS, 2, 3)
const TRIXI_NTHREADS = clamp(Sys.CPU_THREADS, 2, 3)

@time @testset verbose=true showtiming=true "TrixiAtmo.jl tests" begin
    # This is placed first since tests error out otherwise if `TRIXI_TEST == "all"`,
    # at least on some systems.
    @time if TRIXI_TEST == "all" || TRIXI_TEST == "mpi"
        # Do a dummy `@test true`:
        # If the process errors out the testset would error out as well,
        # cf. https://github.com/JuliaParallel/MPI.jl/pull/391
        @test true

        status = run(ignorestatus(`$(mpiexec()) -n $TRIXI_MPI_NPROCS $(Base.julia_cmd()) --threads=1 --check-bounds=yes $(abspath("test_mpi.jl"))`))
        @test success(status)
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "threaded"
        # Do a dummy `@test true`:
        # If the process errors out the testset would error out as well,
        # cf. https://github.com/JuliaParallel/MPI.jl/pull/391
        @test true

        status = run(ignorestatus(`$(Base.julia_cmd()) --threads=$TRIXI_NTHREADS --check-bounds=yes --code-coverage=none $(abspath("test_threaded.jl"))`))
        @test success(status)
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "trixi_consistency"
        include("test_trixi_consistency.jl")
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "unit_fluxes"
        include("test_unit.jl")
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "type_stable_tests"
        include("test_type.jl")
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "moist_euler"
        include("test_2d_moist_euler.jl")
        include("test_2d_rainy_euler.jl")
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "spherical_advection"
        include("test_spherical_advection.jl")
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "shallow_water_3d"
        include("test_3d_shallow_water.jl")
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "shallow_water_2d_covariant"
        include("test_2d_shallow_water_covariant.jl")
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "euler_potential_temperature_1d"
        include("test_1d_potential_temperature.jl")
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "euler_potential_temperature_2d"
        include("test_2d_potential_temperature.jl")
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "euler_potential_temperature_3d"
        include("test_3d_potential_temperature.jl")
    end

    @time if TRIXI_TEST == "upstream"
        include("test_trixi_consistency.jl")
    end
end
