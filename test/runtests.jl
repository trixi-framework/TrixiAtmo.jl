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
        @testset verbose=true showtiming=true "MPI tests" begin
            # Do a dummy `@test true`:
            # If the process errors out the testset would error out as well,
            # cf. https://github.com/JuliaParallel/MPI.jl/pull/391
            @test true

            status = run(ignorestatus(`$(mpiexec()) -n $TRIXI_MPI_NPROCS $(Base.julia_cmd()) --threads=1 --check-bounds=yes $(abspath("test_mpi.jl"))`))
            @test success(status)
        end
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "threaded"
        @testset verbose=true showtiming=true "Threaded tests" begin
            # Do a dummy `@test true`:
            # If the process errors out the testset would error out as well,
            # cf. https://github.com/JuliaParallel/MPI.jl/pull/391
            @test true

            status = run(ignorestatus(`$(Base.julia_cmd()) --threads=$TRIXI_NTHREADS --check-bounds=yes --code-coverage=none $(abspath("test_threaded.jl"))`))
            @test success(status)
        end
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "trixi_consistency" ||
             TRIXI_TEST == "upstream"
        @testset verbose=true showtiming=true "Trixi.jl consistency tests" begin
            include("test_trixi_consistency.jl")
        end
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "unit_fluxes"
        @testset verbose=true showtiming=true "Unit tests" begin
            include("test_unit.jl")
        end
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "type_stable_tests"
        @testset verbose=true showtiming=true "Type stability tests" begin
            include("test_type.jl")
        end
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "moist_euler"
        # Thesis Gemein
        @testset verbose=true showtiming=true "CompressibleMoistEulerEquations2D tests" begin
            include("test_2d_moist_euler.jl")
        end
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "rainy_euler"
        # Thesis Höck
        @testset verbose=true showtiming=true "CompressibleRainyEulerEquations2D tests" begin
            include("test_2d_rainy_euler.jl")
        end
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "spherical_advection"
        @testset verbose=true showtiming=true "Spherical advection tests" begin
            include("test_spherical_advection.jl")
        end
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "shallow_water_3d"
        @testset verbose=true showtiming=true "Spherical SWE Cartesian tests" begin
            include("test_3d_shallow_water.jl")
        end
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "shallow_water_2d_covariant"
        @testset verbose=true showtiming=true "Spherical SWE covariant tests" begin
            include("test_2d_shallow_water_covariant.jl")
        end
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "euler_potential_temperature_1d"
        @testset verbose=true showtiming=true "Euler potential temperatur 1D tests" begin
            include("test_1d_potential_temperature.jl")
        end
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "euler_potential_temperature_2d"
        @testset verbose=true showtiming=true "Euler potential temperatur 2D tests" begin
            include("test_2d_potential_temperature.jl")
        end
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "euler_potential_temperature_3d"
        @testset verbose=true showtiming=true "Euler potential temperatur 3D tests" begin
            include("test_3d_potential_temperature.jl")
        end
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "euler_energy_2d"
        @testset verbose=true showtiming=true "Euler internal energy 2D tests" begin
            include("test_2d_euler_energy.jl")
        end
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "euler_energy_3d"
        @testset verbose=true showtiming=true "Euler internal energy 3D tests" begin
            include("test_3d_euler_energy.jl")
        end
    end
end
