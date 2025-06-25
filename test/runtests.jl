using Test

# Trixi.jl runs tests in parallel with CI jobs setting the `TRIXI_TEST` environment
# variable to determine the subset of tests to execute.
#
# We could do the same once we have a lot of tests
const TRIXI_TEST = get(ENV, "TRIXI_TEST", "all")
const TRIXI_MPI_NPROCS = clamp(Sys.CPU_THREADS, 2, 3)
const TRIXI_NTHREADS = clamp(Sys.CPU_THREADS, 2, 3)

@time @testset verbose=true showtiming=true "TrixiAtmo.jl tests" begin
   
    @time if TRIXI_TEST == "all" || TRIXI_TEST == "shallow_water_3d"
        include("test_3d_shallow_water.jl")
    end
end
