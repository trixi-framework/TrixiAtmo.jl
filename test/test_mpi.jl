module TestMPI

using Trixi: Trixi

include("test_trixiatmo.jl")

# Start with a clean environment: remove Trixi.jl output directory if it exists
outdir = "out"
Trixi.mpi_isroot() && isdir(outdir) && rm(outdir, recursive = true)
Trixi.MPI.Barrier(Trixi.mpi_comm())

@testset verbose=true showtiming=true "MPI tests" begin
#! format: noindent
@trixi_testset "elixir_gemein_bubble" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "euler/dry_air/buoyancy",
                                 "elixir_gemein_bubble.jl"),
                        l2=[
                            9.104437114458848e-7,
                            1.8210536975490044e-5,
                            0.0004707887343135412,
                            0.0063400898518523935,
                            0.0,
                            0.0
                        ],
                        linf=[
                            1.0258941581242631e-5,
                            0.00020520634691933992,
                            0.006392782691233334,
                            0.07637640493339859,
                            0.0,
                            0.0
                        ],
                        tspan=(0.0, 0.1))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 1000)
end

@trixi_testset "elixir_potential_temperature_vortex_shedding with Sleve" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "euler", "dry_air", "global_circulation",
                                 "elixir_potential_temperature_vortex_shedding.jl"),
                        l2=[
                            0.0001056726895239208,
                            0.03612902927134097,
                            0.036128771203962594,
                            0.047764018298670316,
                            0.03181247660774069,
                            0.5931596702203235
                        ],
                        linf=[
                            0.0011918734879285964,
                            0.15502631302821748,
                            0.15502631302823136,
                            0.25776642079474976,
                            0.3442945234858712,
                            18.51483958443782
                        ],
                        rtol=1e-9,
                        tspan=(0.0, 0.0001 * SECONDS_PER_DAY),
                        trees_per_cube_face=(3, 2), adapt_vertical_grid=Sleve(0.7, 0.8))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 1000)
    # Check partitioning (a total of 108 elements split into 4 partitions)
    local_nelems = nelements(solver, semi.cache)

    # Perform the parallel reductions and assign the unwrapped results back!
    nelems_min = Trixi.MPI.Allreduce!(Ref(local_nelems), Base.min, Trixi.mpi_comm())[]
    nelems_max = Trixi.MPI.Allreduce!(Ref(local_nelems), Base.max, Trixi.mpi_comm())[]

    @assert nelems_min == 26
    @assert nelems_max == 28
end
end

# Clean up afterwards: delete Trixi.jl output directory
Trixi.mpi_isroot() && @test_nowarn rm(outdir, recursive = true)
Trixi.MPI.Barrier(Trixi.mpi_comm())

end # module
