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
end

# Clean up afterwards: delete Trixi.jl output directory
Trixi.mpi_isroot() && @test_nowarn rm(outdir, recursive = true)
Trixi.MPI.Barrier(Trixi.mpi_comm())

end # module
