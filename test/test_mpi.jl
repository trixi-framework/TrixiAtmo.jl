module TestMPI

using Trixi: Trixi

include("test_trixiatmo.jl")

# Start with a clean environment: remove Trixi.jl output directory if it exists
outdir = "out"
Trixi.mpi_isroot() && isdir(outdir) && rm(outdir, recursive = true)
Trixi.MPI.Barrier(Trixi.mpi_comm())

@testset verbose=true showtiming=true "MPI tests" begin
#! format: noindent

@trixi_testset "elixir_moist_euler_moist_bubble" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_moist_euler_dry_bubble.jl"),
                        l2=[
                            9.103834949215729e-7,
                            1.8209333828866736e-5,
                            0.0004709417153612775,
                            0.006342004383628925,
                            0.0,
                            0.0
                        ],
                        linf=[
                            1.0258282803210506e-5,
                            0.0002051932980897675,
                            0.006394867661494521,
                            0.076401537633501,
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
