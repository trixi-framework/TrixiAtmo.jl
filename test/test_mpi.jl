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
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_moist_euler_moist_bubble.jl"),
                        l2=[
                            7.3515680983123215e-6,
                            1.1067008939664827e-7,
                            0.0006971968385493199,
                            1.715939603224655,
                            8.832720252288771e-7,
                            1.025736269959355e-6
                        ],
                        linf=[
                            8.056395313560394e-5,
                            1.1981461033088162e-6,
                            0.0058959697735631155,
                            19.248694115842227,
                            1.0043632092967755e-5,
                            1.1439573103299433e-5
                        ],
                        polydeg=3,
                        cells_per_dimension=(16, 8),
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
