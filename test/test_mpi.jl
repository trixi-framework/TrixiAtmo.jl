module TestMPI

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

# Start with a clean environment: remove Trixi.jl output directory if it exists
outdir = "out"
Trixi.mpi_isroot() && isdir(outdir) && rm(outdir, recursive = true)
Trixi.MPI.Barrier(Trixi.mpi_comm())

# CI with MPI and some tests fails often on Windows. Thus, we check whether this
# is the case here. We use GitHub Actions, so we can check whether we run CI
# in the cloud with Windows as follows, see also
# https://docs.github.com/en/actions/learn-github-actions/environment-variables
CI_ON_WINDOWS = (get(ENV, "GITHUB_ACTIONS", false) == "true") && Sys.iswindows()

EXAMPLES_DIR = pkgdir(TrixiAtmo, "examples")

@testset "MPI tests" begin
#! format: noindent

@trixiatmo_testset "elixir_moist_euler_moist_bubble" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_moist_euler_dry_bubble.jl"),
                        l2=[
                            7.351043427240923e-6,
                            1.1070342432069074e-7,
                            0.0006974588377288118,
                            1.715668353329522,
                            8.831269143134121e-7,
                            1.025579538944668e-6
                        ],
                        linf=[
                            8.055695643149896e-5,
                            1.1985203677080201e-6,
                            0.005897639251024437,
                            19.24776030163048,
                            1.0043133039065386e-5,
                            1.1439046776775402e-5
                        ],
                        tspan=(0.0, 0.1))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 2000
    end
end
end

# Clean up afterwards: delete Trixi.jl output directory
Trixi.mpi_isroot() && @test_nowarn rm(outdir, recursive = true)
Trixi.MPI.Barrier(Trixi.mpi_comm())

end # module
