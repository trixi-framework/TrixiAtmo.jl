module TestExamples2DRainyEuler

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

EXAMPLES_DIR = joinpath(TrixiAtmo.examples_dir(), "moist_euler")

@trixiatmo_testset "elixir_rainy_euler_rainy_bubble_diffusion" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "moist_bubble",
                                 "elixir_rainy_euler_moist_bubble.jl"),
                        l2=[
                            1.300428671901329e-6,
                            2.601090012108739e-5,
                            0.0006660124630171347,
                            0.008969786054960861,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0
                        ],
                        linf=[
                            1.0312042909910168e-5,
                            0.00020625488871672815,
                            0.006392107590872236,
                            0.07612038028310053,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0
                        ],
                        cells_per_dimension = (64, 32),
                        tspan=(0.0, 10.0))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "elixir_rainy_euler_rainy_bubble_diffusion" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "rainy_bubble",
                                 "elixir_rainy_euler_rainy_bubble_diffusion.jl"),
                        l2=[
                            1.300428671901329e-6,
                            2.601090012108739e-5,
                            0.0006660124630171347,
                            0.008969786054960861,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0
                        ],
                        linf=[
                            1.0312042909910168e-5,
                            0.00020625488871672815,
                            0.006392107590872236,
                            0.07612038028310053,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0
                        ],
                        initial_refinement_level=3,
                        tspan=(0.0, 10.0))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end
end # module
