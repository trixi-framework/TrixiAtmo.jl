module TestExamples2DRainyEuler

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

EXAMPLES_DIR = joinpath(TrixiAtmo.examples_dir(), "moist_euler")

@trixiatmo_testset "convergence_test_rainy" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "convergence_test",
                                 "convergence_test_rainy.jl"),
                        l2=[
                            2.39895785368954e-6,
                            3.892427386570826e-5,
                            1.5186141589998252e-5,
                            0.007532307120716242,
                            0.02085654062623111,
                            51.988569260688045,
                            1.7470730680079724e-8,
                            3.892376665583286e-5,
                            0.0002156633044466944
                        ],
                        linf=[
                            1.675032502916618e-5,
                            0.0002740479581908595,
                            8.691475272920579e-5,
                            0.023414703432404593,
                            0.06248365429859781,
                            385.10196906188503,
                            6.003313225070965e-8,
                            0.000274026861750043,
                            0.0007433951467703537
                        ])
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "elixir_rainy_euler_moist_bubble" begin
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
                        cells_per_dimension=(64, 32),
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

@trixiatmo_testset "elixir_rainy_euler_rainy_bubble" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "rainy_bubble",
                                 "elixir_rainy_euler_rainy_bubble.jl"),
                        l2=[
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0
                        ],
                        linf=[
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0
                        ],
                        cells_per_dimension = (16, 16),
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
                        initial_refinement_level=4,
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
