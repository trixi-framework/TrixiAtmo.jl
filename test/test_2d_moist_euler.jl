module TestExamples2DMoistEuler

using Test
using TrixiAtmo

include("test_trixiatmo.jl") # TODO - This is a repetition from Trixi.jl

EXAMPLES_DIR = pkgdir(TrixiAtmo, "examples") # TODO - Do we need a subdirectory for examples?

@trixiatmo_testset "elixir_moist_euler_dry_bubble.jl" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_moist_euler_dry_bubble.jl"),
                        l2=[
                            1.300428671901329e-6,
                            2.601090012108739e-5,
                            0.0006660124630171347,
                            0.008969786054960861,
                            0.0,
                            0.0,
                        ],
                        linf=[
                            1.0312042909910168e-5,
                            0.00020625488871672815,
                            0.006392107590872236,
                            0.07612038028310053,
                            0.0,
                            0.0,
                        ],
                        polydeg=3,
                        tspan=(0.0, 0.1))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        analysis_callback(sol)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "elixir_moist_euler_EC_bubble" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_moist_euler_EC_bubble.jl"),
                        l2=[
                            0.01345154393018332,
                            0.8070926361417218,
                            7.938812668709457,
                            4500.747616411578,
                            0.00015592413050814787,
                            0.00014163475049532796,
                        ],
                        linf=[
                            0.1427479052298466,
                            8.564879578662357,
                            91.56822550162855,
                            49528.383866247605,
                            0.0019364397602254623,
                            0.0013259689889851285,
                        ],
                        polydeg=3,
                        cells_per_dimension=(16, 16),
                        tspan=(0.0, 0.1))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        analysis_callback(sol)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

end # module
