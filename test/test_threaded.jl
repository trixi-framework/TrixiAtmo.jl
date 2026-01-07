module TestThreaded

include("test_trixiatmo.jl")

@testset verbose=true showtiming=true "Threaded tests" begin
#! format: noindent

@trixi_testset "elixir_gemein_bubble moist" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "cartesian_euler/moist_air/buoyancy",
                                 "elixir_gemein_bubble.jl"),
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
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 2000)
end
end
end # module
