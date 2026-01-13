module TestExamples2DMoistEuler

include("test_trixiatmo.jl")

EXAMPLES_DIR = joinpath(EXAMPLES_DIR, "euler")

@trixi_testset "elixir_gemein_bubble dry" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "dry_air/buoyancy",
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
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_gemein_bubble moist" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "moist_air/buoyancy",
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
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_gemein_nonhydrostatic_gravity_waves" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "dry_air/buoyancy",
                                 "elixir_gemein_nonhydrostatic_gravity_waves.jl"),
                        l2=[
                            3.54205348345642e-5,
                            0.00212660009365362,
                            0.011728423281347147,
                            9.898615652683135,
                            0.0,
                            0.0
                        ],
                        linf=[
                            0.0017602402189631494,
                            0.1494206825191764,
                            0.5945421207691134,
                            489.9007739148801,
                            0.0,
                            0.0
                        ],
                        polydeg=3,
                        cells_per_dimension=(10, 8),
                        tspan=(0.0, 0.1))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_gemein_source_terms dry" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "dry_air/tests",
                                 "elixir_gemein_source_terms.jl"),
                        l2=[
                            1.3992076791281227e-5,
                            1.4486417486907815e-5,
                            2.609465609739115e-5,
                            6.323484379066432e-5,
                            0.0,
                            0.0
                        ],
                        linf=[
                            7.549984224430872e-5,
                            0.00010065352517929504,
                            0.00015964938767742964,
                            0.0005425860570698049,
                            0.0,
                            0.0
                        ],
                        polydeg=3,
                        cells_per_dimension=(10, 8),
                        tspan=(0.0, 0.1))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_gemein_source_terms moist" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "moist_air/tests",
                                 "elixir_gemein_source_terms.jl"),
                        l2=[
                            0.0001480393714768623,
                            0.11945481031503036,
                            0.07460345535073129,
                            5.9430508465000385,
                            4.471792794168725e-9,
                            7.10320253652373e-10
                        ],
                        linf=[
                            0.0007084175808387272,
                            0.504796333566,
                            0.3697160820183014,
                            27.843165651072923,
                            2.1168438904322837e-8,
                            3.691699932047233e-9
                        ],
                        polydeg=3,
                        cells_per_dimension=(10, 8),
                        tspan=(0.0, 0.1))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

end # module
