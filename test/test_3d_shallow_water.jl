module TestShallowWaterCartesian

include("test_trixiatmo.jl")

EXAMPLES_DIR = joinpath(EXAMPLES_DIR, "shallow_water/cartesian_manifold")

@trixi_testset "elixir_unsteady_solid_body_rotation_EC_correction" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_unsteady_solid_body_rotation_EC_correction.jl"),
                        l2=[
                            1.1385840313142226,
                            464.8283750621118,
                            469.14017838339083,
                            311.31706832161564,
                            0.0
                        ],
                        linf=[
                            5.159323444358506,
                            3303.749191315932,
                            3420.5536004616565,
                            3730.0452131952625,
                            0.0
                        ],
                        polydeg=3,
                        cells_per_dimension=(5, 5),
                        tspan=(0.0, 1.0 * SECONDS_PER_DAY))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_unsteady_solid_body_rotation_EC_projection" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_unsteady_solid_body_rotation_EC_projection.jl"),
                        l2=[
                            1.271506524857498,
                            598.6835303675092,
                            605.7594037155094,
                            460.9643716042415,
                            0.0
                        ],
                        linf=[
                            4.23874046947094,
                            5466.431268695043,
                            5083.837234738506,
                            3502.6899531773233,
                            0.0
                        ],
                        polydeg=3,
                        cells_per_dimension=(5, 5),
                        tspan=(0.0, 1.0 * SECONDS_PER_DAY))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_unsteady_solid_body_rotation_EC_projection (ES)" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_unsteady_solid_body_rotation_EC_projection.jl"),
                        l2=[
                            0.27440876588211627,
                            280.22773491124406,
                            294.071829622588,
                            187.92193710627467,
                            0.0
                        ],
                        linf=[
                            1.4332269086135057,
                            1255.4454482832807,
                            1470.615003655199,
                            1249.359787903988,
                            0.0
                        ],
                        surface_flux=(FluxPlusDissipation(flux_wintermeyer_etal,
                                                          DissipationLaxFriedrichsEntropyVariables()),
                                      flux_nonconservative_wintermeyer_etal),
                        polydeg=3,
                        cells_per_dimension=(5, 5),
                        tspan=(0.0, 1.0 * SECONDS_PER_DAY))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_well_balanced" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_well_balanced.jl"),
                        l2=[0.0, 0.0, 0.0, 0.0, 0.0],
                        linf=[0.0, 0.0, 0.0, 0.0, 0.0],
                        atol=8.0e-11) # Needs a slightly larger tolerance for linf
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_geostrophic_balance (naive)" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_geostrophic_balance.jl"),
                        l2=[0.27676841776660904,
                            103.39838614468599,
                            103.39838614468121,
                            47.51727318373426, 0.0],
                        linf=[1.238368114471541,
                            610.2955303677882,
                            610.2955303679337,
                            276.44949261002847,
                            0.0],
                        polydeg=3,
                        cells_per_dimension=(5, 5),
                        tspan=(0.0, 1.0 * SECONDS_PER_DAY),
                        surface_flux=(FluxPlusDissipation(flux_wintermeyer_etal,
                                                          DissipationLocalLaxFriedrichs(max_abs_speed_naive)),
                                      flux_nonconservative_wintermeyer_etal)) # use "naive" wave speed estimate for coverage
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_isolated_mountain" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_isolated_mountain.jl"),
                        l2=[
                            13.189868962884406,
                            4656.890871865292,
                            4027.7846474475473,
                            6275.998570998393,
                            0.0
                        ],
                        linf=[
                            115.53214502067749,
                            37970.29034857702,
                            42646.8517588789,
                            65362.34875198507,
                            0.0
                        ],
                        polydeg=3,
                        cells_per_dimension=(5, 5),
                        tspan=(0.0, 1.0 * SECONDS_PER_DAY))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

end # module
