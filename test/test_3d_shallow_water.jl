module TestShallowWaterCartesian

include("test_trixiatmo.jl")

@trixi_testset "elixir_shallowwater_cartesian_unsteady_solid_body_rotation_EC_correction" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cartesian_unsteady_solid_body_rotation_EC_correction.jl"),
                        l2=[
                            1.1385840313143962,
                            464.82837506217066,
                            469.1401783834081,
                            311.3170683216476,
                            0.0
                        ],
                        linf=[
                            5.159323444366919,
                            3303.7491913049016,
                            3420.553600463765,
                            3730.0452131916827,
                            0.0
                        ],
                        polydeg=3,
                        cells_per_dimension=(5, 5),
                        tspan=(0.0, 1.0 * SECONDS_PER_DAY))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixi_testset "elixir_shallowwater_cartesian_unsteady_solid_body_rotation_EC_projection" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cartesian_unsteady_solid_body_rotation_EC_projection.jl"),
                        l2=[
                            1.2715065248576713,
                            598.683530367431,
                            605.759403715426,
                            460.96437160416355,
                            0.0
                        ],
                        linf=[
                            4.238740469409095,
                            5466.431268687156,
                            5083.837234735748,
                            3502.6899531778763,
                            0.0
                        ],
                        polydeg=3,
                        cells_per_dimension=(5, 5),
                        tspan=(0.0, 1.0 * SECONDS_PER_DAY))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixi_testset "elixir_shallowwater_cartesian_unsteady_solid_body_rotation_EC_projection (ES)" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cartesian_unsteady_solid_body_rotation_EC_projection.jl"),
                        l2=[
                            0.2744086984644598,
                            280.22833657858405,
                            294.07258247717635,
                            187.92205847355822,
                            0.0
                        ],
                        linf=[
                            1.4332199421835412,
                            1255.449038614228,
                            1470.6155024602194,
                            1249.3580783745856,
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
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixi_testset "elixir_shallowwater_cartesian_well_balanced" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cartesian_well_balanced.jl"),
                        l2=[0.0, 0.0, 0.0, 0.0, 0.0],
                        linf=[0.0, 0.0, 0.0, 0.0, 0.0],
                        atol=8.0e-11) # Needs a slightly larger tolerance for linf
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixi_testset "elixir_shallowwater_cartesian_geostrophic_balance" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cartesian_geostrophic_balance.jl"),
                        l2=[
                            0.27676841776660416,
                            103.39838614468358,
                            103.39838614468256,
                            47.517273183733906,
                            0.0
                        ],
                        linf=[
                            1.2383681144717684,
                            610.2955303677882,
                            610.2955303680574,
                            276.4494926100049,
                            0.0
                        ],
                        polydeg=3,
                        cells_per_dimension=(5, 5),
                        tspan=(0.0, 1.0 * SECONDS_PER_DAY))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixi_testset "elixir_shallowwater_cartesian_isolated_mountain" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cartesian_isolated_mountain.jl"),
                        l2=[
                            13.189867835225384,
                            4656.890929855556,
                            4027.784683604144,
                            6275.998709859527,
                            0.0
                        ],
                        linf=[
                            115.53215616900434,
                            37970.28060001574,
                            42646.814315962474,
                            65362.28474927765,
                            0.0
                        ],
                        polydeg=3,
                        cells_per_dimension=(5, 5),
                        tspan=(0.0, 1.0 * SECONDS_PER_DAY))
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
