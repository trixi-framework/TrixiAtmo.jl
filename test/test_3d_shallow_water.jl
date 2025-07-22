module TestShallowWaterCartesian

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

EXAMPLES_DIR = TrixiAtmo.examples_dir()

@trixiatmo_testset "elixir_shallowwater_cartesian_unsteady_solid_body_rotation_EC_correction" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cartesian_unsteady_solid_body_rotation_EC_correction.jl"),
                        l2=[
                            0.1385840313143962,
                            0.82837506217066,
                            0.1401783834081,
                            0.3170683216476,
                            0.0
                        ],
                        linf=[
                            0.159323444366919,
                            0.7491913049016,
                            0.553600463765,
                            0.0452131916827,
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

@trixiatmo_testset "elixir_shallowwater_cartesian_unsteady_solid_body_rotation_EC_projection" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cartesian_unsteady_solid_body_rotation_EC_projection.jl"),
                        l2=[
                            0.2715065248576713,
                            0.683530367431,
                            0.759403715426,
                            0.96437160416355,
                            0.0
                        ],
                        linf=[
                            0.238740469409095,
                            0.431268687156,
                            0.837234735748,
                            0.6899531778763,
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

@trixiatmo_testset "elixir_shallowwater_cartesian_unsteady_solid_body_rotation_EC_projection (ES)" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cartesian_unsteady_solid_body_rotation_EC_projection.jl"),
                        l2=[
                            0.2744086984644598,
                            0.22833657858405,
                            0.07258247717635,
                            0.92205847355822,
                            0.0
                        ],
                        linf=[
                            0.4332199421835412,
                            0.449038614228,
                            0.6155024602194,
                            0.3580783745856,
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

@trixiatmo_testset "elixir_shallowwater_cartesian_well_balanced" begin
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

@trixiatmo_testset "elixir_shallowwater_cartesian_geostrophic_balance" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cartesian_geostrophic_balance.jl"),
                        l2=[
                            0.0,
                            0.39838614468358,
                            0.39838614468256,
                            0.517273183733906,
                            0.0
                        ],
                        linf=[
                            0.2383681144717684,
                            0.2955303677882,
                            0.2955303680574,
                            0.4494926100049,
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

@trixiatmo_testset "elixir_shallowwater_cartesian_isolated_mountain" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cartesian_isolated_mountain.jl"),
                        l2=[
                            0.189867835225384,
                            0.890929855556,
                            0.784683604144,
                            0.998709859527,
                            0.0
                        ],
                        linf=[
                            0.53215616900434,
                            0.28060001574,
                            0.814315962474,
                            0.28474927765,
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
