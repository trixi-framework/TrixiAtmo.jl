module TestShallowWaterCartesian

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

EXAMPLES_DIR = TrixiAtmo.examples_dir()

@trixiatmo_testset "elixir_shallowwater_cubed_sphere_shell_EC_correction" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cubed_sphere_shell_EC_correction.jl"),
                        l2=[
                            0.07271743465440743,
                            0.07435870820933804,
                            0.05059253372287277,
                            0.07712245696955002,
                            0.002166321001178188
                        ],
                        linf=[
                            1.057593830959098,
                            0.7420535585104321,
                            0.33165510628037276,
                            0.6538582658635804,
                            0.009155117525905546
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

@trixiatmo_testset "elixir_shallowwater_cubed_sphere_shell_EC_projection" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cubed_sphere_shell_EC_projection.jl"),
                        l2=[
                            0.07281376793171955,
                            0.07425222671946745,
                            0.05069528390329348,
                            0.07741804767929857,
                            0.002166321001178188
                        ],
                        linf=[
                            1.0576298065256564,
                            0.7245020984321977,
                            0.3273888680444672,
                            0.6494576573261298,
                            0.009155117525905546
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

@trixiatmo_testset "elixir_shallowwater_cubed_sphere_shell_ES_projection" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cubed_sphere_shell_EC_projection.jl"),
                        l2=[
                            0.06942328774989055,
                            0.06815337757250019,
                            0.038564087277833003,
                            0.07127736553644314,
                            0.002166321001178188
                        ],
                        linf=[
                            1.046573189544382,
                            0.7379677511543662,
                            0.1443340255302651,
                            0.661648435853761,
                            0.009155117525905546
                        ],
                        surface_flux=(FluxPlusDissipation(flux_wintermeyer_etal,
                                                          DissipationLaxFriedrichsEntropyVariables()),
                                      flux_nonconservative_wintermeyer_etal))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "elixir_shallowwater_cubed_sphere_shell_standard" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cubed_sphere_shell_standard.jl"),
                        l2=[
                            0.06866722672214584,
                            0.09439341808798976,
                            0.03259642977310208,
                            0.020263919758088205,
                            0.0
                        ],
                        linf=[
                            1.046030865500346,
                            1.0593845602505878,
                            0.11825902590307225,
                            0.0555545985943406,
                            0.0
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

@trixiatmo_testset "elixir_shallowwater_cartesian_isolated_mountain" begin
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
