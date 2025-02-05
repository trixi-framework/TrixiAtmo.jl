module TestSphericalAdvection

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

EXAMPLES_DIR = pkgdir(TrixiAtmo, "examples")

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
                            0.06866722672214583,
                            0.09439341808798994,
                            0.03259642977310212,
                            0.02026391975808816,
                            0.0
                        ],
                        linf=[
                            1.046030865500347,
                            1.0593845602505891,
                            0.11825902590307091,
                            0.055554598594340664,
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

@trixiatmo_testset "elixir_shallowwater_cubed_sphere_shell_well_balancing" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cubed_sphere_shell_well_balancing.jl"),
                        l2=[
                            0.00000000e+00,
                            0.00000000e+00,
                            0.00000000e+00,
                            0.00000000e+00,
                            0.00000000e+00
                        ],
                        linf=[
                            0.00000000e+00,
                            0.00000000e+00,
                            0.00000000e+00,
                            0.00000000e+00,
                            0.00000000e+00
                        ], atol=8.0e-13) # Needs a slightly larger tolerance for linf
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
