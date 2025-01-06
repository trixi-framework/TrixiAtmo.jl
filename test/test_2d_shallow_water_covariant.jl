module TestSphericalAdvection

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

EXAMPLES_DIR = pkgdir(TrixiAtmo, "examples")

@trixiatmo_testset "elixir_shallowwater_covariant_geostrophic_balance" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_covariant_geostrophic_balance.jl"),
                        l2=[
                            0.3065314463985936,
                            0.00019984467582352902,
                            8.767819502807507e-5
                        ],
                        linf=[
                            1.4786544678290738,
                            0.0013754600033514114,
                            0.0007564014737144256
                        ], tspan=(0.0, 1.0 * SECONDS_PER_DAY))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "elixir_shallowwater_covariant_rossby_haurwitz_EC" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_covariant_rossby_haurwitz_EC.jl"),
                        l2=[
                            265.98581316017135,
                            0.17690835250329687,
                            0.25383326240368487
                        ],
                        linf=[
                            579.0744773822007,
                            0.5054767269383867,
                            0.5628603103981209
                        ], tspan=(0.0, 1.0 * SECONDS_PER_DAY))
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
