module TestSphericalAdvection

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

EXAMPLES_DIR = pkgdir(TrixiAtmo, "examples")

@trixiatmo_testset "elixir_shallowwater_covariant_geostrophic_balance" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_covariant_geostrophic_balance.jl"),
                        l2=[
                            0.30653144639858293,
                            0.00019984467582353295,
                            8.767819502806826e-5
                        ],
                        linf=[
                            1.4786544678290738,
                            0.001375460003351342,
                            0.0007564014737142868
                        ],
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

@trixiatmo_testset "elixir_shallowwater_covariant_rossby_haurwitz_EC" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_covariant_rossby_haurwitz_EC.jl"),
                        l2=[265.9858131601718, 0.1769083525032968, 0.2538332624036849],
                        linf=[579.0744773821989, 0.5054767269383715, 0.5628603103981238],
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

@trixiatmo_testset "elixir_shallowwater_covariant_rossby_haurwitz_EC (LLF surface flux)" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_covariant_rossby_haurwitz_EC.jl"),
                        l2=[265.9818409531758, 0.17644362975403768, 0.2535621643485035],
                        linf=[574.6722628324133, 0.5155374806433987, 0.5497046315783312],
                        tspan=(0.0, 1.0 * SECONDS_PER_DAY),
                        surface_flux=(flux_lax_friedrichs, flux_nonconservative_ec))
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
