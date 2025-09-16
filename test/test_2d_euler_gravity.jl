module TestSphericalAdvection

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

EXAMPLES_DIR = pkgdir(TrixiAtmo, "examples")

@trixiatmo_testset "elixir_euler_gravity_equilibrium (Shima + LLF)" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_gravity_equilibrium.jl"),
                        l2=[0.0, 0.0, 0.0, 0.0, 0.0],
                        linf=[0.0, 0.0, 0.0, 0.0, 0.0])
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "elixir_euler_gravity_equilibrium (Kennedy-Gruber)" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_gravity_equilibrium.jl"),
                        l2=[0.0, 0.0, 0.0, 0.0, 0.0],
                        linf=[0.0, 0.0, 0.0, 0.0, 0.0],
                        volume_flux=(flux_kennedy_gruber, flux_nonconservative_waruszewski),
                        surface_flux=(flux_kennedy_gruber, flux_nonconservative_waruszewski))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "elixir_euler_gravity_equilibrium (Ranocha)" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_gravity_equilibrium.jl"),
                        l2=[0.0, 0.0, 0.0, 0.0, 0.0],
                        linf=[0.0, 0.0, 0.0, 0.0, 0.0],
                        volume_flux=(flux_ranocha, flux_nonconservative_waruszewski),
                        surface_flux=(flux_ranocha, flux_nonconservative_waruszewski))
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
