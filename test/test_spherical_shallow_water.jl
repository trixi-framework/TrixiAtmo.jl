module TestSphericalShallowWater

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

EXAMPLES_DIR = pkgdir(TrixiAtmo, "examples")

@trixiatmo_testset "Spherical shallow water equations, covariant weak form" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_spherical_shallow_water_covariant.jl"),
                        l2=[
                            0.30653144639857915,
                            0.00019984467582354125,
                            8.767819502807496e-5,
                        ],
                        linf=[
                            1.478654467828619,
                            0.0013754600033516612,
                            0.0007564014737141897,
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

@trixiatmo_testset "Spherical shallow water equations, EC split-covariant form" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_spherical_shallow_water_covariant_ec.jl"),
                        l2=[
                            0.9995330728181165,
                            0.00011159271336455212,
                            9.122690167302845e-5,
                        ],
                        linf=[
                            3.659066044232077,
                            0.0012784722821572864,
                            0.0012784722821592154,
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

end # module
