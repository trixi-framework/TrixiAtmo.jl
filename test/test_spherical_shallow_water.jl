module TestSphericalShallowWater

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

EXAMPLES_DIR = pkgdir(TrixiAtmo, "examples")

@trixiatmo_testset "Spherical shallow water equations, covariant form" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_spherical_shallow_water_covariant.jl"),
                        l2=[
                            2.14259253e-01,
                            8.16220315e-05,
                            5.05931104e-05,
                        ],
                        linf=[
                            7.42310936e-01,
                            4.02787433e-04,
                            2.68019703e-04,
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
