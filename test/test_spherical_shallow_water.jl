module TestSphericalShallowWater

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

EXAMPLES_DIR = pkgdir(TrixiAtmo, "examples")

@trixiatmo_testset "Spherical shallow water equations, covariant weak form" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_spherical_shallow_water_covariant.jl"),
                        l2=[
                            3.06531446e-01,
                            1.99844676e-04,
                            8.76781950e-05,
                        ],
                        linf=[
                            1.47865447e+00,
                            1.37546000e-03,
                            7.56401474e-04,
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
                            9.99533073e-01,
                            1.11592713e-04,
                            9.12269017e-05,
                        ],
                        linf=[
                            3.65906604e+00,
                            1.27847228e-03,
                            1.27847228e-03,
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
