module TestSphericalAdvection

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

EXAMPLES_DIR = pkgdir(TrixiAtmo, "examples")

@trixiatmo_testset "Spherical shallow water, covariant weak form, global spherical coords" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_covariant_cubed_sphere.jl"),
                        l2=[
                            0.3065314463985936,
                            0.00019984467582352902,
                            8.767819502807507e-5
                        ],
                        linf=[
                            1.4786544678290738,
                            0.0013754600033514114,
                            0.0007564014737144256
                        ], global_coordinate_system=GlobalSphericalCoordinates())
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "Spherical shallow water, covariant weak form, global Cartesian coords" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_covariant_cubed_sphere.jl"),
                        l2=[
                            0.3065314463985936,
                            0.00019984467582352902,
                            8.767819502807507e-5
                        ],
                        linf=[
                            1.4786544678290738,
                            0.0013754600033514114,
                            0.0007564014737144256
                        ], global_coordinate_system=GlobalCartesianCoordinates())
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "Spherical shallow water, entropy-conservative covariant form" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_covariant_cubed_sphere_EC.jl"),
                        l2=[0.9995330728180628,
                            0.000111592713364556,
                            9.12269016730292e-5],
                        linf=[3.659066044233896,
                            0.0012784722821582717,
                            0.0012784722821565994])
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
