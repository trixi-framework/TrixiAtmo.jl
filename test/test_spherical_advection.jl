module TestSphericalAdvection

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

EXAMPLES_DIR = pkgdir(TrixiAtmo, "examples")

@trixiatmo_testset "Spherical advection, Cartesian form with standard mapping" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_spherical_advection_cartesian.jl"),
                        l2=[
                            8.44498914e-03,
                            8.23407970e-03,
                            1.84210216e-03,
                            0.00000000e+00,
                            6.44302432e-02
                        ],
                        linf=[
                            9.48946345e-02,
                            9.64807833e-02,
                            1.37450721e-02,
                            0.00000000e+00,
                            0.40932295554303133
                        ], element_local_mapping=false)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "Spherical advection, Cartesian form with element-local mapping" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_spherical_advection_cartesian.jl"),
                        l2=[
                            0.00897499647958351,
                            0.008743150662795171,
                            0.001993813300529814,
                            0.0,
                            0.06452666141326718
                        ],
                        linf=[
                            0.10136306017376362,
                            0.10308252591097666,
                            0.014142319489653277,
                            0.0,
                            0.41047171697830986
                        ], element_local_mapping=true)
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
