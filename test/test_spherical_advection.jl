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
                            6.44302432e-02,
                        ],
                        linf=[
                            9.48946345e-02,
                            9.64807833e-02,
                            1.37450721e-02,
                            0.00000000e+00,
                            4.09322999e-01,
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
                            8.97501204e-03,
                            8.74316738e-03,
                            1.99380928e-03,
                            0.00000000e+00,
                            6.45266612e-02,
                        ],
                        linf=[
                            1.01363241e-01,
                            1.03082705e-01,
                            1.41424723e-02,
                            0.00000000e+00,
                            4.10471741e-01,
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
