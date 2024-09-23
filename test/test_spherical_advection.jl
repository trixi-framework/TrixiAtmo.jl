module TestSphericalAdvection

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

EXAMPLES_DIR = pkgdir(TrixiAtmo, "examples")

@trixiatmo_testset "elixir_euler_spherical_advection_cartesian" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_spherical_advection_cartesian.jl"),
                        l2=[
                            8.44505073e-03,
                            8.23414117e-03,
                            1.84210648e-03,
                            0.00000000e+00,
                            6.44302430e-02,
                        ],
                        linf=[
                            9.48950488e-02,
                            9.64811952e-02,
                            1.37453400e-02,
                            0.00000000e+00,
                            4.09322999e-01,
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
