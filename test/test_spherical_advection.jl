module TestSphericalAdvection

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

EXAMPLES_DIR = pkgdir(TrixiAtmo, "examples")

@trixiatmo_testset "Spherical advection, Cartesian weak form" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_spherical_advection_cartesian.jl"),
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

@trixiatmo_testset "Spherical advection, covariant weak form" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_spherical_advection_covariant.jl"),
                        l2=[
                            1.00016205e+00,
                            0.00000000e+00,
                            0.00000000e+00,
                        ],
                        linf=[
                            1.42451839e+01,
                            0.00000000e+00,
                            0.00000000e+00,
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

# The covariant flux-differencing form should be equivalent to the weak form when the 
# arithmetic mean is used as the two-point flux
@trixiatmo_testset "Spherical advection, covariant flux-differencing, arithmetic mean" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_spherical_advection_covariant.jl"),
                        l2=[
                            1.00016205e+00,
                            0.00000000e+00,
                            0.00000000e+00,
                        ],
                        linf=[
                            1.42451839e+01,
                            0.00000000e+00,
                            0.00000000e+00,
                        ], volume_integral=VolumeIntegralFluxDifferencing(flux_central))
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
