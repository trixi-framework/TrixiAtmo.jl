module TestShallowWaterCartesian

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

EXAMPLES_DIR = pkgdir(TrixiAtmo, "examples")

@trixiatmo_testset "elixir_shallowwater_cubed_sphere_shell_EC_correction" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cubed_sphere_shell_EC_correction.jl"),
                        l2=[
                            7.23800458e-02,
                            9.98871590e-02,
                            4.55606969e-02,
                            3.17422064e-02,
                            0.00000000e+00
                        ],
                        linf=[
                            1.05686060e+00,
                            1.04413842e+00,
                            3.12374228e-01,
                            3.19064636e-01,
                            0.00000000e+00
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

@trixiatmo_testset "elixir_shallowwater_cubed_sphere_shell_EC_projection" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cubed_sphere_shell_EC_projection.jl"),
                        l2=[
                            7.25131271e-02,
                            1.01017589e-01,
                            4.39697415e-02,
                            3.08898940e-02,
                            0.00000000e+00
                        ],
                        linf=[
                            1.06007536e+00,
                            1.05786719e+00,
                            3.93826726e-01,
                            2.34129714e-01,
                            0.00000000e+00
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

@trixiatmo_testset "elixir_shallowwater_cubed_sphere_shell_standard" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cubed_sphere_shell_standard.jl"),
                        l2=[
                            6.86661496e-02,
                            9.43310712e-02,
                            3.26202302e-02,
                            2.02514293e-02,
                            0.00000000e+00
                        ],
                        linf=[
                            1.04610420e+00,
                            1.05891219e+00,
                            1.19809349e-01,
                            5.54354032e-02,
                            0.00000000e+00
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
