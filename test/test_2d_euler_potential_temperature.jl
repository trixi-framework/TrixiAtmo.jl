module TestExamples2DEulerPotentialTemperature

using Test
using TrixiAtmo

include("test_trixiatmo.jl") # TODO - This is a repetition from Trixi.jl

EXAMPLES_DIR = pkgdir(TrixiAtmo, "examples") # TODO - Do we need a subdirectory for examples?

@trixiatmo_testset "elixir_euler_potential_temperature_dry_bubble" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_dry_bubble.jl"),
                        l2=[
                            1.300427456987984e-6,
                            2.6010873806861595e-5,
                            0.0006660120005093007,
                            9.51074191163579e-6,
                        ],
                        linf=[
                            1.031131432183141e-5,
                            0.00020623855042956052,
                            0.006392070867001616,
                            7.841424786647622e-5,
                        ],
                        polydeg=3,
                        tspan=(0.0, 0.1))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "elixir_euler_potential_temperature_gravity_wave" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_gravity_wave.jl"),
                        l2=[
                            8.92434879991671e-9,
                            1.8131676850378482e-7,
                            2.6374650049135436e-5,
                            1.4620631924879524e-6,
                        ],
                        linf=[
                            7.544232794032268e-8,
                            1.654911628179434e-6,
                            0.00023724858495745153,
                            1.8884994915424613e-5,
                        ],
                        polydeg=3,
                        tspan=(0.0, 1.0))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "elixir_euler_potential_temperature_ec" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_ec.jl"),
                        l2=[
                            0.06174365254988615,
                            0.05018008812040643,
                            0.05018774429827882,
                            0.0057820937341387935,
                        ],
                        linf=[
                            0.2932352942992503,
                            0.3108954213591686,
                            0.31082193857647905,
                            0.027465217019157606,
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
