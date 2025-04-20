module TestShallowWaterCovariant

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

EXAMPLES_DIR = TrixiAtmo.examples_dir()

@trixiatmo_testset "elixir_shallowwater_covariant_geostrophic_balance" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_covariant_geostrophic_balance.jl"),
                        l2=[
                            0.27823189465180176,
                            0.00021164119121947655,
                            8.588414721470682e-5
                        ],
                        linf=[
                            1.4786544678290738,
                            0.001375460003351342,
                            0.0007564014737142868
                        ],
                        polydeg=3,
                        cells_per_dimension=(5, 5),
                        tspan=(0.0, 1.0 * SECONDS_PER_DAY))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "elixir_shallowwater_covariant_rossby_haurwitz" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_covariant_rossby_haurwitz.jl"),
                        l2=[265.981826097756, 0.17644364627357367, 0.25356217267195796],
                        linf=[574.6725801771408, 0.5155385127558587, 0.549704048104133],
                        polydeg=3,
                        cells_per_dimension=(5, 5),
                        tspan=(0.0, 1.0 * SECONDS_PER_DAY))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "elixir_shallowwater_covariant_isolated_mountain" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_covariant_isolated_mountain.jl"),
                        l2=[13.188835117913722, 0.005698389870463649, 0.007624148100358777],
                        linf=[116.64112009453402, 0.05208844726941367, 0.07855581195821103],
                        polydeg=3,
                        cells_per_dimension=(5, 5),
                        tspan=(0.0, 1.0 * SECONDS_PER_DAY))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "elixir_shallowwater_covariant_unsteady_solid_body_rotation_EC" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_covariant_unsteady_solid_body_rotation_EC.jl"),
                        l2=[
                            0.25498720877863845,
                            0.0002689081664431078,
                            0.00015718255978420539
                        ],
                        linf=[
                            2.135538516899487,
                            0.0013076854064064114,
                            0.0013246467546383922
                        ], # For coverage, test with central flux here instead of usual EC
                        volume_flux=(flux_central, flux_nonconservative_ec),
                        surface_flux=(flux_lax_friedrichs, flux_nonconservative_ec),
                        polydeg=3,
                        cells_per_dimension=(5, 5),
                        tspan=(0.0, 1.0 * SECONDS_PER_DAY))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "elixir_shallowwater_covariant_barotropic_instability" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_covariant_barotropic_instability.jl"),
                        l2=[21.085349167645717, 0.030059547976114206, 0.023419991986079282],
                        linf=[122.97647708776276, 0.17997838634063137, 0.16658465951310109],
                        polydeg=3,
                        cells_per_dimension=(5, 5),
                        tspan=(0.0, 1.0 * SECONDS_PER_DAY))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "elixir_shallowwater_covariant_well_balanced" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_covariant_well_balanced.jl"),
                        l2=[0.0, 0.0, 0.0], linf=[0.0, 0.0, 0.0],
                        polydeg=3,
                        cells_per_dimension=(5, 5),
                        tspan=(0.0, 1.0 * SECONDS_PER_DAY), atol=1e-11)
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
