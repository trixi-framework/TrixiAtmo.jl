module TestSphericalAdvection

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

TEST_POLYDEG = 3
TEST_CELLS_PER_DIMENSION = (5, 5)
TEST_TSPAN = (0.0, 1.0 * SECONDS_PER_DAY)

EXAMPLES_DIR = pkgdir(TrixiAtmo, "examples")

@trixiatmo_testset "elixir_shallowwater_covariant_geostrophic_balance" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_covariant_geostrophic_balance.jl"),
                        l2=[
                            0.30653144639858293,
                            0.00019984467582353295,
                            8.767819502806826e-5
                        ],
                        linf=[
                            1.4786544678290738,
                            0.001375460003351342,
                            0.0007564014737142868
                        ],
                        polydeg=TEST_POLYDEG,
                        cells_per_dimension=TEST_CELLS_PER_DIMENSION,
                        tspan=TEST_TSPAN)
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
                        polydeg=TEST_POLYDEG,
                        cells_per_dimension=TEST_CELLS_PER_DIMENSION,
                        tspan=TEST_TSPAN)
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
                        polydeg=TEST_POLYDEG,
                        cells_per_dimension=TEST_CELLS_PER_DIMENSION,
                        tspan=TEST_TSPAN)
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
                            0.015841167275942963,
                            3.470597408352903e-5,
                            1.4875776786557202e-5
                        ],
                        linf=[
                            0.1176405542039447,
                            0.0001442044402930609,
                            0.00013001377736920894
                        ], # For coverage, test with central flux here instead of usual EC
                        volume_flux=(flux_central, flux_nonconservative_ec),
                        surface_flux=(flux_lax_friedrichs, flux_nonconservative_ec),
                        polydeg=TEST_POLYDEG,
                        cells_per_dimension=TEST_CELLS_PER_DIMENSION,
                        tspan=TEST_TSPAN)
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
