module TestShallowWaterCovariant

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

EXAMPLES_DIR = TrixiAtmo.examples_dir()

@trixiatmo_testset "elixir_shallowwater_covariant_geostrophic_balance" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_covariant_geostrophic_balance.jl"),
                        l2=[
                            0.2782318946518096,
                            0.00021164119121947804,
                            8.588414721470694e-5
                        ],
                        linf=[
                            1.818810315677183,
                            0.0011402373824409423,
                            0.0010284231183849968
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
                        l2=[265.98095267805803, 0.1764452313672691, 0.25356292962842114],
                        linf=[574.6708171274768, 0.5155477330239119, 0.5497352934063009],
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
                        l2=[13.18894432799001, 0.005698447961168719, 0.007624217062402512],
                        linf=[116.645494528163, 0.052086295524203324, 0.07855675891709994],
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
                            0.25500246412246835,
                            0.0002703652960705074,
                            0.0001650991084937778
                        ],
                        linf=[
                            2.11951898785469,
                            0.0013109661149978136,
                            0.001328788247668855
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
                        l2=[21.08826693663232, 0.03006187671520436, 0.023421745045307123],
                        linf=[122.9994523425994, 0.17997299389835533, 0.16659612583251238],
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
