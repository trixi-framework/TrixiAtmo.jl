module TestShallowWaterCovariant

include("test_trixiatmo.jl")

EXAMPLES_DIR = joinpath(EXAMPLES_DIR, "shallowwater/covariant")

@trixi_testset "elixir_geostrophic_balance" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_geostrophic_balance.jl"),
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
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_rossby_haurwitz" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_rossby_haurwitz.jl"),
                        l2=[265.9818260977567, 0.17644364627357362, 0.2535621726719579],
                        linf=[574.6725801771354, 0.5155385127558593, 0.5497040481041348],
                        polydeg=3,
                        cells_per_dimension=(5, 5),
                        tspan=(0.0, 1.0 * SECONDS_PER_DAY),
                        metric_terms=MetricTermsCovariantSphere(christoffel_symbols = ChristoffelSymbolsCollocationDerivative()))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_isolated_mountain" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_isolated_mountain.jl"),
                        l2=[13.18894432799001, 0.005698447961168719, 0.007624217062402512],
                        linf=[116.645494528163, 0.052086295524203324, 0.07855675891709994],
                        polydeg=3,
                        cells_per_dimension=(5, 5),
                        tspan=(0.0, 1.0 * SECONDS_PER_DAY))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_unsteady_solid_body_rotation_EC" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_unsteady_solid_body_rotation_EC.jl"),
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
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_barotropic_instability" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_barotropic_instability.jl"),
                        l2=[21.08826693663232, 0.03006187671520436, 0.023421745045307123],
                        linf=[122.9994523425994, 0.17997299389835533, 0.16659612583251238],
                        polydeg=3,
                        cells_per_dimension=(5, 5),
                        tspan=(0.0, 1.0 * SECONDS_PER_DAY))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_well_balanced" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_well_balanced.jl"),
                        l2=[0.0, 0.0, 0.0], linf=[0.0, 0.0, 0.0],
                        polydeg=3,
                        cells_per_dimension=(5, 5),
                        tspan=(0.0, 1.0 * SECONDS_PER_DAY), atol=1e-11)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

end # module
