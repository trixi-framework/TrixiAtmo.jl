@testsnippet EulerPotentialTemperature1D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "euler", "dry_air", "tests")
end

@testitem "Euler potential temperature 1D: elixir_potential_temperature_well_balanced_1d" setup=[
    Setup,
    EulerPotentialTemperature1D
] tags=[:euler_potential_temperature_1d] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_potential_temperature_well_balanced.jl"),
                        l2=[
                            8.475270558027727e-10,
                            2.6389816902456284e-13,
                            2.542581117596382e-7,
                            5.776677968116953e-13
                        ],
                        linf=[
                            1.157208107116503e-9,
                            6.721044580141687e-13,
                            3.471623699624615e-7,
                            1.8189894035458565e-12
                        ], atol=1e-5)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@testitem "Euler potential temperature 1D: elixir_potential_temperature_well_balanced_1d with isothermal/Souza fluxes" setup=[
    Setup,
    EulerPotentialTemperature1D
] tags=[:euler_potential_temperature_1d] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_potential_temperature_well_balanced.jl"),
                        l2=[
                            4.16546935255369e-9,
                            1.2863411480535975e-6,
                            1.2496407677658893e-6,
                            5.776677968116953e-13
                        ],
                        linf=[
                            1.6270381930638678e-8,
                            6.254416734833885e-6,
                            4.881114818999777e-6,
                            1.8189894035458565e-12
                        ],
                        surface_flux=(flux_ec, flux_nonconservative_isothermal),
                        volume_flux=(flux_etec, flux_nonconservative_souza_etal), atol=1e-5)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end
