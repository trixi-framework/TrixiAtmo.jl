module TestExamples1DEulerPotentialTemperature

include("test_trixiatmo.jl")

@trixi_testset "elixir_euler_potential_temperature_ec" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_ec.jl"),
                        l2=[1.536920933776253, 1.536920944343481, 8.204456414321612e-7],
                        linf=[2.283573460625467, 2.2835683135782463, 4.835357438226495e-6],
                        tspan=(0.0, 0.4), atol=3.8e-5)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_euler_potential_temperature_tec" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_ec.jl"),
                        l2=[1.5369209174209315, 1.5369209181196604, 1.9668131414528186e-7],
                        linf=[2.2835833418413567, 2.283588054059723, 1.1052571658176635e-6],
                        tspan=(0.0, 0.4), surface_flux=flux_tec, volume_flux=flux_tec,
                        atol=2e-5)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_euler_potential_temperature_etec" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_ec.jl"),
                        l2=[1.5369209174090952, 1.5369209159246773, 4.6312969205197763e-7],
                        linf=[2.2835901549809057, 2.2836015746559717, 2.547940256003578e-6],
                        tspan=(0.0, 0.4), surface_flux=flux_etec, volume_flux=flux_etec,
                        atol=1.5e-5)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_euler_potential_temperature_well_balanced_1d" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_well_balanced.jl"),
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
                        ],)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_euler_potential_temperature_well_balanced_1d" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_well_balanced.jl"),
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
                        surface_flux=(flux_ec, flux_nonconservative_waruzewski_etal),
                        volume_flux=(flux_etec, flux_nonconservative_souza_etal))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

end
