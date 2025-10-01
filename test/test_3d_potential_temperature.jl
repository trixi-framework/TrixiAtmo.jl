module TestExamples3DEulerPotentialTemperature

include("test_trixiatmo.jl")

@trixi_testset "elixir_euler_potential_temperature_taylor_green_vortex_3d" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_taylor_green_vortex.jl"),
                        l2=[
                            0.007135833158819355,
                            0.05920026666636394,
                            0.058302519710783984,
                            0.09239792853460427,
                            0.0026512482044306874
                        ],
                        linf=[
                            0.025211900914064778,
                            0.1626379859810771,
                            0.1567281830166123,
                            0.2914137341434033,
                            0.010285071822802083
                        ],
                        tspan=(0.0, 1.0))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_euler_potential_temperature_taylor_green_vortex_3d" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_taylor_green_vortex.jl"),
                        l2=[
                            0.00713580882298077,
                            0.05920027226189781,
                            0.05830253041920023,
                            0.09239789909513721,
                            0.0026512387869924467
                        ],
                        linf=[
                            0.025211360864895616,
                            0.1626378819726536,
                            0.15673002075139453,
                            0.2914131870829597,
                            0.010284993042893986
                        ],
                        tspan=(0.0, 1.0), surface_flux=flux_ec, volume_flux=flux_ec)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_euler_potential_temperature_taylor_green_vortex_3d" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_taylor_green_vortex.jl"),
                        l2=[
                            0.007135944983851745,
                            0.05920027807981362,
                            0.05830252313445877,
                            0.09239791128377808,
                            0.0026512477402571037
                        ],
                        linf=[
                            0.02521600221697473,
                            0.16264074629051783,
                            0.15673138446141865,
                            0.29141401520767773,
                            00.010285068999815183
                        ],
                        tspan=(0.0, 1.0), surface_flux=flux_tec, volume_flux=flux_tec)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_euler_potential_temperature_taylor_green_vortex_3d" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_taylor_green_vortex.jl"),
                        l2=[
                            0.006704842021819753,
                            0.018574880297682488,
                            0.018606626958508878,
                            0.010451109496991842,
                            0.0032071378388946444
                        ],
                        linf=[
                            0.025690924689206418,
                            0.04766269256457434,
                            0.050155417906555394,
                            0.030904595875059383,
                            0.012502671674568866
                        ],
                        tspan=(0.0, 0.2), surface_flux=FluxLMARS(340.0),
                        volume_flux=flux_tec, atol=5e-6)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_euler_potential_temperature_baroclinic_instability" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_baroclinic_instability.jl"),
                        l2=[
                            0.0005404243857532021,
                            0.015509341701453983,
                            0.013816888681309869,
                            0.012939539158878471,
                            0.14000483720269236,
                            99.03486863938933
                        ],
                        linf=[
                            0.004654666656827899,
                            1.1039456510593342,
                            0.6663687353552903,
                            1.1276247531737227,
                            1.0016885503585513,
                            333.2487183585763
                        ],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2),
                        atol=3e-9)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_euler_potential_temperature_baroclinic_instability" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_baroclinic_instability.jl"),
                        l2=[
                            0.0005197413655079329,
                            0.01418417025867247,
                            0.012316380118663012,
                            0.012923996341032916,
                            0.13995018762188344,
                            99.03486863938933
                        ],
                        linf=[
                            0.004231440709953516,
                            1.1219699347067773,
                            0.639464645826298,
                            1.1308558982893109,
                            1.0018943890534047,
                            333.2487183585763
                        ],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2),
                        volume_flux=(flux_ec, flux_nonconservative_waruzewski_etal))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_euler_potential_temperature_baroclinic_instability" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_baroclinic_instability.jl"),
                        l2=[
                            0.0005238485705127624,
                            0.014658738174940492,
                            0.012858043298446848,
                            0.012928235735411159,
                            0.13996694572679722,
                            99.03486863938933
                        ],
                        linf=[
                            0.004242590990033657,
                            1.1140990093072491,
                            0.6517474500047321,
                            1.1292857239244998,
                            1.0019906488095103,
                            333.2487183585763
                        ],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2),
                        volume_flux=(flux_etec, flux_nonconservative_artiano_etal))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

end
