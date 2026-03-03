module TestExamples3DEulerPotentialTemperature

include("test_trixiatmo.jl")

EXAMPLES_DIR = joinpath(EXAMPLES_DIR, "euler/dry_air")

@trixi_testset "elixir_potential_temperature_taylor_green_vortex" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "tests",
                                 "elixir_potential_temperature_taylor_green_vortex.jl"),
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

@trixi_testset "elixir_potential_temperature_taylor_green_vortex" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "tests",
                                 "elixir_potential_temperature_taylor_green_vortex.jl"),
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

@trixi_testset "elixir_potential_temperature_taylor_green_vortex" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "tests",
                                 "elixir_potential_temperature_taylor_green_vortex.jl"),
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

@trixi_testset "elixir_potential_temperature_taylor_green_vortex" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "tests",
                                 "elixir_potential_temperature_taylor_green_vortex.jl"),
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

@trixi_testset "elixir_potential_temperature_baroclinic_instability Souza" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "global_circulation",
                                 "elixir_potential_temperature_baroclinic_instability.jl"),
                        l2=[
                            0.0005402395829626791,
                            0.0155041989645114,
                            0.0138107622096068,
                            0.012941059816206165,
                            0.13997998628314068,
                            98.99596252593278
                        ],
                        linf=[
                            0.00465163790895784,
                            1.103975340100896,
                            0.6663546984155281,
                            1.1276305037303587,
                            1.001311664438731,
                            333.1178018152714
                        ],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2),
                        atol=4e-9)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_potential_temperature_baroclinic_instability Waruszewski" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "global_circulation",
                                 "elixir_potential_temperature_baroclinic_instability.jl"),
                        l2=[
                            0.0005196272070117161,
                            0.0141817969707941,
                            0.012313242513833094,
                            0.012925522828705729,
                            0.1399254026439438,
                            98.99596252593278
                        ],
                        linf=[
                            0.0042295649514556555,
                            1.121978987805587,
                            0.6394799808438503,
                            1.1308638979672447,
                            1.0015178951935582,
                            333.1178018152714
                        ],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2),
                        volume_flux=(flux_ec, flux_nonconservative_waruszewski_etal))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_potential_temperature_baroclinic_instability Artiano" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "global_circulation",
                                 "elixir_potential_temperature_baroclinic_instability.jl"),
                        l2=[
                            0.0005237192887527135,
                            0.014655078153251792,
                            0.012853482849295788,
                            0.012929758648242088,
                            0.13994217088888805,
                            98.99596252593278
                        ],
                        linf=[
                            0.00424067775992687,
                            1.1141182921289354,
                            0.6517477285197878,
                            1.1292930037247177,
                            1.0016134087282467,
                            333.1178018152714
                        ],
                        rtol=2e-8,
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2),
                        volume_flux=(flux_etec, flux_nonconservative_artiano_etal))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_potential_temperature_held_suarez" begin
    import ..CI_ON_MACOS
    if CI_ON_MACOS
        global _rtol = 7e-8  # increased error tolerance
    else
        global _rtol = sqrt(eps(Float64))  # default
    end

    @test_trixi_include(joinpath(EXAMPLES_DIR, "global_circulation",
                                 "elixir_potential_temperature_held_suarez.jl"),
                        l2=[
                            0.0031433373482917877,
                            0.00013227403214817446,
                            0.0001322740320211912,
                            0.00014259768410822775,
                            0.7578318727895532,
                            569.8247547308886
                        ],
                        linf=[
                            0.023356419582470034,
                            0.001522627198332827,
                            0.0015226271932242787,
                            0.0005022230908559857,
                            4.88597072706591,
                            1703.946276059638
                        ],
                        rtol=_rtol,
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY),
                        lat_lon_trees_per_dim=2, layers=2)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

end
