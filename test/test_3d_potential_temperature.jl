module TestExamples3DEulerPotentialTemperature

include("test_trixiatmo.jl")

@trixi_testset "elixir_euler_potential_temperature_taylor_green_vortex_3d" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_taylor_green_vortex.jl"),
                        l2=[
                            0.0071357545760608035,
                            0.05920029261132884,
                            0.05830264156640887,
                            0.09239774238361848,
                            0.002651207567127836
                        ],
                        linf=[
                            0.025211616325348185,
                            0.16263745267015323,
                            0.15672821384993219,
                            0.2914129381475647,
                            0.010284829976216603
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
                            0.007135730245594034,
                            0.05920029821961041,
                            0.05830265228349908,
                            0.09239771293546377,
                            0.0026511981540712527
                        ],
                        linf=[
                            0.02521107720818483,
                            0.1626373500658826,
                            0.15673005143018703,
                            0.29141239017096776,
                            0.010284751307188589
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
                            0.007135866414170671,
                            0.059200304019444996,
                            0.058302644968753516,
                            0.09239772516650108,
                            0.002651207110461276
                        ],
                        linf=[
                            0.025215717602879106,
                            0.16264021301187603,
                            0.15673141525377227,
                            0.291413219311949,
                            0.010284827215485737
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
                            0.047663959124221666,
                            0.050155417906555394,
                            0.030904595875059383,
                            0.012502671674568866
                        ],
                        tspan=(0.0, 0.2), surface_flux=FluxLMARS(340.0),
                        volume_flux=flux_tec)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_euler_potential_temperature_baroclinic_instability" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_baroclinic_instability.jl"),
                        l2=[
                            0.000540424089786219,
                            0.015509349290483546,
                            0.013816889458104184,
                            0.012939501902053414,
                            0.1400048323583557,
                            99.03486863938933
                        ],
                        linf=[
                            0.004654672800263215,
                            1.1039457129713743,
                            0.6663687877166378,
                            1.1276246706650594,
                            1.0016884572416416,
                            333.2487183585763
                        ],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_euler_potential_temperature_baroclinic_instability" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_baroclinic_instability.jl"),
                        l2=[
                            0.0005197412073779303,
                            0.014184223075386094,
                            0.01231640074505267,
                            0.012923945625897474,
                            0.1399501917281505,
                            99.03486863938933
                        ],
                        linf=[
                            0.004231445262675715,
                            1.12197009874999,
                            0.639464719150829,
                            1.1308556842535742,
                            1.0018956805357107,
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
                            0.0005238483791325817,
                            0.014658767813775264,
                            0.012858058446996147,
                            0.012928211199993971,
                            0.13996695188603392,
                            99.03486863938933
                        ],
                        linf=[
                            0.004242592531761069,
                            1.1140991059698901,
                            0.6517475169161369,
                            1.129285594724997,
                            1.0019910798637852,
                            333.2487183585763
                        ],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2),
                        volume_flux=(flux_etec, flux_nonconservative_artiano_etal))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

end
