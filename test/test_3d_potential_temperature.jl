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

@trixi_testset "elixir_euler_potential_temperature_baroclinic_instability" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_baroclinic_instability.jl"),
                        l2=[
                            3.6191804848984097e-6,
                            0.0007022602048939841,
                            0.0009478663759665131,
                            0.000772072851465165,
                            0.00045473753176575744,
                            0.0409097168117604
                        ],
                        linf=[
                            0.00018385259000197607,
                            0.08088210121911142,
                            0.11309008780578583,
                            0.07195098637367636,
                            0.04285569464735772,
                            0.3119894564151764
                        ],
                        tspan=(0.0, 0.01*SECONDS_PER_DAY))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

end
