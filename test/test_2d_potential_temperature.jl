module TestExamples2DEulerPotentialTemperature

include("test_trixiatmo.jl")

@trixi_testset "elixir_euler_potential_temperature_inertia_gravity_waves_2d" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_inertia_gravity_waves.jl"),
                        l2=[
                            5.360162314518395e-7,
                            0.0001309743674925418,
                            7.172813093057917e-5,
                            3.560199562285229e-5,
                            9.251160452856857e-12
                        ],
                        linf=[
                            3.734873233907088e-6,
                            0.0005402594041576947,
                            0.0002323771992720964,
                            0.00015331821202835272,
                            4.3655745685100555e-11
                        ],
                        tspan=(0.0, 1800.0))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_euler_potential_temperature_linear_hydrostatic" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_linear_hydrostatic.jl"),
                        l2=[
                            2.6693322294197e-6,
                            0.0006518362209620777,
                            0.00011039657710646231,
                            0.00022353351014774928,
                            2.647509013539741e-11
                        ],
                        linf=[
                            2.7590006778765996e-5,
                            0.01749722523658548,
                            0.0018056125252837996,
                            0.0010785335742866664,
                            1.7462298274040222e-10
                        ],
                        T=0.1, cells_per_dimension=(20, 12))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_euler_potential_temperature_linear_nonhydrostatic" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_linear_nonhydrostatic.jl"),
                        l2=[
                            4.3458258873570067e-7,
                            0.00018887050881119427,
                            0.00012755554165085224,
                            8.57676250982478e-5,
                            2.7230966512302478e-11
                        ],
                        linf=[
                            1.0335035563580064e-5,
                            0.013797661538720973,
                            0.007887441700822817,
                            0.0003273012893600935,
                            1.4551915228366852e-10
                        ],
                        tspan=(0.0, 360.0), cells_per_dimension=(20, 12))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_euler_potential_temperature_schaer_mountain" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_schaer_mountain.jl"),
                        l2=[
                            0.00010948264917559538,
                            0.0990348113475616,
                            0.09822804595816467,
                            0.004246742749254792,
                            1.8456576434476605e-11
                        ],
                        linf=[
                            0.0017570541001279416,
                            2.379306511006467,
                            2.3541527047673907,
                            0.07647979859615361,
                            8.731149137020111e-11
                        ],
                        tspan=(0.0, 360.0), cells_per_dimension=(20, 12))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

end
