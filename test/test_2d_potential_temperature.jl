module TestExamples2DEulerPotentialTemperature

include("test_trixiatmo.jl")

@trixi_testset "elixir_euler_potential_temperature_inertia_gravity_waves_2d" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_inertia_gravity_waves.jl"),
                        l2=[
                            5.360162314066379e-7,
                            0.00013097436752903342,
                            7.172813096420895e-5,
                            3.560199562542724e-5,
                            9.251160452856857e-12
                        ],
                        linf=[
                            3.7348732354614e-6,
                            0.0005402594042358544,
                            0.00023237719913361634,
                            0.00015331821185782246,
                            4.3655745685100555e-11
                        ], tspan=(0.0, 1800.0))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_euler_potential_temperature_linear_hydrostatic" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_linear_hydrostatic.jl"),
                        l2=[
                            2.669332229493225e-6,
                            0.0006518362209690377,
                            0.00011039657805452965,
                            0.00022353351014161769,
                            2.647509013539741e-11
                        ],
                        linf=[
                            2.7590006774769193e-5,
                            0.017497225180171938,
                            0.0018056132125011,
                            0.0010785330882754351,
                            1.7462298274040222e-10
                        ],
                        tspan=(0.0, 360.0), cells_per_dimension=(20, 12))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_euler_potential_temperature_linear_nonhydrostatic" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_linear_nonhydrostatic.jl"),
                        l2=[
                            4.345825921806206e-7,
                            0.00018887050422385925,
                            0.0001275556315244483,
                            8.576762629640317e-5,
                            2.7230966512302478e-11
                        ],
                        linf=[
                            1.0334905581110831e-5,
                            0.013797654432108786,
                            0.007887478371833376,
                            0.0003273023447150081,
                            1.4551915228366852e-10
                        ], tspan=(0.0, 360.0), cells_per_dimension=(20, 12))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_euler_potential_temperature_schaer_mountain" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_schaer_mountain.jl"),
                        l2=[
                            0.00010948264917453087,
                            0.09903481134935169,
                            0.09822804596320808,
                            0.004246742748260301,
                            1.8456576434476605e-11
                        ],
                        linf=[
                            0.0017570540164735249,
                            2.3793065145457373,
                            2.354152711825134,
                            0.07647977352115731,
                            8.731149137020111e-11
                        ],
                        tspan=(0.0, 360.0), cells_per_dimension=(20, 12))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

end
