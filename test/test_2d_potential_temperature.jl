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

end
