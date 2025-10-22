module TestExamples2DEulerEnergy

include("test_trixiatmo.jl")

@trixi_testset "elixir_euler_energy_inertia_gravity_waves_2d" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_energy_inertia_gravity_waves.jl"),
                        l2=[
                            2.3800999105272615e-7,
                            6.190703408721927e-6,
                            3.3821288686132935e-5,
                            0.04382810097184301,
                            9.251160452856857e-12
                        ],
                        linf=[
                            1.2362858294867607e-6,
                            6.0616987290984525e-5,
                            0.00033470830315131473,
                            0.27272386016556993,
                            4.3655745685100555e-11
                        ], tspan=(0.0, 10.0), atol=5e-11)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

end
