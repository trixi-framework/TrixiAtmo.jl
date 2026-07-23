module TestExamples2DEulerEnergy

include("test_trixiatmo.jl")

EXAMPLES_DIR = joinpath(EXAMPLES_DIR, "euler/dry_air")

@trixi_testset "elixir_energy_inertia_gravity_waves" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "buoyancy",
                                 "elixir_energy_inertia_gravity_waves.jl"),
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

@trixi_testset "elixir_covariant_energy_inertia_gravity_waves" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "buoyancy",
                                 "elixir_covariant_energy_inertia_gravity_waves.jl"),
                        l2=[
                            5.988113879663543e-8,
                            1.652979410541983e-9,
                            5.413741859342455e-8,
                            0.016802245219125052
                        ],
                        linf=[
                            7.418985343843332e-7,
                            2.078836479520174e-8,
                            5.333954417018675e-7,
                            0.18759511329699308
                        ], tspan=(0.0, 10.0))
end
end
