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
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@trixi_testset "elixir_covariant_energy_inertia_gravity_waves" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "buoyancy",
                                 "elixir_covariant_energy_inertia_gravity_waves.jl"),
                        l2=[
                            8.03531577207334e-8,
                            1.799854544366393e-9,
                            5.41198198917674e-8,
                            0.01698674060190595
                        ],
                        linf=[
                            1.9031079829812114e-6,
                            3.44587946336794e-8,
                            5.33957472269021e-7,
                            0.18672661460004747
                        ], tspan=(0.0, 10.0))
end
end
