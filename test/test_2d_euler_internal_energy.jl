module TestExamples2DEulerInternalEnergy

include("test_trixiatmo.jl")

EXAMPLES_DIR = joinpath(EXAMPLES_DIR, "euler/dry_air")

@trixi_testset "elixir_internal_energy_inertia_gravity_waves" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "buoyancy",
                                 "elixir_energy_inertia_gravity_waves.jl"),
                        l2=[
                            2.3801001547126282e-7,
                            6.190703738154071e-6,
                            3.382129115733546e-5,
                            0.04382810096988029,
                            9.251160452856857e-12
                        ],
                        linf=[
                            1.2362858290426715e-6,
                            6.0616987330064376e-5,
                            0.00033470830317765867,
                            0.27272386007825844,
                            4.3655745685100555e-11
                        ], tspan=(0.0, 10.0), atol=5e-11)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

end
