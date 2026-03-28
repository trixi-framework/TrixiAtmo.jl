module TestExamples2DEulerInternalEnergy

include("test_trixiatmo.jl")

EXAMPLES_DIR = joinpath(EXAMPLES_DIR, "euler/dry_air")

@trixi_testset "elixir_internal_energy_inertia_gravity_waves" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "buoyancy",
                                 "elixir_energy_inertia_gravity_waves.jl"),
                        l2=[
                            2.380100116660882e-7,
                            6.190675997897145e-6,
                            3.382129329051498e-5,
                            0.04378077657190983,
                            9.251160452856857e-12
                        ],
                        linf=[
                            1.2361079473333092e-6,
                            6.0614141304426994e-5,
                            0.0003347340017250148,
                            0.2723903215082828,
                            4.3655745685100555e-11
                        ], tspan=(0.0, 10.0), atol=5e-11)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

end
