module TestExamples3DEulerEnergy

include("test_trixiatmo.jl")

EXAMPLES_DIR = joinpath(EXAMPLES_DIR, "euler/dry_air")

@trixi_testset "elixir_energy_baroclinic_instability" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "global_circulation",
                                 "elixir_energy_baroclinic_instability.jl"),
                        l2=[
                            6.540827547514395e-7,
                            0.0001972937313811121,
                            0.00042458497371703434,
                            0.00018720446122477726,
                            0.06698234918918748,
                            0.0408936452021324
                        ],
                        linf=[
                            4.42523425103758e-6,
                            0.024276205576123866,
                            0.04570236938186356,
                            0.022711920085996914,
                            0.9490993925137445,
                            0.3118668645620346
                        ], tspan=(0.0, 0.01), tol=1e-15, atol=1e-8)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

end
