module TestExamples3DEulerEnergy

include("test_trixiatmo.jl")

@trixi_testset "elixir_euler_energy_baroclinic_instability" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_energy_baroclinic_instability.jl"),
                        l2=[
                            6.554982393176306e-7,
                            0.00019734512089352276,
                            0.000424563812892216,
                            0.00018719200535314822,
                            0.06706018255246633,
                            0.0409097168117604
                        ],
                        linf=[
                            4.429177855591604e-6,
                            0.024280872072865722,
                            0.04570529810659252,
                            0.022704417775103652,
                            0.9493969184113666,
                            0.3119894564151764
                        ], tspan=(0.0, 0.1), tol=1e-13)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

end
