module TestExamples1DEulerPotentialTemperature

include("test_trixiatmo.jl")

@trixi_testset "elixir_euler_potential_temperature_ec_1d" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_ec.jl"),
                        l2 = [1.536920933776253, 1.536920944343481, 8.204456414321612e-7], 
                        linf = [2.283573460625467, 2.2835683135782463, 4.835357438226495e-6],
                        tspan=(0.0, 0.4))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

end