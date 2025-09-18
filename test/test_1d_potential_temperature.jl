module TestExamples1DEulerPotentialTemperature

include("test_trixiatmo.jl")

@trixi_testset "elixir_euler_potential_temperature_ec" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_ec.jl"),
                        l2=[1.536920933776253, 1.536920944343481, 8.204456414321612e-7],
                        linf=[2.283573460625467, 2.2835683135782463, 4.835357438226495e-6],
                        tspan=(0.0, 0.4), rtol=10 * sqrt(eps(Float64)))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_euler_potential_temperature_tec" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_ec.jl"),
                        l2=[1.5369209174209315, 1.5369209181196604, 1.9668131414528186e-7],
                        linf=[2.2835833418413567, 2.283588054059723, 1.1052571658176635e-6],
                        tspan=(0.0, 0.4), surface_flux=flux_tec, volume_flux=flux_tec,
                        rtol=10 * sqrt(eps(Float64)))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_euler_potential_temperature_etec" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_potential_temperature_ec.jl"),
                        l2=[1.5369209174090952, 1.5369209159246773, 4.6312969205197763e-7],
                        linf=[2.2835901549809057, 2.2836015746559717, 2.547940256003578e-6],
                        tspan=(0.0, 0.4), surface_flux=flux_etec, volume_flux=flux_etec,
                        rtol=10 * sqrt(eps(Float64)))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

end
