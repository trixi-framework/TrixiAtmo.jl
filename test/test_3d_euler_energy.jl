module TestExamples3DEulerEnergy
using Trixi: flux_kennedy_gruber
include("test_trixiatmo.jl")

@trixi_testset "elixir_euler_energy_baroclinic_instability" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_euler_energy_baroclinic_instability.jl"),
                        l2=[
                            6.552421954615372e-7,
                            0.00019731220715342904,
                            0.00042455589284480093,
                            0.0001871831573805277,
                            0.06704750535405761,
                            0.0409097168117604
                        ],
                        linf=[
                            4.428967270930784e-6,
                            0.024276204894145392,
                            0.045702372863804186,
                            0.022711918877967086,
                            0.9494724034739193,
                            0.3119894564151764
                        ], tspan=(0.0, 0.01), tol=1e-15, atol=1e-8)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "elixir_euler_energy_baroclinic_instability" begin
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR,
                           "elixir_euler_energy_baroclinic_instability.jl"),
                  volume_flux = (flux_kennedy_gruber, flux_nonconservative_souza_etal),
                  tspan = (0.0, 0.01))
    u_ode = copy(sol.u[end])
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR,
                           "elixir_euler_energy_baroclinic_instability_turbo.jl"),
                  tspan = (0.0, 0.01))

    u_ode_specialized = copy(sol.u[end])
    @test u_ode_specialized ≈ u_ode
end

end
