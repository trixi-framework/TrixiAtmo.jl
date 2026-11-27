@testsnippet EulerInternalEnergy2D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "euler", "dry_air")
end

@testitem "Euler internal energy 2D: elixir_internal_energy_inertia_gravity_waves" setup=[
    Setup,
    EulerInternalEnergy2D
] tags=[:euler_internal_energy_2d] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "buoyancy",
                                 "elixir_internal_energy_inertia_gravity_waves.jl"),
                        l2=[
                            2.3801001547126282e-7,
                            6.190703738154071e-6,
                            3.382129115733546e-5,
                            0.04378077771265
                        ],
                        linf=[
                            1.2361505137281625e-6,
                            6.0611579968394835e-5,
                            0.0003347311706169429,
                            0.27240109792910516
                        ], tspan=(0.0, 10.0), atol=5e-11)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end
