@testsnippet EulerEnergy2D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "euler", "dry_air")
end

@testitem "Euler energy 2D: elixir_energy_inertia_gravity_waves" setup=[
    Setup,
    EulerEnergy2D
] tags=[:euler_energy_2d] begin
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

@testitem "Euler energy 2D: elixir_covariant_energy_inertia_gravity_waves" setup=[
    Setup,
    EulerEnergy2D
] tags=[:euler_energy_2d] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "buoyancy",
                                 "elixir_covariant_energy_inertia_gravity_waves.jl"),
                        l2=[
                            5.986441694838273e-8,
                            1.6527854310546966e-9,
                            5.412117115516673e-8,
                            0.01678917928317668
                        ],
                        linf=[
                            7.415741256622255e-7,
                            2.0783123345566312e-8,
                            5.332369591815093e-7,
                            0.18741797714028507
                        ], tspan=(0.0, 10.0))
end
