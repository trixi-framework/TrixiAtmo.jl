module TestExamples2DRainyEuler

include("test_trixiatmo.jl")

EXAMPLES_DIR = joinpath(EXAMPLES_DIR, "euler/cartesian")

@trixi_testset "convergence_test" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "precipitation/tests",
                                 "convergence_test.jl"),
                        l2=[
                            2.39895785368954e-6,
                            3.892427386570826e-5,
                            1.5186141589998252e-5,
                            0.007532307120716242,
                            0.02085654062623111,
                            51.988569260688045,
                            1.7470730680079724e-8,
                            3.892376665583286e-5,
                            0.0002156633044466944
                        ],
                        linf=[
                            1.675032502916618e-5,
                            0.0002740479581908595,
                            8.691475272920579e-5,
                            0.023414703432404593,
                            0.06248365429859781,
                            385.10196906188503,
                            6.003313225070965e-8,
                            0.000274026861750043,
                            0.0007433951467703537
                        ])
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixi_testset "elixir_hoeck_bubble moist" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "moist_air/buoyancy",
                                 "elixir_hoeck_bubble.jl"),
                        l2=[
                            0.0031469268543095233,
                            6.293853704548511e-5,
                            0.0,
                            0.030085941482251792,
                            1.2721823851437415,
                            1040.2702923456666,
                            0.0011556742959893704,
                            0.0011991088277167412,
                            2.080063951506559
                        ],
                        linf=[
                            0.010772626708125621,
                            0.00021545552388836306,
                            0.0,
                            0.17969745708245988,
                            2.469797669783883,
                            3725.8780578381848,
                            0.0030306571979925243,
                            0.0032258070568836963,
                            3.242378389076862
                        ],
                        cells_per_dimension=(64, 32),
                        tspan=(0.0, 10.0))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixi_testset "elixir_hoeck_bubble rainy" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "precipitation/buoyancy",
                                 "elixir_hoeck_bubble.jl"),
                        l2=[
                            7.959523735914366e-5,
                            4.8983178144498016e-5,
                            1.2695889115092138e-9,
                            0.019951374355214414,
                            0.033830586959875855,
                            126.25361431033272,
                            4.899300381610471e-5,
                            3.4512985541873764e-7,
                            0.0016455875546178209
                        ],
                        linf=[
                            0.0011817049726382534,
                            0.0007153644712404369,
                            2.0653953346431095e-8,
                            0.12229462313083396,
                            0.1760334548347,
                            1831.4029281113471,
                            0.0007153644712404369,
                            6.113680304613196e-6,
                            0.008505948056836132
                        ],
                        cells_per_dimension=(16, 16),
                        tspan=(0.0, 10.0))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

# For unknown reasons, github's macos runners produce results exceeding the default
# tolerance
@trixi_testset "elixir_hoeck_bubble_diffusion rainy" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "precipitation/buoyancy",
                                 "elixir_hoeck_bubble_diffusion.jl"),
                        l2=[
                            8.025606283886885e-5,
                            4.9286718733941776e-5,
                            1.1772196106012233e-9,
                            0.019937188110969426,
                            0.03381156061689547,
                            126.9903604077802,
                            4.9265000572589546e-5,
                            3.0236037496562846e-7,
                            0.0016285128672919942
                        ],
                        linf=[
                            0.0011872217179234035,
                            0.0007188675819966152,
                            2.6738841642904725e-8,
                            0.12197434250118369,
                            0.17602879576952665,
                            1840.4575813068077,
                            0.0007188675819966152,
                            7.048542362862916e-6,
                            0.00874081765715573
                        ],
                        rtol=6e-6,
                        initial_refinement_level=4,
                        tspan=(0.0, 10.0))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end
end # module
