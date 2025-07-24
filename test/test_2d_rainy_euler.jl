module TestExamples2DRainyEuler

include("test_trixiatmo.jl")

EXAMPLES_DIR = joinpath(EXAMPLES_DIR, "moist_euler")

@trixi_testset "convergence_test_rainy" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "convergence_test",
                                 "convergence_test_rainy.jl"),
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

@trixi_testset "elixir_rainy_euler_moist_bubble" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "moist_bubble",
                                 "elixir_rainy_euler_moist_bubble.jl"),
                        l2=[
                            0.003147716446423741,
                            6.295432888771655e-5,
                            0.0,
                            0.030093419727259596,
                            1.2721952635118658,
                            1040.525033262258,
                            0.0011553981102543557,
                            0.0011988489480496496,
                            2.079644612712912
                        ],
                        linf=[
                            0.01077582487256823,
                            0.00021551948971416587,
                            0.0,
                            0.17975465277331062,
                            2.4701908567341215,
                            3726.9826530935243,
                            0.003030525617948497,
                            0.003225736615268841,
                            3.2422474910159167
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

@trixi_testset "elixir_rainy_euler_rainy_bubble" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "rainy_bubble",
                                 "elixir_rainy_euler_rainy_bubble.jl"),
                        l2=[
                            7.958084854135629e-5,
                            4.8974003075933214e-5,
                            1.2699485312744632e-9,
                            0.01995521719771058,
                            0.03383695227022311,
                            126.22998495250405,
                            4.8983831176389294e-5,
                            3.4524261261367737e-7,
                            0.0016458343657569546
                        ],
                        linf=[
                            0.001181534085815561,
                            0.0007152546342687649,
                            2.065750070675263e-8,
                            0.12231969676871683,
                            0.176070286893087,
                            1831.1175623952004,
                            0.0007152546342687649,
                            6.114801650025656e-6,
                            0.008508910564898997
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
@trixi_testset "elixir_rainy_euler_rainy_bubble_diffusion" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "rainy_bubble",
                                 "elixir_rainy_euler_rainy_bubble_diffusion.jl"),
                        l2=[
                            8.024161547037353e-5,
                            4.927749867325615e-5,
                            1.1775025026817308e-9,
                            0.019941027438108397,
                            0.03381791952397996,
                            126.96661069448412,
                            4.925578049099529e-5,
                            3.0245852299389407e-7,
                            0.0016287477694402524
                        ],
                        linf=[
                            0.001187049337148749,
                            0.0007187571931339255,
                            2.67415192177334e-8,
                            0.12199924771843951,
                            0.1760655719865429,
                            1840.1706161463226,
                            0.0007187571931339255,
                            7.049541733574276e-6,
                            0.008742438145986853
                        ],
                        rtol=9e-8,
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
