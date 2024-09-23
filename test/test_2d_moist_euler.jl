module TestExamples2DMoistEuler

using Test
using TrixiAtmo

include("test_trixiatmo.jl") # TODO - This is a repetition from Trixi.jl

EXAMPLES_DIR = pkgdir(TrixiAtmo, "examples") # TODO - Do we need a subdirectory for examples?

@trixiatmo_testset "elixir_moist_euler_dry_bubble" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_moist_euler_dry_bubble.jl"),
                        l2=[
                            1.300428671901329e-6,
                            2.601090012108739e-5,
                            0.0006660124630171347,
                            0.008969786054960861,
                            0.0,
                            0.0,
                        ],
                        linf=[
                            1.0312042909910168e-5,
                            0.00020625488871672815,
                            0.006392107590872236,
                            0.07612038028310053,
                            0.0,
                            0.0,
                        ],
                        polydeg=3,
                        tspan=(0.0, 0.1))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "elixir_moist_euler_EC_bubble" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_moist_euler_EC_bubble.jl"),
                        l2=[
                            0.01345154393018332,
                            0.8070926361417218,
                            7.938812668709457,
                            4500.747616411578,
                            0.00015592413050814787,
                            0.00014163475049532796,
                        ],
                        linf=[
                            0.1427479052298466,
                            8.564879578662357,
                            91.56822550162855,
                            49528.383866247605,
                            0.0019364397602254623,
                            0.0013259689889851285,
                        ],
                        polydeg=3,
                        cells_per_dimension=(16, 16),
                        tspan=(0.0, 0.1))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "elixir_moist_euler_moist_bubble" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_moist_euler_moist_bubble.jl"),
                        l2=[
                            7.351043427240923e-6,
                            1.1070342432069074e-7,
                            0.0006974588377288118,
                            1.715668353329522,
                            8.831269143134121e-7,
                            1.025579538944668e-6,
                        ],
                        linf=[
                            8.055695643149896e-5,
                            1.1985203677080201e-6,
                            0.005897639251024437,
                            19.24776030163048,
                            1.0043133039065386e-5,
                            1.1439046776775402e-5,
                        ],
                        polydeg=3,
                        cells_per_dimension=(16, 8),
                        tspan=(0.0, 0.1))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "elixir_moist_euler_nonhydrostatic_gravity_waves" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_moist_euler_nonhydrostatic_gravity_waves.jl"),
                        l2=[
                            3.5420405147937216e-5,
                            0.0021265774152361538,
                            0.01172830034500581,
                            9.898560584459009,
                            0.0,
                            0.0,
                        ],
                        linf=[
                            0.0017602202683439927,
                            0.14941973735192882,
                            0.5945327351674782,
                            489.89171406268724,
                            0.0,
                            0.0,
                        ],
                        polydeg=3,
                        cells_per_dimension=(10, 8),
                        tspan=(0.0, 0.1))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "elixir_moist_euler_source_terms_dry" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_moist_euler_source_terms_dry.jl"),
                        l2=[
                            1.3992076791281227e-5,
                            1.4486417486907815e-5,
                            2.609465609739115e-5,
                            6.323484379066432e-5,
                            0.0,
                            0.0,
                        ],
                        linf=[
                            7.549984224430872e-5,
                            0.00010065352517929504,
                            0.00015964938767742964,
                            0.0005425860570698049,
                            0.0,
                            0.0,
                        ],
                        polydeg=3,
                        cells_per_dimension=(10, 8),
                        tspan=(0.0, 0.1))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "elixir_moist_euler_source_terms_moist" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_moist_euler_source_terms_moist.jl"),
                        l2=[
                            1.8307663991129928e-5,
                            0.04008077097727512,
                            0.015104690877128945,
                            0.5242098451814421,
                            5.474006801215573e-10,
                            1.1103883907226752e-10,
                        ],
                        linf=[
                            0.00013219484616722177,
                            0.10771224937484702,
                            0.03789645369775574,
                            3.90888311638264,
                            3.938382289041286e-9,
                            6.892033377287209e-10,
                        ],
                        polydeg=3,
                        cells_per_dimension=(10, 8),
                        tspan=(0.0, 0.1))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "elixir_moist_euler_source_terms_split_moist" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_moist_euler_source_terms_split_moist.jl"),
                        l2=[
                            0.0001480393848825987,
                            0.11945481031503036,
                            0.07460345535073129,
                            5.943049264963717,
                            4.471792794168725e-9,
                            7.10320253652373e-10,
                        ],
                        linf=[
                            0.0007084183215528839,
                            0.5047962996690205,
                            0.3697160082709827,
                            27.843155286573165,
                            2.1168438904322837e-8,
                            3.691699932047233e-9,
                        ],
                        polydeg=3,
                        cells_per_dimension=(10, 8),
                        tspan=(0.0, 0.1))
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
