module TestTrixiConsistency

include("test_trixiatmo.jl")

# Start with a clean environment: remove Trixi.jl output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

@trixi_testset "Euler bubble" begin
    using Trixi
    using TrixiAtmo: CompressibleMoistEulerEquations2D

    # Dry air warm bubble test case in Trixi.jl
    maxiters = 100
    trixi_elixir = joinpath(Trixi.examples_dir(), "tree_2d_dgsem",
                            "elixir_euler_warm_bubble.jl")

    # Override fluxes, polydeg, cfl, maxiters
    trixi_include(trixi_elixir;
                  volume_flux = Trixi.flux_chandrashekar,
                  surface_flux = Trixi.FluxLMARS(360.0),
                  polydeg = 4,
                  stepsize_callback = Trixi.StepsizeCallback(cfl = 0.2),
                  maxiters = maxiters)

    # Save errors
    errors_trixi = Main.analysis_callback(Main.sol)

    # Create an instance of Trixi's equations, just used for dispatch below
    equations_trixi = Trixi.CompressibleEulerEquations2D(Main.warm_bubble_setup.gamma)
    add_zeros = SVector(zero(eltype(Main.sol)), zero(eltype(Main.sol)))

    # Redefine initial condition in Trixi.jl for CompressibleMoistEulerEquations2D
    # Different formulae were used!
    @inline function (setup::Main.WarmBubbleSetup)(x, t,
                                                   ::CompressibleMoistEulerEquations2D)
        ret_trixi = setup(x, t, equations_trixi)
        return vcat(ret_trixi, add_zeros)
    end

    # Now use the elixir in TrixiAtmo
    elixir_atmo = joinpath(EXAMPLES_DIR, "euler/cartesian/dry_air/buoyancy",
                           "elixir_gemein_bubble.jl")

    # Override initial condition, maxiters,
    # gravitational acceleration constant to match Trixi's equations
    trixi_include(elixir_atmo;
                  initial_condition = Main.warm_bubble_setup,
                  gravity = 9.81,
                  maxiters = maxiters)

    # Save errors
    errors_atmo = Main.analysis_callback(Main.sol)

    for (error_trixi, error_atmo) in zip(errors_trixi.l2, errors_atmo.l2)
        @test isapprox(error_trixi, error_atmo, rtol = 1e-12)
    end
    for (error_trixi, error_atmo) in zip(errors_trixi.linf, errors_atmo.linf)
        @test isapprox(error_trixi, error_atmo, rtol = 1.1e-10)
    end
end

end
