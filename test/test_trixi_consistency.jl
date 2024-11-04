module TestTrixiConsistency

include("test_trixiatmo.jl")

EXAMPLES_DIR = TrixiAtmo.examples_dir()
TRIXI_EXAMPLES_DIR = Trixi.examples_dir()

# Start with a clean environment: remove Trixi.jl output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

@trixiatmo_testset "Dry air consistency check" begin
    # Dry air warm bubble test case in Trixi.jl
    maxiters = 100
    trixi_elixir = joinpath(TRIXI_EXAMPLES_DIR, "structured_2d_dgsem",
                            "elixir_euler_warm_bubble.jl")

    # Override maxiter and fluxes
    trixi_include(trixi_elixir;
                  volume_flux = Trixi.flux_chandrashekar,
                  surface_flux = Trixi.FluxLMARS(360.0),
                  maxiters = maxiters)

    # Save errors
    errors_trixi = Main.analysis_callback(Main.sol)

    # Create an instance of Trixi's equations, just used for dispatch below
    equations_trixi = Trixi.CompressibleEulerEquations2D(Main.warm_bubble_setup.gamma)
    add_zeros = SVector(zero(eltype(Main.sol)), zero(eltype(Main.sol)))

    # Now use moist equations instead
    equations_moist = CompressibleMoistEulerEquations2D()

    # Redefine source terms for CompressibleMoistEulerEquations2D
    @inline function (setup::Main.WarmBubbleSetup)(u, x, t,
                                                   equations_moist::CompressibleMoistEulerEquations2D)
        ret_trixi = setup(u, x, t, equations_trixi)
        return vcat(ret_trixi, add_zeros)
    end

    # Redefine initial condition for CompressibleMoistEulerEquations2D
    @inline function (setup::Main.WarmBubbleSetup)(x, t,
                                                   equations_moist::CompressibleMoistEulerEquations2D)
        ret_trixi = setup(x, t, equations_trixi)
        return vcat(ret_trixi, add_zeros)
    end

    # Run again with overrides
    trixi_include(trixi_elixir,
                  equations = equations_moist,
                  volume_flux = flux_chandrashekar,
                  surface_flux = Trixi.FluxLMARS(360.0),
                  maxiters = maxiters)

    errors_atmo = Main.analysis_callback(Main.sol)

    for (error_trixi, error_atmo) in zip(errors_trixi.l2, errors_atmo.l2)
        @test isapprox(error_trixi, error_atmo, rtol = 1e-12)
    end
    for (error_trixi, error_atmo) in zip(errors_trixi.linf, errors_atmo.linf)
        @test isapprox(error_trixi, error_atmo, rtol = 3e-12)
    end
end
end
