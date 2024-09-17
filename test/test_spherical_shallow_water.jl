module TestSphericalShallowWater

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

EXAMPLES_DIR = pkgdir(TrixiAtmo, "examples")

@trixiatmo_testset "Spherical shallow water equations, covariant weak form" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_spherical_shallow_water_covariant.jl"),
                        l2=[
                            3.06531384e-01,
                            1.99844676e-04,
                            8.76782001e-05,
                        ],
                        linf=[
                            1.47864310e+00,
                            1.37545965e-03,
                            7.56402953e-04,
                        ], tspan=(0.0, 1.0 * SECONDS_PER_DAY))
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
