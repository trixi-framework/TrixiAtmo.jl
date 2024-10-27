module TestSphericalAdvection

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

EXAMPLES_DIR = pkgdir(TrixiAtmo, "examples")

@trixiatmo_testset "Spherical advection, Cartesian weak form, LLF surface flux" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_spherical_advection_cartesian.jl"),
                        l2=[
                            0.008445050734443814,
                            0.008234141170177999,
                            0.0018421064842753803,
                            0.0,
                            0.06443024298066,
                        ],
                        linf=[
                            0.09489504882003308,
                            0.09648119517743559,
                            0.013745339970863746,
                            0.0,
                            0.40932299941471895,
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

@trixiatmo_testset "Spherical advection, covariant weak form, LLF surface flux" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_spherical_advection_covariant.jl"),
                        l2=[1.0001621291410059, 0.0, 0.0],
                        linf=[14.245184617343739, 0.0, 0.0])
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

# The covariant flux-differencing form should be equivalent to the weak form when the 
# arithmetic mean is used as the two-point flux
@trixiatmo_testset "Spherical advection, covariant flux-differencing, central/LLF" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_spherical_advection_covariant.jl"),
                        l2=[1.000162129141023, 0.0, 0.0],
                        linf=[14.245184617342034, 0.0, 0.0],
                        volume_integral=VolumeIntegralFluxDifferencing(flux_central))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

# Version with arithmetic mean used for both the volume and surface fluxes
@trixiatmo_testset "Spherical advection, covariant flux-differencing, central/central" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_spherical_advection_covariant.jl"),
                        l2=[2.4999822899157116, 0.0, 0.0],
                        linf=[38.09056110399384, 0.0, 0.0],
                        volume_integral=VolumeIntegralFluxDifferencing(flux_central),
                        surface_flux=flux_central)
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
