module TestSphericalAdvection

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

EXAMPLES_DIR = pkgdir(TrixiAtmo, "examples")

@trixiatmo_testset "Spherical advection (cubed sphere), Cartesian weak form, LLF surface flux" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cubed_sphere_shell_advection.jl"),
                        l2=[
                            0.796321633853675,
                            20.317829852384286,
                            8.810001095524816,
                            20.317829852393054,
                            0.0
                        ],
                        linf=[
                            10.872101731709677,
                            289.6515963524798,
                            95.1288712006542,
                            289.65159635247255,
                            0.0
                        ])
    # and small reference values
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "Spherical advection (quad icosahedron), Cartesian weak form, LLF surface flux" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_quad_icosahedron_shell_advection.jl"),
                        l2=[
                            0.4570227714893776,
                            11.807355540221533,
                            4.311881740776439,
                            11.807355540221321,
                            0.0
                        ],
                        linf=[
                            13.591965583197862,
                            364.7641889538827,
                            93.69731833991227,
                            364.76418895387906,
                            0.0
                        ])
    # and small reference values
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated TrixiAtmo.Trixi.rhs!(du_ode, u_ode, semi, t)) < 100
    end
end

@trixiatmo_testset "Spherical advection, Cartesian weak form, element-local mapping" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cubed_sphere_shell_advection.jl"),
                        l2=[
                            0.893342967294085,
                            22.84887991899238,
                            9.75885058673732,
                            22.84887991899168,
                            0.0
                        ],
                        linf=[
                            14.289456304679561,
                            380.6958334082083,
                            120.59259301648672,
                            380.695833408201,
                            0.0
                        ], element_local_mapping=true)
    # and small reference values
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
                        l2=[1.0007043506351705, 0.0, 0.0],
                        linf=[14.235905681508598, 0.0, 0.0])
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
                        l2=[1.0007043506351412, 0.0, 0.0],
                        linf=[14.23590568150928, 0.0, 0.0],
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
                        l2=[2.499889861385917, 0.0, 0.0],
                        linf=[38.085244441156085, 0.0, 0.0],
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
