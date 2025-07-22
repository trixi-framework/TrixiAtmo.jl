module TestSphericalAdvection

using Test
using TrixiAtmo

include("test_trixiatmo.jl")

EXAMPLES_DIR = TrixiAtmo.examples_dir()

@trixiatmo_testset "Spherical advection (cubed sphere), Cartesian weak form, LLF surface flux" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cartesian_advection_cubed_sphere.jl"),
                        l2=[
                            0.183524336,
                            5.64470440,
                            2.80467667,
                            5.64470440,
                            0.0
                        ],
                        linf=[
                            2.79650356,
                            91.6208734,
                            24.9330974,
                            91.6208734,
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
                                 "elixir_shallowwater_cartesian_advection_quad_icosahedron.jl"),
                        l2=[
                            0.0771325788,
                            2.24899209,
                            1.45355192,
                            2.24899209,
                            0.0
                        ],
                        linf=[
                            1.18929523,
                            34.1427583,
                            17.9760460,
                            34.1427583,
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
                                 "elixir_shallowwater_cartesian_advection_cubed_sphere.jl"),
                        l2=[
                            0.222179190,
                            6.74419737,
                            3.36551000,
                            6.74419737,
                            0.0
                        ],
                        linf=[
                            3.60436608,
                            113.572183,
                            30.1250336,
                            113.572183,
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
                                 "elixir_spherical_advection_covariant_cubed_sphere.jl"),
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

@trixiatmo_testset "Spherical advection, covariant weak form, LLF surface flux, global spherical coords" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_spherical_advection_covariant_cubed_sphere.jl"),
                        l2=[1.0007043506351705, 0.0, 0.0],
                        linf=[14.235905681508598, 0.0, 0.0],
                        global_coordinate_system=GlobalSphericalCoordinates())
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
                                 "elixir_spherical_advection_covariant_cubed_sphere.jl"),
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
                                 "elixir_spherical_advection_covariant_cubed_sphere.jl"),
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

@trixiatmo_testset "Spherical advection on icosahedral grid, covariant weak form, LLF surface flux" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_spherical_advection_covariant_quad_icosahedron.jl"),
                        l2=[0.5183886767005157, 0.0, 0.0],
                        linf=[13.54834739856517, 0.0, 0.0])
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
