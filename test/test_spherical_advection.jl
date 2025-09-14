module TestSphericalAdvection

include("test_trixiatmo.jl")

@trixi_testset "Spherical advection (cubed sphere), Cartesian weak form, LLF surface flux" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cartesian_advection_cubed_sphere.jl"),
                        l2=[
                            0.796321633847963,
                            20.317829852214242,
                            8.810001095522356,
                            20.317829852220424,
                            0.0
                        ],
                        linf=[
                            10.872101732112924,
                            289.6515963627462,
                            95.1288711914458,
                            289.65159636274984,
                            0.0
                        ])
    # and small reference values
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "Spherical advection (quad icosahedron), Cartesian weak form, LLF surface flux" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cartesian_advection_quad_icosahedron.jl"),
                        l2=[
                            0.45702277148770143,
                            11.807355540181147,
                            4.311881740807178,
                            11.807355540181314,
                            0.0
                        ],
                        linf=[
                            13.591965583195247,
                            364.76418895378083,
                            93.69731833987953,
                            364.7641889537881,
                            0.0
                        ])
    # and small reference values
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "Spherical advection, Cartesian weak form, element-local mapping" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_shallowwater_cartesian_advection_cubed_sphere.jl"),
                        l2=[
                            0.893342967293854,
                            22.848879918993113,
                            9.758850586740538,
                            22.848879918993124,
                            0.0
                        ],
                        linf=[
                            14.289456304585201,
                            380.6958334056544,
                            120.59259301568181,
                            380.69583340557074,
                            0.0
                        ], element_local_mapping=true)
    # and small reference values
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "Spherical advection, covariant weak form, LLF surface flux" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_spherical_advection_covariant_cubed_sphere.jl"),
                        l2=[1.0007043506351705, 0.0, 0.0],
                        linf=[14.235905681508598, 0.0, 0.0])
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "Spherical advection, covariant weak form, LLF surface flux, global spherical coords" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_spherical_advection_covariant_cubed_sphere.jl"),
                        l2=[1.0007043506351705, 0.0, 0.0],
                        linf=[14.235905681508598, 0.0, 0.0],
                        global_coordinate_system=GlobalSphericalCoordinates())
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

# The covariant flux-differencing form should be equivalent to the weak form when the 
# arithmetic mean is used as the two-point flux
@trixi_testset "Spherical advection, covariant flux-differencing, central/LLF" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_spherical_advection_covariant_cubed_sphere.jl"),
                        l2=[1.0007043506351412, 0.0, 0.0],
                        linf=[14.23590568150928, 0.0, 0.0],
                        volume_integral=VolumeIntegralFluxDifferencing(flux_central))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

# Version with arithmetic mean used for both the volume and surface fluxes
@trixi_testset "Spherical advection, covariant flux-differencing, central/central" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_spherical_advection_covariant_cubed_sphere.jl"),
                        l2=[2.499889861385917, 0.0, 0.0],
                        linf=[38.085244441156085, 0.0, 0.0],
                        volume_integral=VolumeIntegralFluxDifferencing(flux_central),
                        surface_flux=flux_central)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

@trixi_testset "Spherical advection on icosahedral grid, covariant weak form, LLF surface flux" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_spherical_advection_covariant_quad_icosahedron.jl"),
                        l2=[0.5183886767005157, 0.0, 0.0],
                        linf=[13.54834739856517, 0.0, 0.0])
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs!, semi, sol, 100)
end

end # module
