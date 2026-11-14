@testsnippet Visualization begin
    EXAMPLES_DIR = joinpath(examples_dir(), "euler", "dry_air", "global_circulation")

    # Short simulation of the baroclinic instability on a coarse mesh of a spherical shell
    # with a thickness of 30 km and `trees_per_cube_face[2]` layers of elements
    function run_baroclinic_instability(trees_per_cube_face)
        trixi_include(@__MODULE__,
                      joinpath(EXAMPLES_DIR,
                               "elixir_potential_temperature_baroclinic_instability.jl"),
                      trees_per_cube_face = trees_per_cube_face, polydeg = 3,
                      tspan = (0.0, 100.0),
                      # Avoid intermediate time stops of the `SaveSolutionCallback`
                      dt = 100000)
    end
end

@testitem "Visualization: spherical_layer_data" setup=[Setup, Visualization] tags=[:visualization] begin
    run_baroclinic_instability((2, 2))
    n_layer_elements = 6 * 2^2
    nvisnodes = 5

    geopotential(u, equations) = u[6]
    data = spherical_layer_data(geopotential, sol; nvisnodes)
    @test size(data.lon) == size(data.lat) == size(data.values) ==
          (nvisnodes^2, n_layer_elements)
    @test all(lon -> 0 <= lon < 360, data.lon)
    @test all(lat -> -90 <= lat <= 90, data.lat)
    @test spherical_layer_data(geopotential, sol; layer = 1, nvisnodes) == data
    @test size(spherical_layer_data(geopotential, sol).lon) ==
          ((2 * polydeg)^2, n_layer_elements)

    # The geopotential is constant on spheres. Thus, it shows whether the correct layer is
    # chosen. The second layer of elements starts at half the thickness of 30 km.
    g = EARTH_GRAVITATIONAL_ACCELERATION
    @test all(isapprox(-EARTH_RADIUS * g), data.values)
    data_2 = spherical_layer_data(geopotential, sol; layer = 2, nvisnodes)
    @test all(isapprox(-EARTH_RADIUS^2 * g / (EARTH_RADIUS + 15000)), data_2.values)
    # The second layer lies above the first one
    @test sind.(data_2.lon) ≈ sind.(data.lon)
    @test cosd.(data_2.lon) ≈ cosd.(data.lon)
    @test data_2.lat ≈ data.lat

    @test spherical_layer_data(pressure, sol.u[end], semi; layer = 2) ==
          spherical_layer_data(pressure, sol; layer = 2)

    # Check that the values belong to the longitude and latitude by comparing with the
    # initial condition. At the bottom of the second layer, the velocity of the initial
    # condition has no perturbation, and its Cartesian components depend on both the
    # longitude and the latitude. With wrong coordinates (e.g., shifted by 30° in longitude
    # and mirrored at the equator), the relative error of the momentum is about 0.5.
    for variable in (2, 3)
        data_u0 = spherical_layer_data((u, equations) -> u[variable], ode.u0, semi;
                                       layer = 2, nvisnodes)
        values_exact = map(data_u0.lon, data_u0.lat) do lon, lat
            x = (EARTH_RADIUS + 15000) *
                SVector(cosd(lat) * cosd(lon), cosd(lat) * sind(lon), sind(lat))
            return initial_condition_baroclinic_instability(x, 0.0, equations)[variable]
        end
        # The deviation results from the interpolation of the polynomials of degree 3
        @test isapprox(data_u0.values, values_exact; rtol = 1.0e-2)
    end

    @test_throws ArgumentError spherical_layer_data(pressure, sol; layer = :top)
    @test_throws ArgumentError spherical_layer_data(pressure, sol; layer = 0)
    @test_throws ArgumentError spherical_layer_data(pressure, sol; layer = 1.5)
    # There are only two layers of elements
    @test_throws ArgumentError spherical_layer_data(pressure, sol; layer = 3)
end

@testitem "Visualization: spherical_layer_data errors" setup=[Setup] tags=[:visualization] begin
    # A periodic mesh of a cube has no boundary `:inside` at the bottom of a spherical shell
    equations = CompressibleEulerEquations3D(1.4)
    solver = DGSEM(polydeg = 1, surface_flux = flux_lax_friedrichs)
    mesh = P4estMesh((1, 1, 1); polydeg = 1, coordinates_min = (0.0, 0.0, 0.0),
                     coordinates_max = (1.0, 1.0, 1.0))
    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_constant,
                                        solver;
                                        boundary_conditions = boundary_condition_periodic)
    u_ode = compute_coefficients(0.0, semi)
    @test_throws ArgumentError spherical_layer_data(pressure, u_ode, semi)
end

@testitem "Visualization: triangulation in longitude-latitude" setup=[
    Setup,
    Visualization
] tags=[:visualization] begin
    # Signed area of the triangle with vertices `a`, `b`, `c`
    area(a, b, c) = ((b[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (b[2] - a[2])) / 2

    vertices(triangulation, triangle) = (SVector(triangulation.lon[vertex],
                                                 triangulation.lat[vertex])
                                         for vertex in triangle)

    function is_covered(triangulation, point)
        return any(eachcol(triangulation.triangles)) do triangle
            a, b, c = vertices(triangulation, triangle)
            areas = (area(point, a, b), area(point, b, c), area(point, c, a))
            return all(>=(-1.0e-12), areas) || all(<=(1.0e-12), areas)
        end
    end

    # Points in the interior of cells of 5° x 5° in longitude and latitude
    sample_points = [SVector(lon, lat) for lon in 2.5:5:357.5, lat in -87.5:5:87.5]

    function test_triangulation(data)
        triangulation = TrixiAtmo._spherical_layer_triangulation(data)
        (; lon, lat, values, triangles) = triangulation
        @test length(lon) == length(lat) == length(values)
        @test size(triangles, 1) == 3
        @test all(in(eachindex(lon)), triangles)
        # The triangulation covers the whole range of longitudes
        @test minimum(lon[triangles]) <= 0
        @test maximum(lon[triangles]) >= 360
        # No triangle spans more than 180° in longitude
        @test all(t -> maximum(lon[t]) - minimum(lon[t]) <= 180, eachcol(triangles))
        return triangulation
    end

    # With an even number of trees per cube face dimension, the poles are vertices of
    # elements and the periodic boundary at 0°/360° longitude lies at element interfaces.
    # Thus, the triangulation covers exactly the rectangle of all longitudes and latitudes.
    run_baroclinic_instability((2, 2))
    triangulation = test_triangulation(spherical_layer_data(pressure, sol; nvisnodes = 5))
    @test extrema(triangulation.lon) == (0, 360)
    @test extrema(triangulation.lat) == (-90, 90)
    total_area = sum(eachcol(triangulation.triangles)) do triangle
        return abs(area(vertices(triangulation, triangle)...))
    end
    @test total_area ≈ 360 * 180
    @test all(point -> is_covered(triangulation, point), sample_points)

    # With an odd number, the poles lie inside of elements and some elements cross the
    # periodic boundary at 0°/360° longitude. These elements are added a second time.
    run_baroclinic_instability((3, 2))
    data = spherical_layer_data(pressure, sol; nvisnodes = 5)
    triangulation = test_triangulation(data)
    @test minimum(triangulation.lon) < 0
    @test maximum(triangulation.lon) > 360
    # The poles are visualization points
    @test all(point -> is_covered(triangulation, point), sample_points)

    # The poles lie inside of cells of visualization points, which are removed
    triangulation = test_triangulation(spherical_layer_data(pressure, sol; nvisnodes = 4))
    @test !all(point -> is_covered(triangulation, point), sample_points)
    @test all(point -> is_covered(triangulation, point),
              filter(point -> abs(point[2]) < 80, sample_points))

    # A range of longitudes including 0° (from 60°W to 300°E), where the periodic boundary
    # of the longitude lies at its lower limit inside of elements
    lon_min = -60
    triangulation = TrixiAtmo._spherical_layer_triangulation(data; lon_min)
    (; lon, triangles) = triangulation
    @test minimum(lon[triangles]) <= lon_min
    @test maximum(lon[triangles]) >= lon_min + 360
    @test all(t -> maximum(lon[t]) - minimum(lon[t]) <= 180, eachcol(triangles))
    @test all(point -> is_covered(triangulation, point),
              [point + SVector(lon_min, 0) for point in sample_points])
end

@testitem "Visualization: plot_spherical_layer" setup=[Setup, Visualization] tags=[:visualization] begin
    using CairoMakie: CairoMakie, Makie

    run_baroclinic_instability((2, 2))

    pressure_hpa(u, equations) = pressure(u, equations) / 100
    fap = plot_spherical_layer(pressure_hpa, sol;
                               lon_limits = (0, 240), lat_limits = (0, 90),
                               levels = 5, colormap = :plasma,
                               title = "Surface Pressure [hPa]")
    @test fap isa Makie.FigureAxisPlot
    fig, ax, plt = fap
    @test fig isa Makie.Figure
    @test ax isa Makie.Axis
    @test plt isa Makie.Tricontourf
    @test length(plt.levels[]) == 6
    mktempdir() do dir
        filename = joinpath(dir, "surface_pressure.png")
        CairoMakie.save(filename, fap)
        @test isfile(filename)
    end

    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1])
    data = spherical_layer_data(pressure_hpa, sol; layer = 2)
    levels = range(extrema(data.values)...; length = 6)
    # Additional keywords are passed to `Makie.tricontourf!`
    plt = plot_spherical_layer!(ax, pressure_hpa, sol.u[end], semi; layer = 2,
                                nvisnodes = 6, levels, extendhigh = :auto)
    @test plt isa Makie.Tricontourf
    @test plt.extendhigh[] === :auto

    ext = Base.get_extension(TrixiAtmo, :TrixiAtmoMakieExt)

    # The default levels only use the values at the vertices of triangles overlapping the
    # plotted range. Here, the second triangle lies in the southern hemisphere.
    lon = [10.0, 20.0, 10.0, 200.0, 250.0, 200.0]
    lat = [10.0, 10.0, 20.0, -50.0, -50.0, -40.0]
    values = [1.0, 2.0, 3.0, -100.0, 100.0, 50.0]
    triangles = [1 4; 2 5; 3 6]
    @test ext.default_levels(5, lon, lat, values, triangles, (0, 240), (0, 90)) ≈
          range(1, 3; length = 6)
    # Triangles crossing the boundary of the plotted range are used as well
    @test ext.default_levels(5, lon, lat, values, triangles, (0, 240), (-45, 90)) ≈
          range(-100, 100; length = 6)
    # If no triangle overlaps the plotted range, all values are used
    @test ext.default_levels(2, lon, lat, values, triangles, (100, 180), (0, 90)) ≈
          range(-100, 100; length = 3)

    # Ticks and their labels
    @test ext.degree_ticks((0, 240)) == 0:30:240
    @test ext.degree_ticks((-90, 90)) == -90:30:90
    @test ext.degree_ticks((-60, 180)) == -60:30:180
    @test ext.degree_ticks((10, 20)) === Makie.automatic
    @test ext.format_longitude.([-60, 0, 30, 180, 240, 360, 22.5]) ==
          ["60W", "0", "30E", "180", "120W", "0", "22.5E"]
    @test ext.format_latitude.([-30, 0, 60, 90]) == ["30S", "0", "60N", "90N"]
end
