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

@testitem "Visualization: PlotDataShell" setup=[Setup, Visualization] tags=[:visualization] begin
    run_baroclinic_instability((2, 2))
    nvisnodes = 5

    # Values of a single variable
    series_values(pds) = pds.plot_data.data[:, pds.variable_id]

    pd = PlotDataShell(sol; solution_variables = cons2cons, nvisnodes)
    @test keys(pd) == ("rho", "rho_v1", "rho_v2", "rho_v3", "rho_theta", "phi")
    @test length(pd) == 6
    @test length(pd.x) == length(pd.y) == size(pd.data, 1)
    @test size(pd.data, 2) == 6
    @test size(pd.triangles, 1) == 3
    @test all(in(eachindex(pd.x)), pd.triangles)
    @test pd.projection == LonLat(lon_min = 0)
    # With an even number of trees per cube face dimension, the triangulation covers exactly
    # the rectangle of all longitudes and latitudes (see also the test of the triangulation)
    @test extrema(pd.x) == (0, 360)
    @test extrema(pd.y) == (-90, 90)
    @test pd["phi"] isa TrixiAtmo.Trixi.PlotDataSeries
    @test pd["phi"].variable_id == 6
    @test_throws KeyError pd["p"]
    @test [name for (name, _) in pd] == collect(keys(pd))
    @test startswith(repr(pd), "PlotDataShell{")

    @test PlotDataShell(sol.u[end], semi; solution_variables = cons2cons, layer = 1,
                        nvisnodes).data == pd.data
    # By default, the primitive variables are visualized with `nvisnodes = 2 * polydeg`
    pd_prim = PlotDataShell(sol)
    @test keys(pd_prim) == ("rho", "v1", "v2", "v3", "p1", "phi")
    @test pd_prim.data ==
          PlotDataShell(sol; solution_variables = cons2prim,
                        nvisnodes = 2 * polydeg).data
    @test length(pd_prim.x) > length(pd.x)

    # Scalar quantities use the same triangulation as all other variables
    pd_pressure = PlotDataShell(sol; solution_variables = pressure, nvisnodes)
    @test keys(pd_pressure) == ("pressure",)
    @test pd_pressure.x == pd.x
    @test pd_pressure.y == pd.y
    @test pd_pressure.triangles == pd.triangles
    @test vec(pd_pressure.data) ≈ series_values(PlotDataShell(sol; nvisnodes)["p1"])
    pd_pressure = PlotDataShell(sol; solution_variables = pressure, nvisnodes,
                                variable_names = "p [Pa]")
    @test keys(pd_pressure) == ("p [Pa]",)
    two_variables(u, equations) = SVector(u[1], u[5])
    @test_throws ArgumentError PlotDataShell(sol; solution_variables = two_variables)
    pd_two = PlotDataShell(sol; solution_variables = two_variables,
                           variable_names = ("rho", "rho_theta"), nvisnodes)
    @test pd_two.data == pd.data[:, [1, 5]]
    @test_throws ArgumentError PlotDataShell(sol; variable_names = ("a", "b"))

    # The geopotential is constant on spheres. Thus, it shows whether the correct layer is
    # chosen. The second layer of elements starts at half the thickness of 30 km.
    g = EARTH_GRAVITATIONAL_ACCELERATION
    @test all(isapprox(-EARTH_RADIUS * g), series_values(pd["phi"]))
    @test PlotDataShell(sol; solution_variables = cons2cons, layer = :surface,
                        nvisnodes).data == pd.data
    pd_2 = PlotDataShell(sol; solution_variables = cons2cons, layer = 2, nvisnodes)
    @test all(isapprox(-EARTH_RADIUS^2 * g / (EARTH_RADIUS + 15000)),
              series_values(pd_2["phi"]))

    # The `n`-th element of the second layer lies above the `n`-th element of the first layer
    u = TrixiAtmo.Trixi.wrap_array(sol.u[end], semi)
    function cartesian_coordinates(layer)
        coordinates, _ = TrixiAtmo._layer_visualization_points(u, equations, solver,
                                                               semi.cache, cons2cons,
                                                               layer, nvisnodes)
        return coordinates
    end
    for (x_1, x_2) in zip(cartesian_coordinates(1), cartesian_coordinates(2))
        @test x_2 ≈ (EARTH_RADIUS + 15000) / EARTH_RADIUS * x_1
    end

    # Check that the values belong to the longitude and latitude by comparing with the
    # initial condition. At the bottom of the second layer, the velocity of the initial
    # condition has no perturbation, and its Cartesian components depend on both the
    # longitude and the latitude.
    pd_u0 = PlotDataShell(ode.u0, semi; solution_variables = cons2cons, layer = 2,
                          nvisnodes)
    for (variable, name) in ((2, "rho_v1"), (3, "rho_v2"))
        values_exact = map(pd_u0.x, pd_u0.y) do lon, lat
            x = (EARTH_RADIUS + 15000) *
                SVector(cosd(lat) * cosd(lon), cosd(lat) * sind(lon), sind(lat))
            return initial_condition_baroclinic_instability(x, 0.0, equations)[variable]
        end
        # The deviation results from the interpolation of the polynomials of degree 3
        @test isapprox(series_values(pd_u0[name]), values_exact; rtol = 1.0e-2)
    end

    @test_throws ArgumentError PlotDataShell(sol; layer = :top)
    @test_throws ArgumentError PlotDataShell(sol; layer = 0)
    @test_throws ArgumentError PlotDataShell(sol; layer = 1.5)
    # There are only two layers of elements
    @test_throws ArgumentError PlotDataShell(sol; layer = 3)

    # Swapping the primary and secondary sides of all interfaces describes the same mesh.
    # Thus, the layers must not change. In the mesh above, the lower element is the
    # primary element of each interface between two layers, so this checks that the upper
    # neighbors are found if the upper element is the primary element.
    # Since this modifies the cache, it must be the last test of this test item.
    (; interfaces) = semi.cache
    interfaces.neighbor_ids .= interfaces.neighbor_ids[[2, 1], :]
    interfaces.node_indices .= interfaces.node_indices[[2, 1], :]
    @test PlotDataShell(sol; solution_variables = cons2cons, layer = 2, nvisnodes).data ==
          pd_2.data
end

@testitem "Visualization: PlotDataShell errors" setup=[Setup] tags=[:visualization] begin
    # A periodic mesh of a cube has no boundary `:inside` at the bottom of a spherical shell
    equations = CompressibleEulerEquations3D(1.4)
    solver = DGSEM(polydeg = 1, surface_flux = flux_lax_friedrichs)
    mesh = P4estMesh((1, 1, 1); polydeg = 1, coordinates_min = (0.0, 0.0, 0.0),
                     coordinates_max = (1.0, 1.0, 1.0))
    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_constant,
                                        solver;
                                        boundary_conditions = boundary_condition_periodic)
    u_ode = compute_coefficients(0.0, semi)
    @test_throws ArgumentError PlotDataShell(u_ode, semi)

    # A non-periodic mesh of a cube whose boundaries are all renamed to `:inside` has the
    # boundary `:inside` also at faces other than the bottom faces of the elements
    mesh = P4estMesh((1, 1, 1); polydeg = 1, coordinates_min = (0.0, 0.0, 0.0),
                     coordinates_max = (1.0, 1.0, 1.0), periodicity = false)
    boundary_condition = BoundaryConditionDirichlet(initial_condition_constant)
    boundary_conditions = Dict(name => boundary_condition
                               for name in (:x_neg, :x_pos, :y_neg, :y_pos, :z_neg, :z_pos))
    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_constant,
                                        solver; boundary_conditions)
    semi.cache.boundaries.name .= :inside
    u_ode = compute_coefficients(0.0, semi)
    @test_throws ArgumentError("expected the `:inside` boundary at the bottom face " *
                               "of each element") PlotDataShell(u_ode, semi)
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

    # Longitude and latitude of the visualization points of each element at the surface
    function visualization_points(nvisnodes)
        u = TrixiAtmo.Trixi.wrap_array(sol.u[end], semi)
        coordinates, _ = TrixiAtmo._layer_visualization_points(u, equations, solver,
                                                               semi.cache, pressure, 1,
                                                               nvisnodes)
        return TrixiAtmo._cartesian_to_lon_lat(coordinates...)
    end

    function test_triangulation(lon_vis, lat_vis; lon_min = 0)
        @test size(lon_vis) == size(lat_vis)
        # Points on the meridian at 0° with a tiny negative `y` coordinate (due to round-off
        # errors) have a longitude of `360` instead of `0`
        @test all(lon -> 0 <= lon <= 360, lon_vis)
        @test all(lat -> -90 <= lat <= 90, lat_vis)

        triangulation = TrixiAtmo._spherical_layer_triangulation(lon_vis, lat_vis;
                                                                 lon_min)
        (; lon, lat, triangles, point_indices) = triangulation
        @test length(lon) == length(lat) == length(point_indices)
        @test size(triangles, 1) == 3
        @test all(in(eachindex(lon)), triangles)
        # The triangulation covers the whole range of longitudes
        @test minimum(lon[triangles]) <= lon_min
        @test maximum(lon[triangles]) >= lon_min + 360
        # No triangle spans more than 180° in longitude
        @test all(t -> maximum(lon[t]) - minimum(lon[t]) <= 180, eachcol(triangles))

        # Each point of the triangulation is an original visualization point, possibly with
        # the longitude shifted by a multiple of 360°. Only the longitude at a pole is
        # arbitrary, so the values of all quantities are given by the original points.
        @test lat == lat_vis[point_indices]
        not_at_pole = findall(lat -> !isapprox(abs(lat), 90), lat)
        lon_original = lon_vis[point_indices[not_at_pole]]
        @test sind.(lon[not_at_pole]) ≈ sind.(lon_original)
        @test cosd.(lon[not_at_pole]) ≈ cosd.(lon_original)
        return triangulation
    end

    # With an even number of trees per cube face dimension, the poles are vertices of
    # elements and the periodic boundary at 0°/360° longitude lies at element interfaces.
    # Thus, the triangulation covers exactly the rectangle of all longitudes and latitudes.
    run_baroclinic_instability((2, 2))
    lon_vis, lat_vis = visualization_points(5)
    n_layer_elements = 6 * 2^2
    @test size(lon_vis) == (5^2, n_layer_elements)
    triangulation = test_triangulation(lon_vis, lat_vis)
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
    lon_vis, lat_vis = visualization_points(5)
    triangulation = test_triangulation(lon_vis, lat_vis)
    @test minimum(triangulation.lon) < 0
    @test maximum(triangulation.lon) > 360
    # The poles are visualization points
    @test all(point -> is_covered(triangulation, point), sample_points)

    # The poles lie inside of cells of visualization points, which are removed
    triangulation = test_triangulation(visualization_points(4)...)
    @test !all(point -> is_covered(triangulation, point), sample_points)
    @test all(point -> is_covered(triangulation, point),
              filter(point -> abs(point[2]) < 80, sample_points))

    # A range of longitudes including 0° (from 60°W to 300°E), where the periodic boundary
    # of the longitude lies at its lower limit inside of elements
    lon_min = -60
    triangulation = test_triangulation(lon_vis, lat_vis; lon_min)
    @test all(point -> is_covered(triangulation, point),
              [point + SVector(lon_min, 0) for point in sample_points])

    # The `PlotDataShell` uses this triangulation
    pd = PlotDataShell(sol; nvisnodes = 5, projection = LonLat(; lon_min))
    @test pd.x == triangulation.lon
    @test pd.y == triangulation.lat
    @test pd.triangles == triangulation.triangles
end

@testitem "Visualization: Makie.contourf" setup=[Setup, Visualization] tags=[:visualization] begin
    using CairoMakie: CairoMakie, Makie

    run_baroclinic_instability((2, 2))
    ext = Base.get_extension(TrixiAtmo, :TrixiAtmoMakieExt)

    # Contour lines drawn by `Makie.contourf!` in the axis `ax`
    contour_lines(ax) = only(filter(plot -> plot isa Makie.Tricontour, ax.scene.plots))

    # Several variables of the same `PlotDataShell`
    pd = PlotDataShell(sol)
    fap = Makie.contourf(pd["p1"]; lon_limits = (0, 240), lat_limits = (0, 90),
                         levels = 5, colormap = :plasma)
    @test fap isa Makie.FigureAxisPlot
    fig, ax, plt = fap
    @test fig isa Makie.Figure
    @test ax isa Makie.Axis
    @test plt isa Makie.Tricontourf
    # An integer `levels` is passed to Makie.jl, which computes the levels
    @test plt.levels[] == 5
    @test length(plt.computed_levels[]) == 6
    # The contour lines are drawn at the edges of the bands of the filled contours
    @test contour_lines(ax).levels[] == plt.computed_levels[]
    # The title is the name of the variable by default
    @test ax.title[] == "p1"
    @test ax.limits[] == ((0, 240), (0, 90))
    mktempdir() do dir
        filename = joinpath(dir, "surface_pressure.png")
        CairoMakie.save(filename, fap)
        @test isfile(filename)
    end
    # By default, the whole range of the projection is plotted
    fig, ax, plt = Makie.contourf(pd["v1"]; title = "v1 [m/s]")
    @test plt isa Makie.Tricontourf
    @test ax.title[] == "v1 [m/s]"
    @test ax.limits[] == ((0, 360), (-90, 90))
    levels = plt.computed_levels[]
    @test first(levels) ≈ minimum(pd.data[:, 2])
    @test last(levels) ≈ maximum(pd.data[:, 2])

    # Makie.jl computes the levels only from the values of the triangles overlapping the
    # plotted range. The surface pressure decreases towards the poles, so that its minimum
    # between 0° and 30°N is larger than on the whole sphere.
    fig, ax, plt = Makie.contourf(pd["p1"]; lat_limits = (0, 30))
    _, _, values, _ = ext.crop_triangulation(pd.x, pd.y, pd.data[:, 5], pd.triangles,
                                             (0, 360), (0, 30))
    levels = plt.computed_levels[]
    @test first(levels) ≈ minimum(values)
    @test last(levels) ≈ maximum(values)
    @test !(first(levels) ≈ minimum(pd.data[:, 5]))

    # With a `colorscale`, the levels of the contour lines are transformed back as in the
    # colorbar of Makie.jl
    fig, ax, plt = Makie.contourf(pd["p1"]; levels = 5, colorscale = log10)
    @test contour_lines(ax).levels[] ≈ exp10.(plt.computed_levels[])

    # The plotted range of longitude must lie in the range of the projection
    @test_throws ArgumentError Makie.contourf(pd["p1"]; lon_limits = (-60, 180))
    @test_throws ArgumentError Makie.contourf(pd["p1"]; lon_limits = (0, 400))
    pd = PlotDataShell(sol; projection = LonLat(lon_min = -60))
    fig, ax, plt = Makie.contourf(pd["p1"]; lon_limits = (-60, 180))
    @test ax.limits[] == ((-60, 180), (-90, 90))
    fig, ax, plt = Makie.contourf(pd["p1"])
    @test ax.limits[] == ((-60, 300), (-90, 90))

    # A `PlotDataShell` with a single variable can be plotted directly
    pressure_hpa(u, equations) = pressure(u, equations) / 100
    pd = PlotDataShell(sol; solution_variables = pressure_hpa,
                       projection = LonLat(lon_min = -60))
    fig, ax, plt = Makie.contourf(pd; lon_limits = (-60, 180), lat_limits = (0, 90),
                                  levels = 5, colormap = :plasma)
    @test plt isa Makie.Tricontourf
    @test length(plt.computed_levels[]) == 6
    @test ax.title[] == "pressure_hpa"
    @test ax.limits[] == ((-60, 180), (0, 90))
    fig, ax, plt = Makie.contourf(pd; title = "Surface Pressure [hPa]")
    @test ax.title[] == "Surface Pressure [hPa]"
    @test_throws ArgumentError Makie.contourf(PlotDataShell(sol))
    @test_throws ArgumentError Makie.contourf!(ax, PlotDataShell(sol))

    # For the nearly constant surface pressure of the initial condition, Makie.jl puts all
    # values into a single band instead of showing round-off errors, also for the default
    # (even) number of bands
    pd = PlotDataShell(ode.u0, semi; solution_variables = pressure_hpa)
    values = vec(pd.data)
    @test isapprox(minimum(values), maximum(values))
    fig, ax, plt = Makie.contourf(pd)
    @test !any(level -> minimum(values) <= level <= maximum(values),
               plt.computed_levels[])

    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1])
    pd = PlotDataShell(sol; solution_variables = pressure_hpa, variable_names = "p",
                       layer = 2, nvisnodes = 6)
    levels = range(extrema(pd.data)...; length = 6)
    # Additional keywords are passed to `Makie.tricontourf!`
    plt = Makie.contourf!(ax, pd["p"]; levels, extendhigh = :auto)
    @test plt isa Makie.Tricontourf
    @test plt.extendhigh[] === :auto
    @test plt.computed_levels[] ≈ levels
    plt = Makie.contourf!(ax, pd; levels, extendhigh = :auto)
    @test plt isa Makie.Tricontourf
    @test plt.extendhigh[] === :auto

    # Only the triangles overlapping the plotted range are passed to Makie.jl, with their
    # vertices numbered consecutively. Here, the first triangle lies in the southern
    # hemisphere.
    lon = [200.0, 250.0, 200.0, 10.0, 20.0, 10.0]
    lat = [-50.0, -50.0, -40.0, 10.0, 10.0, 20.0]
    values = [-100.0, 100.0, 50.0, 1.0, 2.0, 3.0]
    triangles = [1 6; 2 4; 3 5]
    @test ext.crop_triangulation(lon, lat, values, triangles, (0, 240), (0, 90)) ==
          ([10.0, 20.0, 10.0], [10.0, 10.0, 20.0], [1.0, 2.0, 3.0],
           reshape([3, 1, 2], 3, 1))
    # Triangles crossing the boundary of the plotted range are kept as well
    @test ext.crop_triangulation(lon, lat, values, triangles, (0, 240), (-45, 90)) ==
          (lon, lat, values, triangles)
    # If no triangle overlaps the plotted range, all triangles are kept
    @test ext.crop_triangulation(lon, lat, values, triangles, (100, 180), (0, 90)) ==
          (lon, lat, values, triangles)

    # Ticks and their labels
    @test ext.degree_ticks((0, 240)) == 0:30:240
    @test ext.degree_ticks((-90, 90)) == -90:30:90
    @test ext.degree_ticks((-60, 180)) == -60:30:180
    @test ext.degree_ticks((10, 20)) === Makie.automatic
    @test ext.format_longitude.([-60, 0, 30, 180, 240, 360, 22.5]) ==
          ["60W", "0", "30E", "180", "120W", "0", "22.5E"]
    @test ext.format_latitude.([-30, 0, 60, 90]) == ["30S", "0", "60N", "90N"]
end
