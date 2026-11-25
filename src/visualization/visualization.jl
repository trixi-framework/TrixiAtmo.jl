# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Abstract supertype of the map projections of a `PlotDataShell`, which determine the
# coordinates of the visualization points in the plot and their triangulation. Further
# projections (e.g., other maps or a three-dimensional globe) can be added as subtypes with a
# method of `_triangulate`.
abstract type AbstractProjection end

"""
    LonLat(; lon_min = 0)

Equirectangular map projection for [`PlotDataShell`](@ref), which uses the longitude and
latitude in degrees as horizontal and vertical coordinates. The periodic boundary of the
longitude is placed at `lon_min`, i.e., the longitude ranges from `lon_min` to
`lon_min + 360`. For example, `LonLat(lon_min = -180)` represents the longitudes from 180°W
to 180°E, with 0° in the center.
"""
struct LonLat{RealT <: Real} <: AbstractProjection
    lon_min::RealT
end

LonLat(; lon_min = 0) = LonLat(lon_min)

"""
    PlotDataShell(u_ode, semi; solution_variables = nothing, variable_names = nothing,
                  layer = :surface, nvisnodes = nothing, projection = LonLat())
    PlotDataShell(sol; kwargs...)

Create a `PlotDataShell` object for visualizing a horizontal layer of a simulation on a
spherical shell, e.g., on a mesh created by [`P4estMeshCubedSphereTopography`](@ref) or
`Trixi.P4estMeshCubedSphere`, in the map `projection`. The bottom of the spherical shell must
be the boundary `:inside` at the bottom face (third local coordinate `-1`) of the elements.
For meshes with topography, the layers of elements and thus the visualized layers follow the
terrain.

All variables given by `solution_variables` are evaluated at once, and the visualization
points are triangulated only once in the map projection, so that all variables can be plotted
without repeating these steps. Single variables can be extracted via `pd[variable_name]` and
plotted with `Makie.contourf` and `Makie.contourf!` (see "Plotting" below). As for the plot
data types of Trixi.jl, `keys(pd)` returns the variable names, and `pd` can be iterated over
pairs of the variable names and the single variables.

For simulations on GPUs, only the data at the nodes of the visualized layer is copied to the
CPU. Parallel simulations with MPI are not supported.

# Keywords
- `solution_variables`: function `(u, equations) -> SVector(...)` converting the conservative
  variables `u` at a node to the variables to be visualized, which may also return a scalar
  (e.g., `pressure`). By default, the primitive variables (if `cons2prim` is available) or the
  conservative variables (otherwise) are used.
- `variable_names`: names of the visualized variables. By default, they are given by
  `varnames(solution_variables, equations)`. If this method does not exist, the name of the
  function `solution_variables` is used for a scalar quantity, e.g., `"pressure"`.
- `layer`: the horizontal layer to extract, either `:surface` for the bottom of the spherical
  shell or a positive integer `n` for the bottom face of the `n`-th layer of elements counted
  from the bottom, i.e., `layer = 1` is equivalent to `layer = :surface`.
- `nvisnodes`: number of visualization nodes per element and horizontal direction; defaults to
  twice the polynomial degree of the solver as for 2D plots in Trixi.jl.
- `projection`: map projection of the visualization points; currently only [`LonLat`](@ref).

# Plotting
With Makie.jl loaded, a single variable `pds = pd[variable_name]` with the projection
[`LonLat`](@ref) can be plotted as filled contours with contour lines in longitude-latitude
coordinates via
```julia
fig, ax, plt = Makie.contourf(pds; kwargs...)
plt = Makie.contourf!(ax, pds; kwargs...)
```
`Makie.contourf` creates a figure with an axis and a colorbar and returns a
`Makie.FigureAxisPlot`, which can be destructured as above. `Makie.contourf!` adds the plot to
the existing axis `ax` and returns the filled contour plot, e.g., to create a
`Makie.Colorbar`. A `PlotDataShell` with a single variable, e.g., created with
`solution_variables = pressure`, can also be passed directly instead of `pds`.

The visualization points of each element are triangulated separately in the
longitude-latitude plane. If a pole lies inside a cell of visualization points (e.g., for
an odd number of trees per cube face dimension and an even `nvisnodes`), this cell cannot
be represented and a small region around the pole remains empty.

Keywords for plotting:
- `lon_limits`, `lat_limits`: plotted range of longitude and latitude in degrees. The range of
  longitude must lie in the range from `lon_min` to `lon_min + 360` of the projection
  `LonLat(; lon_min)` and defaults to this range. Thus, a range of longitude including 0°,
  e.g., `lon_limits = (-60, 180)` for the range from 60°W to 180°, requires a projection with
  `lon_min <= -60`, e.g., `LonLat(lon_min = -60)` or `LonLat(lon_min = -180)`. The range of
  latitude defaults to `(-90, 90)`.
- `title`: title of the plot (only for `Makie.contourf`); defaults to the name of the variable.
- `colormap`: colormap of the filled contours; defaults to `:viridis`.
- `linecolor`, `linewidth`, `linestyle`: color, width, and style (e.g., `:dash`) of the
  contour lines at the edges of the bands of the filled contours; default to `:black`, `0.5`,
  and `:solid`. With `linewidth = 0`, no contour lines are drawn.

All other keywords are passed to `Makie.tricontourf!`, e.g., `levels` for the number of
equally spaced contour bands (10 by default) or a vector of contour levels, and
`extendlow = :auto` and `extendhigh = :auto` to also fill the regions with values outside of
the given `levels`. Only the triangles overlapping the plotted range of longitude and latitude
are passed to Makie.jl. Thus, the contour levels computed by Makie.jl for an integer `levels`
(or for `mode = :relative`) depend only on the values in the plotted range.

Other plotting functions of Makie.jl (e.g., `Makie.plot` or `Makie.contour`) and plots with
Plots.jl are not supported and throw an error.

# Examples
```julia
julia> pd = PlotDataShell(sol) # primitive variables at the surface

julia> keys(pd)
("rho", "v1", "v2", "v3", "p1", "phi")

julia> pd = PlotDataShell(sol; solution_variables = pressure, layer = 2,
                          projection = LonLat(lon_min = -180))
```
The surface pressure in hPa of the baroclinic instability (see
[`initial_condition_baroclinic_instability`](@ref)) on the northern hemisphere from 0° to
120°W can be plotted after a simulation with the elixir
`elixir_potential_temperature_baroclinic_instability.jl` in
`examples/euler/dry_air/global_circulation` via
```julia
using CairoMakie
pressure_hpa(u, equations) = pressure(u, equations) / 100
pd = PlotDataShell(sol; solution_variables = pressure_hpa,
                   variable_names = "Surface pressure [hPa]")
fig, ax, plt = contourf(pd; lon_limits = (0, 240), lat_limits = (0, 90),
                        levels = 11, colormap = :plasma)
save("surface_pressure.png", fig)
```
Several variables, e.g., the pressure and the first Cartesian component of the velocity at
the surface, can be plotted using the same triangulation via
```julia
pd = PlotDataShell(sol) # primitive variables
contourf(pd["p1"]; lon_limits = (0, 240), lat_limits = (0, 90))
contourf(pd["v1"]; lon_limits = (0, 240), lat_limits = (0, 90))
```
"""
struct PlotDataShell{Coordinates, Data, VariableNames, Projection} <:
       Trixi.AbstractPlotData{2}
    x::Coordinates # coordinates of the visualization points in the map projection,
    y::Coordinates # e.g., longitude and latitude in degrees
    triangles::Matrix{Int} # `3 × n_triangles` indices of the vertices of each triangle
    data::Data # `n_points × n_variables` values of the variables at the points
    variable_names::VariableNames
    projection::Projection
end

# Show only a truncated output for convenience
function Base.show(io::IO, pd::PlotDataShell)
    @nospecialize pd # reduce precompilation time

    print(io, "PlotDataShell{",
          typeof(pd.x), ", ",
          typeof(pd.data), ", ",
          typeof(pd.variable_names), ", ",
          typeof(pd.projection),
          "}(<x>, <y>, <triangles>, <data>, <variable_names>, <projection>)")
    return nothing
end

# Since `PlotDataShell <: Trixi.AbstractPlotData{2}`, the generic plotting methods of Trixi.jl
# for two-dimensional plot data (e.g., the Plots.jl recipes and `Makie.plot` as heatmap) would
# apply to it. However, they assume the data structure of `Trixi.PlotData2DCartesian` and
# fail with obscure errors. Thus, we throw a clear error instead (see also the method of
# `Makie.convert_arguments` in the Makie.jl extension).
function _throw_unsupported_plot()
    throw(ArgumentError("a `PlotDataShell` and its variables can only be plotted with " *
                        "`Makie.contourf` and `Makie.contourf!` (with Makie.jl loaded), " *
                        "e.g., `contourf(pd[variable_name])`"))
end

function RecipesBase.apply_recipe(::AbstractDict{Symbol, Any}, ::PlotDataShell)
    return _throw_unsupported_plot()
end

function RecipesBase.apply_recipe(::AbstractDict{Symbol, Any},
                                  ::Trixi.PlotDataSeries{<:PlotDataShell})
    return _throw_unsupported_plot()
end

function PlotDataShell(sol::Trixi.TrixiODESolution; kwargs...)
    return PlotDataShell(sol.u[end], sol.prob.p; kwargs...)
end

function PlotDataShell(u_ode, semi; solution_variables = nothing,
                       variable_names = nothing,
                       layer = :surface, nvisnodes = nothing, projection = LonLat())
    if layer === :surface
        layer = 1
    elseif !(layer isa Integer && layer >= 1)
        throw(ArgumentError("`layer` must be `:surface` or a positive integer, " *
                            "got `layer = $(repr(layer))`"))
    end

    if Trixi.mpi_isparallel()
        throw(ArgumentError("`PlotDataShell` does not support parallel simulations " *
                            "with MPI"))
    end

    u = Trixi.wrap_array(u_ode, semi)
    _, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)

    if nvisnodes === nothing
        nvisnodes = 2 * Trixi.polydeg(solver)
    end

    solution_variables_ = Trixi.digest_solution_variables(equations, solution_variables)

    # `u` wraps the memory of `u_ode` (see `Trixi.wrap_array`)
    coordinates, values = GC.@preserve u_ode _layer_visualization_points(u, equations,
                                                                         solver, cache,
                                                                         solution_variables_,
                                                                         layer,
                                                                         nvisnodes)
    n_variables = size(values, 3)
    variable_names_ = _variable_names(solution_variables_, equations, variable_names,
                                      n_variables)

    # The triangulation depends only on the coordinates of the visualization points.
    # Points that are duplicated in the triangulation (e.g., at the periodic boundary of the
    # longitude) get the values of the original points via `point_indices`.
    (; x, y, triangles, point_indices) = _triangulate(projection, coordinates...)
    data = reshape(values, :, n_variables)[point_indices, :]

    return PlotDataShell(x, y, triangles, data, variable_names_, projection)
end

# Return the Cartesian coordinates `(x, y, z)` and the values of `solution_variables` at the
# visualization points of the bottom faces of the elements of the `layer`-th layer of elements.
# The coordinates are matrices of size `(nvisnodes^2, n_layer_elements)`, and the values are
# an array of size `(nvisnodes^2, n_layer_elements, n_variables)`. The visualization points of
# each element are equally spaced in the reference element with the first index running
# fastest (see `Trixi.plotting_interpolation_matrix`).
function _layer_visualization_points(u, equations, solver, cache, solution_variables,
                                     layer, nvisnodes)
    (; node_coordinates) = cache.elements
    layer_elements = _layer_elements(layer, solver, cache)
    n_nodes_2d = Trixi.nnodes(solver)^2
    n_layer_elements = length(layer_elements)

    # Select the conservative variables and the coordinates at the nodes at the bottom face
    # of the elements of the layer. For simulations on GPUs, they are selected on the GPU,
    # and only these nodes are copied to the CPU.
    u_layer = Array(u[:, :, :, 1, layer_elements])
    layer_coordinates = reshape(Array(node_coordinates[:, :, :, 1, layer_elements]),
                                3, n_nodes_2d, n_layer_elements)

    # Evaluate `solution_variables` at these nodes
    nodes = Iterators.product(eachnode(solver), eachnode(solver), 1:n_layer_elements)
    layer_values = map(nodes) do (i, j, idx)
        u_node = Trixi.get_node_vars(u_layer, equations, solver, i, j, idx)
        return _as_svector(solution_variables(u_node, equations))
    end

    n_variables = length(first(layer_values))
    plotting_interp_matrix = Trixi.plotting_interpolation_matrix(solver; nvisnodes)

    # Interpolate the coordinates of the nodes at the bottom face of the elements of the layer
    coordinates = ntuple(i -> plotting_interp_matrix * layer_coordinates[i, :, :], 3)

    # Interpolate each variable separately
    values_2d = reshape(layer_values, n_nodes_2d, n_layer_elements)
    values = similar(coordinates[1], eltype(eltype(layer_values)),
                     size(coordinates[1])..., n_variables)
    for v in 1:n_variables
        values[:, :, v] = plotting_interp_matrix * getindex.(values_2d, v)
    end

    return coordinates, values
end

# Scalar quantities are treated as a single variable
_as_svector(x::Number) = SVector(x)
_as_svector(x::SVector) = x

function _variable_names(solution_variables, equations, variable_names, n_variables)
    if variable_names isa AbstractString
        variable_names = (variable_names,)
    elseif variable_names === nothing
        if hasmethod(varnames, Tuple{typeof(solution_variables), typeof(equations)})
            variable_names = varnames(solution_variables, equations)
        elseif n_variables == 1 && solution_variables isa Function
            # Use the name of the function, e.g., `"pressure"`, for a scalar quantity
            variable_names = (string(nameof(solution_variables)),)
        else
            throw(ArgumentError("`variable_names` must be passed since there is no " *
                                "method `varnames(solution_variables, equations)`"))
        end
    end

    if length(variable_names) != n_variables
        throw(ArgumentError("got $(length(variable_names)) `variable_names` for " *
                            "$n_variables variables"))
    end

    return SVector{n_variables}(variable_names...)
end

# Return the elements of the `layer`-th layer of elements counted from the bottom of the
# spherical shell. The `n`-th of these elements lies above the `n`-th element at the bottom.
function _layer_elements(layer, solver, cache)
    # Copy the connectivity to the CPU, since the containers are stored on the GPU for
    # simulations on GPUs
    boundary_names = Array(cache.boundaries.name)
    boundary_neighbor_ids = Array(cache.boundaries.neighbor_ids)
    boundary_node_indices = Array(cache.boundaries.node_indices)

    # Find the elements at the bottom of the spherical shell via the `:inside` boundary
    inside_boundaries = findall(==(:inside), boundary_names)
    if isempty(inside_boundaries)
        throw(ArgumentError("expected a mesh of a spherical shell with the boundary " *
                            "`:inside` at its bottom"))
    end
    layer_elements = boundary_neighbor_ids[inside_boundaries]

    # The `:inside` boundary must be the bottom face of each element, i.e., the face
    # with the first index in the third local coordinate direction
    expected_node_indices = (NodeIndex.IForward, NodeIndex.JForward, NodeIndex.Begin)
    if !all(boundary_node_indices[boundary] == expected_node_indices
            for boundary in inside_boundaries)
        throw(ArgumentError("expected the `:inside` boundary at the bottom face " *
                            "of each element"))
    end

    if layer == 1
        return layer_elements
    end

    # Find the element above each element, i.e., its neighbor across its top face.
    # At such an interface, the index in the third local coordinate direction is
    # `NodeIndex.End` in the lower element and `NodeIndex.Begin` in the upper element. At
    # all other interfaces of a spherical shell, these indices are of the form
    # `NodeIndex.IForward`, `NodeIndex.JBackward`, etc.
    interface_neighbor_ids = Array(cache.interfaces.neighbor_ids)
    interface_node_indices = Array(cache.interfaces.node_indices)
    upper_neighbors = zeros(Int, Trixi.nelements(solver, cache))
    for interface in Trixi.eachinterface(solver, cache)
        primary_element = interface_neighbor_ids[1, interface]
        secondary_element = interface_neighbor_ids[2, interface]
        primary_indices = interface_node_indices[1, interface]
        secondary_indices = interface_node_indices[2, interface]
        if primary_indices[3] == NodeIndex.End &&
           secondary_indices[3] == NodeIndex.Begin
            upper_neighbors[primary_element] = secondary_element
        elseif primary_indices[3] == NodeIndex.Begin &&
               secondary_indices[3] == NodeIndex.End
            upper_neighbors[secondary_element] = primary_element
        end
    end

    # Move upwards from the bottom layer by layer
    for _ in 2:layer
        for (idx, element) in enumerate(layer_elements)
            upper_element = upper_neighbors[element]
            if upper_element == 0
                throw(ArgumentError("`layer = $layer` exceeds the number of layers of " *
                                    "conforming elements"))
            end
            layer_elements[idx] = upper_element
        end
    end

    return layer_elements
end

# Return the coordinates `x` and `y` of the triangulated visualization points in the map
# `projection`, the `3 × n_triangles` matrix `triangles` of the indices of the vertices of
# each triangle, and the vector `point_indices` of the (linear) indices of the original
# visualization points of each point of the triangulation. The arguments are the Cartesian
# coordinates of the visualization points (see `_layer_visualization_points`).
function _triangulate(projection::LonLat, x, y, z)
    lon, lat = _cartesian_to_lon_lat(x, y, z)
    (; lon, lat, triangles, point_indices) = _spherical_layer_triangulation(lon, lat;
                                                                            projection.lon_min)
    return (; x = lon, y = lat, triangles, point_indices)
end

# Convert Cartesian coordinates to longitude in `[0, 360)` and latitude in degrees
function _cartesian_to_lon_lat(x, y, z)
    lon = similar(x)
    lat = similar(x)
    for n in eachindex(x, y, z)
        lon_n, lat_n, _ = cartesian_to_sphere(SVector(x[n], y[n], z[n]))
        lon[n] = rad2deg(lon_n)
        lat[n] = rad2deg(lat_n)
    end
    return lon, lat
end

# Triangulate the visualization points with longitude `lon_vis` and latitude `lat_vis` (in
# degrees, matrices of size `(nvisnodes^2, n_layer_elements)`) in the longitude-latitude
# plane for contour plots with `Makie.tricontourf`.
#
# As for contour plots of DG solutions in Trixi.jl, the visualization points of each element
# are triangulated separately using the triangulation of the visualization points on the
# reference element. In contrast to a Delaunay triangulation of all points, this does not
# connect points of different elements across discontinuities of the solution and avoids
# duplicate points, since the points at element interfaces belong to several elements.
# The longitude is shifted to the range from `lon_min` to `lon_min + 360°` such that the
# periodic boundary of the longitude is at the boundary of the plotted range starting at
# `lon_min`. Elements crossing this periodic boundary are made contiguous in longitude. If
# they extend beyond `lon_min + 360°`, they are added a second time shifted by -360° such
# that both ends of the range are covered.
# A pole is a line in longitude-latitude coordinates. Thus, a triangle with a vertex at a
# pole, where the longitude is not defined, is replaced by the quadrilateral between this
# line and the opposite edge of the triangle.
# In an element containing a pole in its interior, the longitude jumps by 360° somewhere.
# The longitude of the vertices of the triangles crossing this jump is shifted accordingly.
# Only a triangle enclosing a pole cannot be represented and is removed.
#
# The triangulation depends only on the coordinates. Since points are duplicated, the vector
# `point_indices` contains the linear index in `lon_vis` and `lat_vis` of the original
# visualization point of each point of the triangulation. Thus, the values of any quantity
# at the points of the triangulation are given by `values_vis[point_indices]`.
#
# Return a named tuple `(; lon, lat, triangles, point_indices)` of vectors of the
# coordinates of the points, a `3 × n_triangles` matrix of the indices of the vertices of
# each triangle, and the indices of the original visualization points.
function _spherical_layer_triangulation(lon_vis, lat_vis; lon_min = 0)
    # Triangulation of the equally spaced visualization points on the reference element
    # (see `Trixi.plotting_interpolation_matrix`), with the first index running fastest
    nvisnodes = isqrt(size(lon_vis, 1))
    reference_nodes = LinRange(-1, 1, nvisnodes)
    r = vec([reference_nodes[i] for i in 1:nvisnodes, j in 1:nvisnodes])
    s = vec([reference_nodes[j] for i in 1:nvisnodes, j in 1:nvisnodes])
    reference_triangles = Trixi.reference_plotting_triangulation((r, s))

    lon = eltype(lon_vis)[]
    lat = eltype(lat_vis)[]
    point_indices = Int[]
    triangles = Int[]
    linear_indices = LinearIndices(lon_vis)

    for element in axes(lon_vis, 2)
        lon_element = @. lon_min + mod(lon_vis[:, element] - lon_min, 360)
        if maximum(lon_element) - minimum(lon_element) > 180
            @. lon_element = ifelse(lon_element < lon_min + 180, lon_element + 360,
                                    lon_element)
        end
        shifts = maximum(lon_element) > lon_min + 360 ? (0, -360) : (0,)

        for shift in shifts
            offset = length(lon)
            append!(lon, lon_element .+ shift)
            append!(lat, view(lat_vis, :, element))
            append!(point_indices, view(linear_indices, :, element))

            for triangle in axes(reference_triangles, 1)
                vertices = ntuple(k -> offset + Int(reference_triangles[triangle, k]),
                                  3)
                _push_triangle!(triangles, lon, lat, point_indices, vertices)
            end
        end
    end

    return (; lon, lat, triangles = reshape(triangles, 3, :), point_indices)
end

function _push_triangle!(triangles, lon, lat, point_indices, vertices)
    at_pole = map(vertex -> isapprox(abs(lat[vertex]), 90), vertices)
    n_vertices_at_pole = count(at_pole)

    if n_vertices_at_pole == 0
        _push_triangle_contiguous!(triangles, lon, lat, point_indices, vertices)
    elseif n_vertices_at_pole == 1
        # Replace the triangle by the quadrilateral between the pole and the opposite edge
        # from `vertex_1` to `vertex_2`, using the longitudes of `vertex_1` and `vertex_2`
        # for two copies of the vertex at the pole
        k = findfirst(at_pole)
        pole = vertices[k]
        vertex_1 = vertices[mod1(k + 1, 3)]
        vertex_2 = vertices[mod1(k + 2, 3)]
        pole_1 = _push_point!(lon, lat, point_indices, lon[vertex_1], lat[pole],
                              point_indices[pole])
        pole_2 = _push_point!(lon, lat, point_indices, lon[vertex_2], lat[pole],
                              point_indices[pole])
        _push_triangle_contiguous!(triangles, lon, lat, point_indices,
                                   (pole_1, vertex_1, vertex_2))
        _push_triangle_contiguous!(triangles, lon, lat, point_indices,
                                   (pole_1, vertex_2, pole_2))
    end
    # Triangles with more than one vertex at a pole are degenerate and removed

    return triangles
end

function _push_triangle_contiguous!(triangles, lon, lat, point_indices, vertices)
    lon_min, lon_max = extrema(vertex -> lon[vertex], vertices)
    if lon_max - lon_min <= 180
        return append!(triangles, vertices)
    end

    # The triangle crosses the jump of the longitude by 360° in an element containing a
    # pole. Use copies of the vertices with the longitude shifted by multiples of 360° to
    # the range of 180° around the longitude of the first vertex.
    lon_reference = lon[first(vertices)]
    shifted_vertices = map(vertices) do vertex
        lon_shifted = lon_reference + mod(lon[vertex] - lon_reference + 180, 360) - 180
        return _push_point!(lon, lat, point_indices, lon_shifted, lat[vertex],
                            point_indices[vertex])
    end

    # Otherwise, the triangle encloses a pole and cannot be represented
    lon_min, lon_max = extrema(vertex -> lon[vertex], shifted_vertices)
    if lon_max - lon_min <= 180
        append!(triangles, shifted_vertices)
    end

    return triangles
end

function _push_point!(lon, lat, point_indices, lon_point, lat_point, point_index)
    push!(lon, lon_point)
    push!(lat, lat_point)
    push!(point_indices, point_index)
    return length(lon)
end
end # @muladd
