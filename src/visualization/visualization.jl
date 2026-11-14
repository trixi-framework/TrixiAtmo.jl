# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

"""
    spherical_layer_data(func, u_ode, semi; layer = :surface, nvisnodes = nothing)
    spherical_layer_data(func, sol; kwargs...)

Extract a scalar quantity on a horizontal layer of a simulation on a spherical shell
for visualization, e.g., on a mesh created by [`P4estMeshCubedSphereTopography`](@ref) or
`Trixi.P4estMeshCubedSphere`. The bottom of the spherical shell must be the boundary
`:inside` at the bottom face (third local coordinate `-1`) of the elements.

The function `func(u, equations)` takes the conservative variables `u` at a node and the
`equations` and returns the scalar quantity to be visualized, e.g., `pressure`.

Return a named tuple `(; lon, lat, values)` of matrices of size
`(nvisnodes^2, n_layer_elements)` with the longitude and latitude (in degrees) of the
visualization points and the corresponding values of `func`. Each column contains the
visualization points of one element of the layer.

Parallel simulations with MPI are not supported.

# Keywords
- `layer`: the horizontal layer to extract, either `:surface` for the bottom of the
  spherical shell or a positive integer `n` for the bottom face of the `n`-th layer of
  elements counted from the bottom, i.e., `layer = 1` is equivalent to `layer = :surface`.
- `nvisnodes`: number of visualization nodes per element and horizontal direction;
  defaults to twice the polynomial degree of the solver as for 2D plots in Trixi.jl.
"""
function spherical_layer_data(func, u_ode, semi; layer = :surface, nvisnodes = nothing)
    if layer === :surface
        layer = 1
    elseif !(layer isa Integer && layer >= 1)
        throw(ArgumentError("`layer` must be `:surface` or a positive integer, " *
                            "got `layer = $(repr(layer))`"))
    end

    if Trixi.mpi_isparallel()
        throw(ArgumentError("`spherical_layer_data` does not support parallel " *
                            "simulations with MPI"))
    end

    # Setup
    u = Trixi.wrap_array(u_ode, semi)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    (; node_coordinates) = cache.elements

    if nvisnodes === nothing
        nvisnodes = 2 * Trixi.polydeg(solver)
    end

    layer_elements = _layer_elements(layer, solver, cache)

    # Coordinates of the nodes at the bottom face of the elements of the layer
    layer_coordinates = node_coordinates[:, :, :, 1, layer_elements]

    # Evaluate `func` at the nodes at the bottom face of the elements of the layer.
    # `element` is the global element index in `u`, `idx` the position in `layer_elements`.
    layer_values = zeros(eltype(u), size(layer_coordinates)[2:4])
    for (idx, element) in enumerate(layer_elements)
        for j in eachnode(solver), i in eachnode(solver)
            u_node = Trixi.get_node_vars(u, equations, solver, i, j, 1, element)
            layer_values[i, j, idx] = func(u_node, equations)
        end
    end

    n_nodes_2d = Trixi.nnodes(solver)^2
    n_layer_elements = length(layer_elements)
    plotting_interp_matrix = Trixi.plotting_interpolation_matrix(solver; nvisnodes)

    values_2d = reshape(layer_values, n_nodes_2d, n_layer_elements)
    coordinates_2d = reshape(layer_coordinates, 3, n_nodes_2d, n_layer_elements)

    values_vis = plotting_interp_matrix * values_2d
    x_vis = plotting_interp_matrix * coordinates_2d[1, :, :]
    y_vis = plotting_interp_matrix * coordinates_2d[2, :, :]
    z_vis = plotting_interp_matrix * coordinates_2d[3, :, :]

    # Convert the Cartesian coordinates of the visualization points to longitude and
    # latitude in degrees
    lon = similar(x_vis)
    lat = similar(x_vis)
    for n in eachindex(x_vis)
        lon_n, lat_n, _ = cartesian_to_sphere(SVector(x_vis[n], y_vis[n], z_vis[n]))
        lon[n] = rad2deg(lon_n)
        lat[n] = rad2deg(lat_n)
    end

    return (; lon, lat, values = values_vis)
end

function spherical_layer_data(func, sol::Trixi.TrixiODESolution; kwargs...)
    return spherical_layer_data(func, sol.u[end], sol.prob.p; kwargs...)
end

# Return the elements of the `layer`-th layer of elements counted from the bottom of the
# spherical shell. The `n`-th of these elements lies above the `n`-th element at the bottom.
function _layer_elements(layer, solver, cache)
    (; boundaries, interfaces) = cache

    # Find the elements at the bottom of the spherical shell via the `:inside` boundary
    inside_boundaries = findall(==(:inside), boundaries.name)
    if isempty(inside_boundaries)
        throw(ArgumentError("expected a mesh of a spherical shell with the boundary " *
                            "`:inside` at its bottom"))
    end
    layer_elements = boundaries.neighbor_ids[inside_boundaries]

    # The `:inside` boundary must be the bottom face of each element, i.e., the face
    # with the first index in the third local coordinate direction
    expected_node_indices = (:i_forward, :j_forward, :begin)
    if !all(boundaries.node_indices[boundary] == expected_node_indices
            for boundary in inside_boundaries)
        throw(ArgumentError("expected the `:inside` boundary at the bottom face " *
                            "of each element"))
    end

    if layer == 1
        return layer_elements
    end

    # Find the element above each element, i.e., its neighbor across its top face.
    # At such an interface, the index in the third local coordinate direction is `:end`
    # in the lower element and `:begin` in the upper element. At all other interfaces of
    # a spherical shell, these indices are of the form `:i_forward`, `:j_backward`, etc.
    upper_neighbors = zeros(Int, Trixi.nelements(solver, cache))
    for interface in Trixi.eachinterface(solver, cache)
        primary_element = interfaces.neighbor_ids[1, interface]
        secondary_element = interfaces.neighbor_ids[2, interface]
        primary_indices = interfaces.node_indices[1, interface]
        secondary_indices = interfaces.node_indices[2, interface]
        if primary_indices[3] === :end && secondary_indices[3] === :begin
            upper_neighbors[primary_element] = secondary_element
        elseif primary_indices[3] === :begin && secondary_indices[3] === :end
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

# Triangulate the visualization points returned by `spherical_layer_data` in the
# longitude-latitude plane for contour plots, e.g., with `Makie.tricontourf`.
#
# As for contour plots of DG solutions in Trixi.jl, the visualization points of each element
# are triangulated separately using the triangulation of the visualization points on the
# reference element. In contrast to a Delaunay triangulation of all points, this does not
# connect points of different elements across discontinuities of the solution and avoids
# duplicate points, since the points at element interfaces belong to several elements.
# As in the plots of the baroclinic instability by Artiano and Ranocha (see
# https://github.com/MarcoArtiano/2026_nonconservative_sbp), the longitude is shifted to the
# range from `lon_min` to `lon_min + 360°` such that the periodic boundary of the longitude
# is at the boundary of the plotted range starting at `lon_min`. Elements crossing this
# periodic boundary are made contiguous in longitude. If they extend beyond `lon_min + 360°`,
# they are added a second time shifted by -360° such that both ends of the range are
# covered.
# A pole is a line in longitude-latitude coordinates. Thus, a triangle with a vertex at a
# pole, where the longitude is not defined, is replaced by the quadrilateral between this
# line and the opposite edge of the triangle.
# In an element containing a pole in its interior, the longitude jumps by 360° somewhere.
# The longitude of the vertices of the triangles crossing this jump is shifted accordingly.
# Only a triangle enclosing a pole cannot be represented and is removed.
#
# Return a named tuple `(; lon, lat, values, triangles)` of vectors of the coordinates and
# values of the points and a `3 × n_triangles` matrix of the indices of the vertices of
# each triangle.
function _spherical_layer_triangulation(data; lon_min = 0)
    # Triangulation of the equally spaced visualization points on the reference element
    # (see `Trixi.plotting_interpolation_matrix`), with the first index running fastest
    nvisnodes = isqrt(size(data.lon, 1))
    reference_nodes = LinRange(-1, 1, nvisnodes)
    r = vec([reference_nodes[i] for i in 1:nvisnodes, j in 1:nvisnodes])
    s = vec([reference_nodes[j] for i in 1:nvisnodes, j in 1:nvisnodes])
    reference_triangles = Trixi.reference_plotting_triangulation((r, s))

    lon = eltype(data.lon)[]
    lat = eltype(data.lat)[]
    values = eltype(data.values)[]
    triangles = Int[]

    for element in axes(data.lon, 2)
        lon_element = @. lon_min + mod(data.lon[:, element] - lon_min, 360)
        if maximum(lon_element) - minimum(lon_element) > 180
            @. lon_element = ifelse(lon_element < lon_min + 180, lon_element + 360,
                                    lon_element)
        end
        shifts = maximum(lon_element) > lon_min + 360 ? (0, -360) : (0,)

        for shift in shifts
            offset = length(lon)
            append!(lon, lon_element .+ shift)
            append!(lat, view(data.lat, :, element))
            append!(values, view(data.values, :, element))

            for triangle in axes(reference_triangles, 1)
                vertices = ntuple(k -> offset + Int(reference_triangles[triangle, k]),
                                  3)
                _push_triangle!(triangles, lon, lat, values, vertices)
            end
        end
    end

    return (; lon, lat, values, triangles = reshape(triangles, 3, :))
end

function _push_triangle!(triangles, lon, lat, values, vertices)
    at_pole = map(vertex -> isapprox(abs(lat[vertex]), 90), vertices)
    n_vertices_at_pole = count(at_pole)

    if n_vertices_at_pole == 0
        _push_triangle_contiguous!(triangles, lon, lat, values, vertices)
    elseif n_vertices_at_pole == 1
        # Replace the triangle by the quadrilateral between the pole and the opposite edge
        # from `vertex_1` to `vertex_2`, using the longitudes of `vertex_1` and `vertex_2`
        # for two copies of the vertex at the pole
        k = findfirst(at_pole)
        pole = vertices[k]
        vertex_1 = vertices[mod1(k + 1, 3)]
        vertex_2 = vertices[mod1(k + 2, 3)]
        pole_1 = _push_point!(lon, lat, values, lon[vertex_1], lat[pole], values[pole])
        pole_2 = _push_point!(lon, lat, values, lon[vertex_2], lat[pole], values[pole])
        _push_triangle_contiguous!(triangles, lon, lat, values,
                                   (pole_1, vertex_1, vertex_2))
        _push_triangle_contiguous!(triangles, lon, lat, values,
                                   (pole_1, vertex_2, pole_2))
    end
    # Triangles with more than one vertex at a pole are degenerate and removed

    return triangles
end

function _push_triangle_contiguous!(triangles, lon, lat, values, vertices)
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
        return _push_point!(lon, lat, values, lon_shifted, lat[vertex], values[vertex])
    end

    # Otherwise, the triangle encloses a pole and cannot be represented
    lon_min, lon_max = extrema(vertex -> lon[vertex], shifted_vertices)
    if lon_max - lon_min <= 180
        append!(triangles, shifted_vertices)
    end

    return triangles
end

function _push_point!(lon, lat, values, lon_point, lat_point, value)
    push!(lon, lon_point)
    push!(lat, lat_point)
    push!(values, value)
    return length(lon)
end

"""
    plot_spherical_layer(func, u_ode, semi; kwargs...)
    plot_spherical_layer(func, sol; kwargs...)

Create a filled contour plot in longitude-latitude coordinates of the scalar quantity
`func(u, equations)` on a horizontal layer of a simulation on a spherical shell, using
[`spherical_layer_data`](@ref). Return a `Makie.Figure` with the plot and a colorbar.

This requires Makie.jl (e.g., via CairoMakie.jl or GLMakie.jl) to be loaded.

The visualization points of each element are triangulated separately in the
longitude-latitude plane. If a pole lies inside a cell of visualization points (e.g., for
an odd number of trees per cube face dimension and an even `nvisnodes`), this cell cannot
be represented and a small region around the pole remains empty.

# Keywords
- `layer`, `nvisnodes`: passed to [`spherical_layer_data`](@ref).
- `lon_limits`, `lat_limits`: plotted range of longitude and latitude in degrees. The
  range of longitude can include 0°, e.g., `lon_limits = (-60, 180)` for the range from
  60°W to 180°.
- `levels`: number of equally spaced contour bands or vector of contour levels.
- `colormap`: colormap of the filled contours.
- `title`: title of the plot.

# Examples
The surface pressure in hPa of the baroclinic instability (see
[`initial_condition_baroclinic_instability`](@ref)) on the northern hemisphere from 0° to
120°W can be plotted after a simulation with the elixir
`elixir_potential_temperature_baroclinic_instability.jl` in
`examples/euler/dry_air/global_circulation` via
```julia
using CairoMakie
fig = plot_spherical_layer((u, equations) -> pressure(u, equations) / 100, sol;
                           lon_limits = (0, 240), lat_limits = (0, 90),
                           levels = 920:10:1030, colormap = :plasma,
                           title = "Surface pressure [hPa]")
```
"""
function plot_spherical_layer end

"""
    plot_spherical_layer!(ax, func, u_ode, semi; kwargs...)
    plot_spherical_layer!(ax, func, sol; kwargs...)

Add a filled contour plot with contour lines of `func(u, equations)` on a horizontal layer
of a simulation on a spherical shell to the existing Makie.jl axis `ax` and return the
filled contour plot, e.g., to create a `Makie.Colorbar`. The keywords are the same as for
[`plot_spherical_layer`](@ref) except for `title`. This requires Makie.jl to be loaded.
"""
function plot_spherical_layer! end
end # @muladd
