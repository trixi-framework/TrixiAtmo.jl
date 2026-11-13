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
for visualization.

The function `func(u, equations)` takes the conservative variables `u` at a node and the
`equations` and returns the scalar quantity to be visualized, e.g., `pressure`.

Return a named tuple `(; lon, lat, values)` of matrices of size
`(nvisnodes^2, n_surface_elements)` with the longitude and latitude (in degrees) of the
visualization points and the corresponding values of `func`. Each column contains the
visualization points of one element.

# Keywords
- `layer`: the horizontal layer to extract; only `:surface` is supported.
- `nvisnodes`: number of visualization nodes per element and horizontal direction;
  defaults to the number of nodes of the solver.
"""
function spherical_layer_data(func, u_ode, semi; layer = :surface, nvisnodes = nothing)
    if layer !== :surface
        throw(ArgumentError("`layer = $(repr(layer))` is not supported; " *
                            "currently only `layer = :surface` is available"))
    end

    # Setup
    u = Trixi.wrap_array(u_ode, semi)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    (; node_coordinates) = cache.elements

    if nvisnodes === nothing
        nvisnodes = Trixi.nnodes(solver)
    end

    # Find the elements at the bottom of the spherical shell via the `:inside` boundary
    (; boundaries) = cache
    inside_boundaries = findall(==(:inside), boundaries.name)
    surface_elements = boundaries.neighbor_ids[inside_boundaries]

    # The `:inside` boundary must be the bottom face of each element, i.e., consist of the nodes
    expected_node_indices = (:i_forward, :j_forward, :begin)
    if !all(boundaries.node_indices[boundary] == expected_node_indices
            for boundary in inside_boundaries)
        throw(ArgumentError("expected the `:inside` boundary at the bottom face" *
                            "of each element"))
    end

    # Coordinates of the nodes at the bottom of the surface elements
    surface_coordinates = node_coordinates[:, :, :, 1, surface_elements]

    # Evaluate `func` at the nodes at the bottom of the surface elements
    # `element` is the global element index in `u`, `idx` the position in `surface_elements`.
    surface_values = zeros(eltype(u), size(surface_coordinates)[2:4])
    for (idx, element) in enumerate(surface_elements)
        for j in eachnode(solver), i in eachnode(solver)
            u_node = Trixi.get_node_vars(u, equations, solver, i, j, 1, element)
            surface_values[i, j, idx] = func(u_node, equations)
        end
    end

    n_nodes_2d = Trixi.nnodes(solver)^2
    n_surface_elements = length(surface_elements)
    plotting_interp_matrix = Trixi.plotting_interpolation_matrix(solver; nvisnodes)

    values_2d = reshape(surface_values, n_nodes_2d, n_surface_elements)
    coordinates_2d = reshape(surface_coordinates, 3, n_nodes_2d, n_surface_elements)

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

"""
    plot_spherical_layer(func, u_ode, semi; kwargs...)
    plot_spherical_layer(func, sol; kwargs...)

Create a filled contour plot in longitude-latitude coordinates of the scalar quantity
`func(u, equations)` on a horizontal layer of a simulation on a spherical shell, using
[`spherical_layer_data`](@ref). Return a `Makie.Figure` with the plot and a colorbar.

This requires Makie.jl (e.g., via CairoMakie.jl or GLMakie.jl) to be loaded.

# Keywords
- `layer`, `nvisnodes`: passed to [`spherical_layer_data`](@ref).
- `lon_limits`, `lat_limits`: plotted range of longitude and latitude in degrees.
- `levels`: number of equally spaced contour bands or vector of contour levels.
- `colormap`: colormap of the filled contours.
- `title`: title of the plot.
"""
function plot_spherical_layer end

"""
    plot_spherical_layer!(ax, func, u_ode, semi; kwargs...)
    plot_spherical_layer!(ax, func, sol; kwargs...)

Add a filled contour plot with contour lines of `func(u, equations)` on a horizontal layer
of a simulation on a spherical shell to the existing Makie.jl axis `ax` and return the
filled contour plot, e.g., to create a `Makie.Colorbar`. This requires Makie.jl to be
loaded.
"""
function plot_spherical_layer! end
end # @muladd
