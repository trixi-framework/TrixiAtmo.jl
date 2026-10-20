# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

"""
    PlotDataSphere(u_ode, semi; level = 1, solution_variables = cons2prim,
                   nvisnodes = 2 * nnodes(semi.solver))
    PlotDataSphere(sol; kwargs...)

Data of a solution on a [`P4estMeshCubedSphereTopography`](@ref) on the horizontal level
`level` for plotting in longitude and latitude, e.g., with
`contourf(pd["p1"])` of the Makie.jl extension of TrixiAtmo.jl. The levels of nodes are
counted from the surface (`level = 1`) to the top, where the nodes on the interfaces between
two layers of elements belong to one level. With topography, the levels follow the terrain.
On each element, the solution is interpolated to `nvisnodes` equidistant nodes in each
horizontal direction and triangulated. The variables are computed by
`solution_variables` and selected by their names, e.g., `pd["p1"]`.
"""
struct PlotDataSphere{RealT <: Real}
    longitude::Vector{RealT} # in degrees
    latitude::Vector{RealT} # in degrees
    data::Matrix{RealT} # `[point, variable]`
    triangles::Matrix{Int} # `[vertex, triangle]`
    variable_names::Vector{String}
    level::Int
end

# One variable of a `PlotDataSphere`, see `getindex(::PlotDataSphere, variable_name)`
struct PlotDataSphereSeries{PlotData <: PlotDataSphere}
    plot_data::PlotData
    variable_id::Int
end

function Base.getindex(pd::PlotDataSphere, variable_name)
    variable_id = findfirst(isequal(variable_name), pd.variable_names)
    if isnothing(variable_id)
        throw(KeyError(variable_name))
    end
    return PlotDataSphereSeries(pd, variable_id)
end

Base.keys(pd::PlotDataSphere) = tuple(pd.variable_names...)

function PlotDataSphere(sol::Trixi.TrixiODESolution; kwargs...)
    return PlotDataSphere(sol.u[end], sol.prob.p; kwargs...)
end

function PlotDataSphere(u_ode, semi; level = 1, solution_variables = cons2prim,
                        nvisnodes = 2 * nnodes(semi.solver))
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    nlayers = number_of_layers_cubed_sphere(mesh)
    if nelements(solver, cache) != size(mesh.tree_node_coordinates, ndims(mesh) + 2)
        throw(ArgumentError("plotting is not supported for refined meshes"))
    end

    # The top node of a layer and the bottom node of the layer above belong to one level
    nnodes_ = nnodes(solver)
    nlevels = nlayers * (nnodes_ - 1) + 1
    if !(1 <= level <= nlevels)
        throw(ArgumentError("`level` must be between 1 and $nlevels"))
    end
    layer = min((level - 1) ÷ (nnodes_ - 1) + 1, nlayers)
    k = level - (layer - 1) * (nnodes_ - 1)

    # Interpolation to the plotting nodes and triangulation of the reference element as in
    # `Trixi.PlotData2DTriangulated`
    interpolation_matrix = Trixi.plotting_interpolation_matrix(solver; nvisnodes)
    r, s = Trixi.reference_node_coordinates_2d(solver)
    reference_triangles = Trixi.reference_plotting_triangulation((interpolation_matrix *
                                                                  r,
                                                                  interpolation_matrix *
                                                                  s))

    # The elements of the cubed sphere are ordered column by column, i.e., the elements of a
    # layer are every `nlayers`-th element
    elements = layer:nlayers:nelements(solver, cache)
    variable_names = collect(String, varnames(solution_variables, equations))
    nvariables_ = length(variable_names)
    nvisnodes_2d = nvisnodes^2
    npoints = nvisnodes_2d * length(elements)
    RealT = eltype(u_ode)
    x = Matrix{RealT}(undef, npoints, 3)
    data = Matrix{RealT}(undef, npoints, nvariables_)
    triangles = Matrix{Int}(undef, 3, size(reference_triangles, 1) * length(elements))
    x_nodes = Matrix{RealT}(undef, nnodes_^2, 3)
    data_nodes = Matrix{RealT}(undef, nnodes_^2, nvariables_)

    u = wrap_array(u_ode, semi)
    (; node_coordinates) = cache.elements
    for (index, element) in enumerate(elements)
        for j in eachnode(solver), i in eachnode(solver)
            node = i + (j - 1) * nnodes_
            for d in 1:3
                x_nodes[node, d] = node_coordinates[d, i, j, k, element]
            end
            u_node = get_node_vars(u, equations, solver, i, j, k, element)
            data_nodes[node, :] .= solution_variables(u_node, equations)
        end
        offset = (index - 1) * nvisnodes_2d
        points = (offset + 1):(offset + nvisnodes_2d)
        x[points, :] .= interpolation_matrix * x_nodes
        data[points, :] .= interpolation_matrix * data_nodes
        for t in axes(reference_triangles, 1)
            triangle = (index - 1) * size(reference_triangles, 1) + t
            for vertex in 1:3
                triangles[vertex, triangle] = offset + reference_triangles[t, vertex]
            end
        end
    end

    longitude = Vector{RealT}(undef, npoints)
    latitude = Vector{RealT}(undef, npoints)
    for point in 1:npoints
        lambda, phi, _ = cartesian_to_sphere(view(x, point, :))
        longitude[point] = lambda > pi ? rad2deg(lambda) - 360 : rad2deg(lambda)
        latitude[point] = rad2deg(phi)
    end

    longitude, latitude, data, triangles = split_triangles_at_date_line(longitude,
                                                                        latitude, data,
                                                                        triangles)

    return PlotDataSphere(longitude, latitude, data, triangles, variable_names, level)
end

# Number of layers of elements of a `P4estMeshCubedSphereTopography`
function number_of_layers_cubed_sphere(mesh::P4estMesh{3})
    key = UInt64(Trixi.P4est.PointerWrappers.pointer(mesh.p4est))
    if !haskey(MeshContextRegistry.REGISTRY, key)
        throw(ArgumentError("plotting on the sphere requires a " *
                            "`P4estMeshCubedSphereTopography`"))
    end
    return Int(MeshContextRegistry.REGISTRY[key].num_layers_per_column)
end

function number_of_layers_cubed_sphere(mesh)
    throw(ArgumentError("plotting on the sphere requires a " *
                        "`P4estMeshCubedSphereTopography`"))
end

# Triangles crossing the date line (longitude ±180°) would span the whole plot, and the
# longitude of points at the poles is arbitrary. Thus, the longitudes of the vertices of each
# triangle are taken relative to its first vertex, and the triangles extending beyond ±180°
# are added a second time shifted by ∓360°, such that both parts are shown within the
# limits of the plot.
function split_triangles_at_date_line(longitude, latitude, data, triangles)
    RealT = eltype(longitude)
    npoints = length(longitude)
    new_points = Int[]
    new_longitudes = RealT[]
    new_triangles = NTuple{3, Int}[]
    for t in axes(triangles, 2)
        ids = (triangles[1, t], triangles[2, t], triangles[3, t])
        longitudes = map(ids) do id
            return longitude[ids[1]] +
                   mod(longitude[id] - longitude[ids[1]] + 180, 360) - 180
        end
        at_pole = map(id -> abs(latitude[id]) > 90 - sqrt(eps(RealT)), ids)
        if count(at_pole) == 1
            # The triangle is split into two triangles, where the point at the pole gets
            # the longitudes of the other two vertices, such that the plot reaches the pole
            pole = findfirst(at_pole)
            v1, v2 = mod1(pole + 1, 3), mod1(pole + 2, 3)
            add_triangle!(new_triangles, new_points, new_longitudes, longitude, npoints,
                          (ids[pole], ids[v1], ids[v2]),
                          (longitudes[v1], longitudes[v1], longitudes[v2]))
            add_triangle!(new_triangles, new_points, new_longitudes, longitude, npoints,
                          (ids[pole], ids[v2], ids[pole]),
                          (longitudes[v1], longitudes[v2], longitudes[v2]))
        else
            add_triangle!(new_triangles, new_points, new_longitudes, longitude, npoints,
                          ids, longitudes)
        end
    end

    longitude = vcat(longitude, new_longitudes)
    latitude = vcat(latitude, latitude[new_points])
    data = vcat(data, data[new_points, :])
    triangles = [triangle[vertex] for vertex in 1:3, triangle in new_triangles]
    return longitude, latitude, data, triangles
end

# Add the triangle with the points `ids` and the (modified) longitudes `longitudes`, a
# second time shifted by ∓360° if it extends beyond ±180°. Points with modified longitudes
# are added as new points.
function add_triangle!(new_triangles, new_points, new_longitudes, longitude, npoints,
                       ids,
                       longitudes)
    shifts = (0, maximum(longitudes) > 180 ? -360 : 0,
              minimum(longitudes) < -180 ? 360 : 0)
    for (n, shift) in enumerate(shifts)
        (n > 1 && iszero(shift)) && continue
        new_ids = ntuple(3) do v
            if longitudes[v] + shift == longitude[ids[v]]
                return ids[v]
            end
            push!(new_points, ids[v])
            push!(new_longitudes, longitudes[v] + shift)
            return npoints + length(new_points)
        end
        push!(new_triangles, new_ids)
    end

    return nothing
end
end # @muladd
