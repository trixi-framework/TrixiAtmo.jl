# Package extension for adding Makie-based features to TrixiAtmo.jl
module TrixiAtmoMakieExt

# Required for visualization code
using Makie: Makie

using TrixiAtmo: TrixiAtmo, @muladd, PlotDataShell, LonLat
using TrixiAtmo.Trixi: PlotDataSeries

# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# A `PlotDataShell` in longitude-latitude coordinates and a single variable of it
const PlotDataShellLonLat = PlotDataShell{<:Any, <:Any, <:Any, <:LonLat}
const PlotDataSeriesLonLat = PlotDataSeries{<:PlotDataShellLonLat}

# As for the plot data of Trixi.jl, filled contour plots of a single variable are created
# with `Makie.contourf` and `Makie.contourf!`. The keywords are documented in the docstring
# of `PlotDataShell`, since methods defined in this extension do not appear in the API
# reference of TrixiAtmo.jl.
function Makie.contourf(pds::PlotDataSeriesLonLat;
                        lon_limits = default_lon_limits(pds.plot_data.projection),
                        lat_limits = (-90, 90),
                        title = pds.plot_data.variable_names[pds.variable_id],
                        kwargs...)
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1]; title, xlabel = "longitude", ylabel = "latitude",
                    xticks = degree_ticks(lon_limits),
                    xtickformat = lons -> format_longitude.(lons),
                    yticks = degree_ticks(lat_limits),
                    ytickformat = lats -> format_latitude.(lats))
    plt = Makie.contourf!(ax, pds; lon_limits, lat_limits, kwargs...)
    Makie.Colorbar(fig[1, 2], plt)

    # Make the width of the axis proportional to the plotted longitude range and its
    # height proportional to the plotted latitude range.
    lon_lat_ratio = (lon_limits[2] - lon_limits[1]) / (lat_limits[2] - lat_limits[1])
    Makie.colsize!(fig.layout, 1, Makie.Aspect(1, lon_lat_ratio))
    Makie.resize_to_layout!(fig)

    # As for single plots in Trixi.jl, return the figure together with the axis and the plot
    return Makie.FigureAxisPlot(fig, ax, plt)
end

function Makie.contourf!(ax, pds::PlotDataSeriesLonLat;
                         lon_limits = default_lon_limits(pds.plot_data.projection),
                         lat_limits = (-90, 90), colormap = :viridis, kwargs...)
    (; plot_data, variable_id) = pds
    (; x, y, triangles, data, projection) = plot_data
    check_lon_limits(projection, lon_limits)

    # The visualization points are already triangulated in the longitude-latitude plane
    # with the periodic boundary of the longitude at the lower limit `lon_min` of the projection.
    # We use this triangulation instead of letting Makie compute a Delaunay triangulation of all
    # points, but pass only the triangles in the plotted range (see `crop_triangulation`).
    x, y, values, triangles = crop_triangulation(x, y, view(data, :, variable_id),
                                                 triangles, lon_limits, lat_limits)

    # All other keywords (e.g., `levels`, `extendlow`, and `extendhigh`) are passed only to
    # the filled contours, which are returned
    plt = Makie.tricontourf!(ax, x, y, values; triangulation = triangles, colormap,
                             kwargs...)

    # Draw the contour lines at the edges of the bands of the filled contours. For an integer
    # `levels`, `Makie.tricontour!` would choose other levels than `Makie.tricontourf!`.
    # As in the colorbar of Makie.jl, the computed levels are transformed back if a
    # `colorscale` (e.g., `log10`) is used.
    levels = Makie.apply_scale(Makie.inverse_transform(plt.colorscale[]),
                               plt.computed_levels[])
    Makie.tricontour!(ax, x, y, values; triangulation = triangles, levels,
                      color = :black, linewidth = 0.5)
    Makie.limits!(ax, lon_limits..., lat_limits...)

    return plt
end

# A `PlotDataShell` with a single variable, e.g., a scalar quantity such as the pressure,
# can be plotted directly
function Makie.contourf(pd::PlotDataShellLonLat; kwargs...)
    return Makie.contourf(single_variable(pd); kwargs...)
end

function Makie.contourf!(ax, pd::PlotDataShellLonLat; kwargs...)
    return Makie.contourf!(ax, single_variable(pd); kwargs...)
end

function single_variable(pd::PlotDataShell)
    if length(pd) != 1
        throw(ArgumentError("the `PlotDataShell` contains $(length(pd)) variables. " *
                            "Select a single variable via `pd[variable_name]`, e.g., " *
                            "`pd[\"$(first(keys(pd)))\"]`."))
    end
    return PlotDataSeries(pd, 1)
end

# The full range of longitudes of the projection
default_lon_limits(projection::LonLat) = (projection.lon_min, projection.lon_min + 360)

# The plotted range of longitudes must lie in the range covered by the triangulation
function check_lon_limits(projection::LonLat, lon_limits)
    lon_min, lon_max = default_lon_limits(projection)
    if lon_limits[1] < lon_min || lon_limits[2] > lon_max
        throw(ArgumentError("`lon_limits = $lon_limits` must lie in the range from " *
                            "$lon_min to $lon_max of the projection " *
                            "`LonLat(lon_min = $lon_min)`. Create the `PlotDataShell` " *
                            "with, e.g., `projection = LonLat(lon_min = " *
                            "$(lon_limits[1]))` to plot this range."))
    end
    return nothing
end

# Return the coordinates `lon` and `lat`, the `values`, and the `triangles` restricted to the
# triangles whose bounding box overlaps the plotted range of longitude and latitude, with the
# vertices of these triangles numbered consecutively.
# Makie.jl computes the contour levels for an integer `levels` (and for `mode = :relative`)
# from all values passed to it, independently of the limits of the axis. Since a
# `PlotDataShell` always contains the whole sphere, we pass only these triangles to Makie.jl
# such that the levels are chosen according to the values in the plotted range. Since the
# values in a triangle are interpolated linearly from its vertices, the levels cover all
# visible values, and there are no gaps in the filled contours at the boundary of the axis.
# If no triangle overlaps the plotted range, all triangles are returned.
function crop_triangulation(lon, lat, values, triangles, lon_limits, lat_limits)
    is_visible = map(eachcol(triangles)) do triangle
        lon_min, lon_max = extrema(vertex -> lon[vertex], triangle)
        lat_min, lat_max = extrema(vertex -> lat[vertex], triangle)
        return (lon_max >= lon_limits[1] && lon_min <= lon_limits[2] &&
                lat_max >= lat_limits[1] && lat_min <= lat_limits[2])
    end
    if !any(is_visible)
        return lon, lat, values, triangles
    end
    visible_triangles = triangles[:, is_visible]

    # Keep only the vertices of the visible triangles and number them consecutively
    is_vertex = falses(length(lon))
    is_vertex[visible_triangles] .= true
    vertices = findall(is_vertex)
    new_indices = zeros(Int, length(lon))
    new_indices[vertices] .= eachindex(vertices)

    return lon[vertices], lat[vertices], values[vertices],
           new_indices[visible_triangles]
end

# Ticks at multiples of 30° within the plotted range. For small
# ranges, Makie.jl chooses the ticks automatically.
function degree_ticks(limits)
    ticks = (30 * ceil(Int, limits[1] / 30)):30:(30 * floor(Int, limits[2] / 30))
    return length(ticks) >= 2 ? ticks : Makie.automatic
end

# Labels of the ticks in degrees east or west, e.g., `30E` or `120W` (for 240°)
function format_longitude(lon)
    lon = mod(lon, 360)
    if iszero(lon) || lon == 180
        return format_degrees(lon)
    elseif lon < 180
        return format_degrees(lon) * "E"
    else
        return format_degrees(360 - lon) * "W"
    end
end

# Labels of the ticks in degrees north or south, e.g., `60N` or `30S`
function format_latitude(lat)
    if iszero(lat)
        return format_degrees(lat)
    elseif lat > 0
        return format_degrees(lat) * "N"
    else
        return format_degrees(-lat) * "S"
    end
end

format_degrees(x) = isinteger(x) ? string(Int(x)) : string(x)
end # @muladd

end # module TrixiAtmoMakieExt
