# Package extension for adding Makie-based features to TrixiAtmo.jl
module TrixiAtmoMakieExt

# Required for visualization code
using Makie: Makie

using TrixiAtmo: TrixiAtmo, @muladd, PlotDataShell, LonLat
using TrixiAtmo.Trixi: PlotDataSeries

# Import functions such that they can be extended with new methods
import TrixiAtmo: plot_spherical_layer, plot_spherical_layer!

# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Single variable of a `PlotDataShell` in longitude-latitude coordinates
const PlotDataSeriesLonLat = PlotDataSeries{<:PlotDataShell{<:Any, <:Any, <:Any,
                                                            <:LonLat}}

function plot_spherical_layer(pds::PlotDataSeriesLonLat;
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
    plt = plot_spherical_layer!(ax, pds; lon_limits, lat_limits, kwargs...)
    Makie.Colorbar(fig[1, 2], plt)

    # Make the width of the axis proportional to the plotted longitude range and its
    # height proportional to the plotted latitude range.
    lon_lat_ratio = (lon_limits[2] - lon_limits[1]) / (lat_limits[2] - lat_limits[1])
    Makie.colsize!(fig.layout, 1, Makie.Aspect(1, lon_lat_ratio))
    Makie.resize_to_layout!(fig)

    # As for single plots in Trixi.jl, return the figure together with the axis and the plot
    return Makie.FigureAxisPlot(fig, ax, plt)
end

function plot_spherical_layer!(ax, pds::PlotDataSeriesLonLat;
                               lon_limits = default_lon_limits(pds.plot_data.projection),
                               lat_limits = (-90, 90),
                               levels = 10, colormap = :viridis, kwargs...)
    (; plot_data, variable_id) = pds
    (; x, y, triangles, data, projection) = plot_data
    check_lon_limits(projection, lon_limits)

    # The visualization points are already triangulated in the longitude-latitude plane
    # with the periodic boundary of the longitude at the lower limit `lon_min` of the projection.
    # We use this triangulation instead of letting Makie compute a Delaunay triangulation of all points.
    values = data[:, variable_id]

    if levels isa Integer
        levels = default_levels(levels, x, y, values, triangles, lon_limits, lat_limits)
    end

    # Additional keywords (e.g., `extendlow` and `extendhigh`) are passed only to the
    # filled contours, which are returned
    plt = Makie.tricontourf!(ax, x, y, values; triangulation = triangles, levels,
                             colormap, kwargs...)
    Makie.tricontour!(ax, x, y, values; triangulation = triangles, levels,
                      color = :black, linewidth = 0.5)
    Makie.limits!(ax, lon_limits..., lat_limits...)

    return plt
end

# Plot a scalar quantity `func(u, equations)` directly from the solution, which creates a
# `PlotDataShell` with the periodic boundary of the longitude at the lower limit of the
# plotted range
function plot_spherical_layer(func, args...; layer = :surface, nvisnodes = nothing,
                              lon_limits = (0, 360), title = "", kwargs...)
    pds = scalar_plot_data_series(func, args...; layer, nvisnodes, lon_limits)
    return plot_spherical_layer(pds; lon_limits, title, kwargs...)
end

function plot_spherical_layer!(ax, func, args...; layer = :surface, nvisnodes = nothing,
                               lon_limits = (0, 360), kwargs...)
    pds = scalar_plot_data_series(func, args...; layer, nvisnodes, lon_limits)
    return plot_spherical_layer!(ax, pds; lon_limits, kwargs...)
end

function scalar_plot_data_series(func, args...; layer, nvisnodes, lon_limits)
    pd = PlotDataShell(args...; solution_variables = func, layer, nvisnodes,
                       projection = LonLat(lon_min = lon_limits[1]))
    if length(pd) != 1
        throw(ArgumentError("`func` must return a scalar quantity, got " *
                            "$(length(pd)) variables"))
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

# Return `n_bands + 1` equally spaced contour levels between the extrema of the values in
# the plotted range of longitude and latitude. We choose the levels independently of the
# values outside of this range. To this end, we use the vertices of all triangles whose
# bounding box overlaps the plotted range. Since the values in a triangle are interpolated
# linearly from its vertices, the levels cover all visible values, and there are no gaps in
# the filled contours at the boundary of the axis.
function default_levels(n_bands, lon, lat, values, triangles, lon_limits, lat_limits)
    value_min = typemax(eltype(values))
    value_max = typemin(eltype(values))
    for triangle in eachcol(triangles)
        lon_min, lon_max = extrema(vertex -> lon[vertex], triangle)
        lat_min, lat_max = extrema(vertex -> lat[vertex], triangle)
        if (lon_max >= lon_limits[1] && lon_min <= lon_limits[2] &&
            lat_max >= lat_limits[1] && lat_min <= lat_limits[2])
            for vertex in triangle
                value_min = min(value_min, values[vertex])
                value_max = max(value_max, values[vertex])
            end
        end
    end

    # Use all values if no triangle overlaps the plotted range
    if value_min > value_max
        value_min, value_max = extrema(values)
    end

    return range(value_min, value_max; length = n_bands + 1)
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
