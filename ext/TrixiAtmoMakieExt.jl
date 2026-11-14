# Package extension for adding Makie-based features to TrixiAtmo.jl
module TrixiAtmoMakieExt

# Required for visualization code
using Makie: Makie

using TrixiAtmo: TrixiAtmo, @muladd, spherical_layer_data, _spherical_layer_triangulation

# Import functions such that they can be extended with new methods
import TrixiAtmo: plot_spherical_layer, plot_spherical_layer!

# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

function plot_spherical_layer(func, args...;
                              layer = :surface, nvisnodes = nothing,
                              lon_limits = (0, 360), lat_limits = (-90, 90),
                              levels = 10, colormap = :viridis, title = "")
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1]; title, xlabel = "longitude", ylabel = "latitude",
                    xticks = degree_ticks(lon_limits),
                    xtickformat = lons -> format_longitude.(lons),
                    yticks = degree_ticks(lat_limits),
                    ytickformat = lats -> format_latitude.(lats))
    plt = plot_spherical_layer!(ax, func, args...; layer, nvisnodes, lon_limits,
                                lat_limits, levels, colormap)
    Makie.Colorbar(fig[1, 2], plt)

    # Make the width of the axis proportional to the plotted longitude range and its
    # height proportional to the plotted latitude range.
    lon_lat_ratio = (lon_limits[2] - lon_limits[1]) / (lat_limits[2] - lat_limits[1])
    Makie.colsize!(fig.layout, 1, Makie.Aspect(1, lon_lat_ratio))
    Makie.resize_to_layout!(fig)

    return fig
end

function plot_spherical_layer!(ax, func, args...;
                               layer = :surface, nvisnodes = nothing,
                               lon_limits = (0, 360), lat_limits = (-90, 90),
                               levels = 10, colormap = :viridis)
    data = spherical_layer_data(func, args...; layer, nvisnodes)

    # Triangulate the visualization points of each element separately instead of letting
    # Makie compute a Delaunay triangulation of all points. The periodic boundary of the
    # longitude is placed at the lower limit of the plotted range.
    (; lon, lat, values, triangles) = _spherical_layer_triangulation(data;
                                                                     lon_min = lon_limits[1])

    # Equally spaced bands between the extrema
    if levels isa Integer
        levels = range(extrema(values)...; length = levels + 1)
    end

    plt = Makie.tricontourf!(ax, lon, lat, values; triangulation = triangles, levels,
                             colormap)
    Makie.tricontour!(ax, lon, lat, values; triangulation = triangles, levels,
                      color = :black, linewidth = 0.5)
    Makie.limits!(ax, lon_limits..., lat_limits...)

    return plt
end

# Ticks at multiples of 30° within the plotted range as usual for global maps. For small
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
