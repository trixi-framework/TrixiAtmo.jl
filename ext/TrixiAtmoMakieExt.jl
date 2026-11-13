# Package extension for adding Makie-based features to TrixiAtmo.jl
module TrixiAtmoMakieExt

# Required for visualization code
using Makie: Makie

using TrixiAtmo: TrixiAtmo, @muladd, spherical_layer_data

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
    ax = Makie.Axis(fig[1, 1]; title, xlabel = "longitude [°]", ylabel = "latitude [°]")
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

    # Makie triangulates the scattered visualization points
    lon = vec(data.lon)
    lat = vec(data.lat)
    values = vec(data.values)

    # Equally spaced bands between the extrema.
    if levels isa Integer
        levels = range(extrema(values)...; length = levels + 1)
    end

    plt = Makie.tricontourf!(ax, lon, lat, values; levels, colormap)
    Makie.tricontour!(ax, lon, lat, values; levels, color = :black, linewidth = 0.5)
    Makie.limits!(ax, lon_limits..., lat_limits...)

    return plt
end
end # @muladd

end # module TrixiAtmoMakieExt
