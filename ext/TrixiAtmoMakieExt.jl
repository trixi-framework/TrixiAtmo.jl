# Package extension for plotting with Makie.jl
module TrixiAtmoMakieExt

using Makie: Makie, Figure, Axis, Colorbar, DataAspect, tricontourf!
using TrixiAtmo: PlotDataSphereSeries

"""
    contourf!(ax, pds::PlotDataSphereSeries; kwargs...)

Filled contour plot of the variable `pds`, e.g., `pd["p1"]` of a [`PlotDataSphere`](@ref) `pd`,
in longitude and latitude (in degrees) into the axis `ax`. The keyword arguments are passed
to `Makie.tricontourf!`, e.g., `levels` and `colormap`.
"""
function Makie.contourf!(ax, pds::PlotDataSphereSeries; kwargs...)
    pd = pds.plot_data
    return tricontourf!(ax, pd.longitude, pd.latitude, view(pd.data, :, pds.variable_id);
                        triangulation = pd.triangles, kwargs...)
end

"""
    contourf(pds::PlotDataSphereSeries; kwargs...)

Figure with a filled contour plot of the variable `pds`, e.g., `pd["p1"]` of a
[`PlotDataSphere`](@ref) `pd`, in longitude and latitude (in degrees) and a colorbar. The
keyword arguments are passed to `Makie.tricontourf!`, e.g., `levels` and `colormap`.
"""
function Makie.contourf(pds::PlotDataSphereSeries; kwargs...)
    pd = pds.plot_data
    fig = Figure()
    ax = Axis(fig[1, 1]; title = pd.variable_names[pds.variable_id],
              xlabel = "longitude", ylabel = "latitude", aspect = DataAspect(),
              limits = (-180, 180, -90, 90), xticks = -180:60:180, yticks = -90:30:90)
    plot = Makie.contourf!(ax, pds; kwargs...)
    Colorbar(fig[1, 2], plot)
    return Makie.FigureAxisPlot(fig, ax, plot)
end

end # module TrixiAtmoMakieExt
