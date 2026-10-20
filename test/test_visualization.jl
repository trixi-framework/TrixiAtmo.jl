@testsnippet Visualization begin
    EXAMPLES_DIR = joinpath(examples_dir(), "euler", "dry_air", "global_circulation")
end

@testitem "Visualization: PlotDataSphere" setup=[Setup, Visualization] tags=[:visualization] begin
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR,
                           "elixir_potential_temperature_baroclinic_instability.jl"),
                  trees_per_cube_face = (2, 2), polydeg = 3, sol = nothing)

    # The geopotential `phi = -a^2 g / r` of the initial condition is constant on the
    # levels of a mesh without topography
    pd = PlotDataSphere(ode.u0, semi)
    @test pd.level == 1
    @test collect(keys(pd)) == collect(varnames(cons2prim, equations))
    phi = view(pd.data, :, pd["phi"].variable_id)
    @test all(isapprox.(phi, -EARTH_RADIUS * EARTH_GRAVITATIONAL_ACCELERATION))
    @test all(lat -> -90 <= lat <= 90, pd.latitude)

    # The top level is at the height of 30 km
    nlevels = 2 * (nnodes(solver) - 1) + 1
    pd = PlotDataSphere(ode.u0, semi; level = nlevels)
    phi = view(pd.data, :, pd["phi"].variable_id)
    @test all(isapprox.(phi,
                        -EARTH_RADIUS^2 * EARTH_GRAVITATIONAL_ACCELERATION /
                        (EARTH_RADIUS + 30_000)))
    @test_throws ArgumentError PlotDataSphere(ode.u0, semi; level = nlevels + 1)
    @test_throws KeyError pd["temperature"]

    # Away from the poles, no triangle spans the whole plot, e.g., at the date line
    @test all(eachcol(pd.triangles)) do triangle
        longitudes = pd.longitude[triangle]
        return maximum(abs, pd.latitude[triangle]) > 80 ||
               maximum(longitudes) - minimum(longitudes) < 90
    end
end

@testitem "Visualization: contourf with CairoMakie" setup=[Setup, Visualization] tags=[:visualization] begin
    # Using CairoMakie inside the testitem since otherwise the bindings are hidden by the anonymous modules
    using CairoMakie
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR,
                           "elixir_potential_temperature_baroclinic_instability.jl"),
                  trees_per_cube_face = (2, 2), polydeg = 3, sol = nothing)
    pd = PlotDataSphere(ode.u0, semi)

    fig, ax, plot = contourf(pd["p1"]; levels = 10)
    @test fig isa Figure
    filename = joinpath(mktempdir(), "surface_pressure.png")
    save(filename, fig)
    @test isfile(filename)

    fig = Figure()
    ax = Axis(fig[1, 1])
    @test contourf!(ax, pd["v1"]; levels = 10) isa Makie.Tricontourf
end
