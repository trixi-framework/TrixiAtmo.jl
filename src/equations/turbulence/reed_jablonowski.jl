@muladd begin
#! format: noindent

@doc raw"""
    ReedJablonowskiBoundaryLayer(; pressure_top = 85000, pressure_decay = 10000,
                                 exchange_coefficient = 0.0011,
                                 drag_coefficients = (7.0e-4, 6.5e-5),
                                 drag_coefficient_max = 0.002,
                                 wind_speed_threshold = 20)

Vertical model of the planetary boundary layer of Reed and Jablonowski (2012) as in the
simple physics of DCMIP2016, used as vertical model of a
[`HorizontalVerticalTurbulence`](@ref). The eddy viscosity and diffusivity are
```math
\nu = C_d \lVert \vec{v}_a \rVert z_a f(p), \qquad
\kappa = C_E \lVert \vec{v}_a \rVert z_a f(p), \qquad
f(p) = \begin{cases}
  1, & p \geq p_{top}, \\
  \exp\left(-\left(\frac{p_{top} - p}{p_{decay}}\right)^2\right), & p < p_{top},
\end{cases}
```
with the pressure ``p`` of the node, the top of the boundary layer ``p_{top}``
`= pressure_top` and the decay ``p_{decay}`` `= pressure_decay` in Pa, the exchange
coefficient ``C_E`` `= exchange_coefficient` of heat and moisture, and the drag coefficient
``C_d`` of [`drag_coefficient`](@ref). Here, ``\lVert \vec{v}_a \rVert`` is the speed of
the horizontal velocity and ``z_a`` the height above the surface of the lowest level of the
column of the node, see [`TurbulentDiffusion3D`](@ref). The default values are the
ones of DCMIP2016. The corresponding fluxes through the surface are given by
[`BoundaryConditionSurfaceFluxes`](@ref).

## References

- Kevin A. Reed, Christiane Jablonowski (2012)
  Idealized tropical cyclone simulations of intermediate complexity: A test case for AGCMs
  [DOI: 10.1029/2011MS000099](https://doi.org/10.1029/2011MS000099)
- Paul A. Ullrich et al. (2016)
  Dynamical Core Model Intercomparison Project (DCMIP2016) test case document
  [DOI: 10.5281/zenodo.1298671](https://doi.org/10.5281/zenodo.1298671)
"""
struct ReedJablonowskiBoundaryLayer{RealT <: Real}
    pressure_top::RealT # top of the boundary layer in Pa
    pressure_decay::RealT # decay of the diffusivities above the boundary layer in Pa
    exchange_coefficient::RealT # bulk exchange coefficient of heat and moisture
    drag_coefficients::NTuple{2, RealT} # (C_d0, C_d1 in s/m) of `drag_coefficient`
    drag_coefficient_max::RealT # drag coefficient for strong winds
    wind_speed_threshold::RealT # wind speed in m/s above which `drag_coefficient_max` holds
end

function ReedJablonowskiBoundaryLayer(; pressure_top = 85000, pressure_decay = 10000,
                                      exchange_coefficient = 0.0011,
                                      drag_coefficients = (7.0e-4, 6.5e-5),
                                      drag_coefficient_max = 0.002,
                                      wind_speed_threshold = 20)
    values = promote(pressure_top, pressure_decay, exchange_coefficient,
                     drag_coefficients..., drag_coefficient_max, wind_speed_threshold,
                     0.0)
    return ReedJablonowskiBoundaryLayer(values[1], values[2], values[3],
                                        (values[4], values[5]), values[6], values[7])
end

function Base.similar(turbulence::ReedJablonowskiBoundaryLayer,
                      ::Type{NewRealT}) where {NewRealT}
    return ReedJablonowskiBoundaryLayer(convert(NewRealT, turbulence.pressure_top),
                                        convert(NewRealT, turbulence.pressure_decay),
                                        convert(NewRealT,
                                                turbulence.exchange_coefficient),
                                        convert.(NewRealT,
                                                 turbulence.drag_coefficients),
                                        convert(NewRealT,
                                                turbulence.drag_coefficient_max),
                                        convert(NewRealT,
                                                turbulence.wind_speed_threshold))
end

@inline nturbulencevariables(::ReedJablonowskiBoundaryLayer) = 0

@doc raw"""
    drag_coefficient(wind_speed, turbulence::ReedJablonowskiBoundaryLayer)

Drag coefficient ``C_d = C_{d0} + C_{d1} \lVert \vec{v}_a \rVert`` for wind speeds
``\lVert \vec{v}_a \rVert`` below the threshold and ``C_d = C_{d,max}`` above, see
[`ReedJablonowskiBoundaryLayer`](@ref).
"""
@inline function drag_coefficient(wind_speed, turbulence::ReedJablonowskiBoundaryLayer)
    c_d0, c_d1 = turbulence.drag_coefficients
    return ifelse(wind_speed < turbulence.wind_speed_threshold,
                  c_d0 + c_d1 * wind_speed, turbulence.drag_coefficient_max)
end

@inline function eddy_diffusivities(w, gradients, k, column,
                                    turbulence::ReedJablonowskiBoundaryLayer, equations)
    wind_speed, height = column
    p = pressure_specific(w, equations)
    (; pressure_top, pressure_decay) = turbulence
    decay = ifelse(p >= pressure_top, one(p),
                   exp(-((pressure_top - p) / pressure_decay)^2))
    scale = wind_speed * height * decay
    return (drag_coefficient(wind_speed, turbulence) * scale,
            turbulence.exchange_coefficient * scale)
end
end # @muladd
