@muladd begin
#! format: noindent

@doc raw"""
    KesslerMicrophysics(; latent_heat = 2.5e6, rho_surface = 1.2, condensation_time = 1)

Warm-rain microphysics of Kessler (1969) in the form of Klemp and Wilhelmson (1978, KW) as in
DCMIP2016, with the moisture variables water vapor, cloud water, and rain water, i.e., the
densities ``\varrho_v``, ``\varrho_c``, and ``\varrho_r`` with the mixing ratios
``q_x = \varrho_x / \varrho_d``. The rain falls with the terminal velocity of
[`terminal_velocity_rain`](@ref). The right-hand sides are
```math
\begin{aligned}
\partial_t \varrho_v &= \varrho_d (E - C), \\
\partial_t \varrho_c &= \varrho_d (C - A - C_r), \\
\partial_t \varrho_r &= \varrho_d (A + C_r - E),
\end{aligned}
```
with
- the condensation (``C > 0``) or evaporation of cloud water (``C < 0``),
  ```math
  C = \frac{1}{\tau} \max\left(\frac{q_v - q_{vs}}
      {1 + q_{vs} \frac{237 \cdot 17.27 L / c_p}{(T - 36)^2}}, -q_c\right),
  ```
  i.e., the saturation adjustment of KW (3.10) relaxed with the constant time scale
  ``\tau`` `= condensation_time` in s, as the condensation rate of Morrison and Grabowski
  (2008, Eq. A3) for a prescribed phase relaxation time, where ``q_{vs}`` is the saturation
  mixing ratio of [`saturation_mixing_ratio_tetens`](@ref),
- the autoconversion ``A = 10^{-3} \max(q_c - 10^{-3}, 0)`` and the accretion
  ``C_r = 2.2 q_c q_r^{0.875}`` of cloud water to rain (KW 2.13),
- the evaporation of rain (KW 2.14), with the density ``r = 10^{-3} \varrho_d`` in g/cm³ and
  the pressure in mb,
  ```math
  E = \frac{(1.6 + 124.9 (r q_r)^{0.2046}) (r q_r)^{0.525}}
      {2.55 \cdot 10^6 / (p q_{vs}) + 5.4 \cdot 10^5}
      \frac{\max(q_{vs} - q_v, 0)}{r q_{vs}},
  ```
and the isobaric latent heating ``\partial_t T = L / c_p (C - E)`` with the latent heat of
vaporization `latent_heat` ``L`` in J/kg and the specific heat ``c_p`` of dry air at constant
pressure of the equations. The default latent heat is the one of DCMIP2016.

## References

- Edwin Kessler (1969)
  On the distribution and continuity of water substance in atmospheric circulations
  [DOI: 10.1007/978-1-935704-36-2_1](https://doi.org/10.1007/978-1-935704-36-2_1)
- Su-Tzai Soong, Yoshimitsu Ogura (1973)
  A comparison between axisymmetric and slab-symmetric cumulus cloud models
  [DOI: 10.1175/1520-0469(1973)030<0879:ACBAAS>2.0.CO;2](https://doi.org/10.1175/1520-0469(1973)030<0879:ACBAAS>2.0.CO;2)
- Joseph B. Klemp, Robert B. Wilhelmson (1978)
  The simulation of three-dimensional convective storm dynamics
  [DOI: 10.1175/1520-0469(1978)035<1070:TSOTDC>2.0.CO;2](https://doi.org/10.1175/1520-0469(1978)035<1070:TSOTDC>2.0.CO;2)
- Hugh Morrison, Wojciech W. Grabowski (2008)
  Modeling supersaturation and subgrid-scale mixing with two-moment bulk warm microphysics
  [DOI: 10.1175/2007JAS2374.1](https://doi.org/10.1175/2007JAS2374.1)
- Paul A. Ullrich et al. (2016)
  Dynamical Core Model Intercomparison Project (DCMIP2016) test case document
  [DOI: 10.5281/zenodo.1298671](https://doi.org/10.5281/zenodo.1298671)
"""
struct KesslerMicrophysics{RealT <: Real}
    latent_heat::RealT # latent heat of vaporization in J/kg
    rho_surface::RealT # density at the surface in the terminal velocity of rain in kg/m³
    condensation_time::RealT # relaxation time of the condensation in s
    # Coefficient of `rho_rain^1.1364` in `rho_rain * V_r` of `terminal_velocity_rain`
    rain_flux_coefficient::RealT
end

function KesslerMicrophysics(; latent_heat = 2.5e6, rho_surface = 1.2,
                             condensation_time = 1)
    latent_heat, rho_surface, condensation_time = promote(latent_heat, rho_surface,
                                                          condensation_time)
    rain_flux_coefficient = convert(typeof(latent_heat), 36.34 * 1.0e-3^0.1364)
    return KesslerMicrophysics(latent_heat, rho_surface, condensation_time,
                               rain_flux_coefficient)
end

function Base.similar(microphysics::KesslerMicrophysics,
                      ::Type{NewRealT}) where {NewRealT}
    (; latent_heat, rho_surface, condensation_time) = microphysics
    return KesslerMicrophysics(; latent_heat = convert(NewRealT, latent_heat),
                               rho_surface = convert(NewRealT, rho_surface),
                               condensation_time = convert(NewRealT, condensation_time))
end

@inline nmoisturevariables(::KesslerMicrophysics) = 3

function varnames_moisture(::typeof(cons2cons), ::KesslerMicrophysics)
    return ("rho_vapor", "rho_cloud", "rho_rain")
end

function varnames_moisture(::typeof(cons2prim), ::KesslerMicrophysics)
    return ("q_vapor", "q_cloud", "q_rain")
end

@inline function density_water(moisture, ::KesslerMicrophysics)
    return moisture[1] + moisture[2] + moisture[3]
end

@doc raw"""
    terminal_velocity_rain(rho_dry, rho_rain, microphysics::KesslerMicrophysics)

Terminal velocity of rain in m/s of Soong and Ogura (1973) with the density adjustment of
Klemp and Wilhelmson (1978, Eq. 2.15),
```math
V_r = 36.34 \left(10^{-3} \varrho_r\right)^{0.1364} \sqrt{\frac{\varrho_{sfc}}{\varrho_d}},
```
with the rain water density ``\varrho_r = \varrho_d q_r`` in kg/m³, i.e.,
``10^{-3} \varrho_r`` in g/cm³. As the base-state density at the ground of Klemp and
Wilhelmson, the density ``\varrho_{sfc}`` is the constant `rho_surface` of the microphysics.
Its default value of 1.2 kg/m³ is close to the mean density of dry air at the ground in the
moist baroclinic wave of DCMIP2016.
"""
@inline function terminal_velocity_rain(rho_dry, rho_rain,
                                        microphysics::KesslerMicrophysics)
    RealT = typeof(rho_rain)
    rho_rain_positive = max(rho_rain, zero(RealT))
    return convert(RealT, 36.34) *
           (convert(RealT, 1.0e-3) * rho_rain_positive)^convert(RealT, 0.1364) *
           sqrt(microphysics.rho_surface / rho_dry)
end

@inline function sedimentation(moisture, rho_dry, microphysics::KesslerMicrophysics)
    RealT = eltype(moisture)
    rho_rain = max(moisture[3], zero(RealT))
    # `rho_rain * V_r` of `terminal_velocity_rain` as a single power of `rho_rain`, whose
    # derivative is finite (zero) without rain
    rain = microphysics.rain_flux_coefficient *
           rho_rain^convert(RealT, 1.1364) *
           sqrt(microphysics.rho_surface / rho_dry)
    return rain, SVector(zero(RealT), zero(RealT), rain)
end

@inline function max_terminal_velocity(moisture, rho_dry,
                                       microphysics::KesslerMicrophysics)
    return terminal_velocity_rain(rho_dry, moisture[3], microphysics)
end

@doc raw"""
    saturation_mixing_ratio_tetens(T, p)

Saturation mixing ratio of water vapor of Tetens as in Klemp and Wilhelmson
(1978, Eq. 2.11) and DCMIP2016, with the temperature `T` in K and the pressure `p` in Pa,
```math
q_{vs}(T, p) = \frac{380}{p} \exp\left(17.27 \frac{T - 273}{T - 36}\right).
```
"""
@inline function saturation_mixing_ratio_tetens(T, p)
    RealT = typeof(T)
    return convert(RealT, 380) / p * exp(convert(RealT, 17.27) * (T - 273) / (T - 36))
end

@inline function microphysics_tendencies(moisture, rho_dry, T, p, c_p,
                                         microphysics::KesslerMicrophysics)
    RealT = eltype(moisture)
    zero_ = zero(RealT)
    q_v = max(moisture[1], zero_) / rho_dry
    q_c = max(moisture[2], zero_) / rho_dry
    q_r = max(moisture[3], zero_) / rho_dry
    latent_heat_over_c_p = microphysics.latent_heat / c_p

    # Condensation, where the denominator accounts for the increase of the saturation
    # mixing ratio due to the latent heating with the derivative
    # `d q_vs / dT = q_vs * 17.27 * 237 / (T - 36)^2` of `saturation_mixing_ratio_tetens`
    q_vs = saturation_mixing_ratio_tetens(T, p)
    condensation = max((q_v - q_vs) /
                       (1 +
                        q_vs * convert(RealT, 237 * 17.27) * latent_heat_over_c_p /
                        (T - 36)^2), -q_c) / microphysics.condensation_time

    # Autoconversion and accretion
    autoconversion = convert(RealT, 1.0e-3) * max(q_c - convert(RealT, 1.0e-3), zero_)
    accretion = convert(RealT, 2.2) * q_c * q_r^convert(RealT, 0.875)

    # Evaporation of rain
    r = convert(RealT, 1.0e-3) * rho_dry
    rain = r * q_r
    evaporation = ((convert(RealT, 1.6) +
                    convert(RealT, 124.9) * rain^convert(RealT, 0.2046)) *
                   rain^convert(RealT, 0.525)) /
                  (convert(RealT, 2.55e6) / (convert(RealT, 1.0e-2) * p * q_vs) +
                   convert(RealT, 5.4e5)) *
                  max(q_vs - q_v, zero_) / (r * q_vs)

    heating = latent_heat_over_c_p * (condensation - evaporation)
    return heating,
           SVector(rho_dry * (evaporation - condensation),
                   rho_dry * (condensation - autoconversion - accretion),
                   rho_dry * (autoconversion + accretion - evaporation))
end
end # @muladd
