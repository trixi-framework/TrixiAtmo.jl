@muladd begin
#! format: noindent

# Physical constants in SI units (reference values from the Williamson et al. test suite)
const EARTH_RADIUS = 6.37122e6  # m
const EARTH_GRAVITATIONAL_ACCELERATION = 9.80616  # m/s²
const EARTH_ROTATION_RATE = 7.292e-5  # rad/s
const SECONDS_PER_DAY = 8.64e4

@doc raw"""
    initial_condition_gaussian(x, t, equations)

This Gaussian bell case is a smooth initial condition suitable for testing the convergence 
of discretizations of the linear advection equation on a spherical domain of radius $a = 6.
37122 \times 10^3\ \mathrm{m}$, representing the surface of the Earth. Denoting the 
Euclidean norm as $\lVert \cdot \rVert$, the initial height field is given by
```math
h(\vec{x}) = h_0 \exp
\Big(-b_0 \big(\lVert \vec{x} - \vec{x}_0 \rVert / \lVert \vec{x} \rVert\big)^2 \Big),
```
where $h_0 = 1 \times 10^3\ \mathrm{m}$ is the height of the bell, $b_0 = 5$ is the 
width parameter, and $\vec{x}_0$ is the position of the centre of the bell, which is 
initialized at a longitude of $3\pi/2$ and a latitude of zero. The velocity field 
corresponds to a solid body rotation with a period of 12 days at an angle of 
$\alpha = \pi/4$ from the polar axis. Denoting $\vec{\omega}$ as the corresponding angular
velocity vector, the velocity is therefore initialized as
```math
\vec{v}(\vec{x}) = \vec{\omega} \times \vec{x}.
```
This problem is adapted from Case 1 of the test suite described in the following paper:
- D. L. Williamson, J. B. Drake, J. J. Hack, R. Jakob, and P. N. Swarztrauber (1992). A  
  standard test set for numerical approximations to the shallow water equations in
  spherical geometry. Journal of Computational Physics, 102(1):211-224. 
  [DOI: 10.1016/S0021-9991(05)80016-6](https://doi.org/10.1016/S0021-9991(05)80016-6)
"""
@inline function initial_condition_gaussian(x, t, equations)
    RealT = eltype(x)
    a = sqrt(x[1]^2 + x[2]^2 + x[3]^2)  # radius of the sphere
    omega = convert(RealT, 2π) / (12 * SECONDS_PER_DAY)  # angular velocity
    alpha = convert(RealT, π / 4)  # angle of rotation
    h_0 = 1000.0f0  # bump height in metres
    b_0 = 5.0f0  # bump width parameter
    lon_0, lat_0 = convert(RealT, 3π / 2), 0.0f0  # initial bump location

    # axis of rotation
    axis = SVector(-cos(alpha), 0.0f0, sin(alpha))

    # convert initial position to Cartesian coordinates
    x_0 = SVector(a * cos(lat_0) * cos(lon_0),
                  a * cos(lat_0) * sin(lon_0),
                  a * sin(lat_0))

    # apply rotation using Rodrigues' formula 
    axis_cross_x_0 = cross(axis, x_0)
    x_0 = x_0 + sin(omega * t) * axis_cross_x_0 +
          (1 - cos(omega * t)) * cross(axis, axis_cross_x_0)

    # compute Gaussian bump profile
    h = h_0 *
        exp(-b_0 * ((x[1] - x_0[1])^2 + (x[2] - x_0[2])^2 + (x[3] - x_0[3])^2) / (a^2))

    # get Cartesian velocity components
    vx, vy, vz = omega * cross(axis, x)

    # Convert primitive variables from Cartesian coordinates to the chosen global 
    # coordinate system, which depends on the equation type
    return cartesian2global(SVector(h, vx, vy, vz, 0.0f0), x, equations)
end

@doc raw"""
    initial_condition_geostrophic_balance(x, t, equations)

Steady geostrophic balance for the spherical shallow water equations, corresponding to a 
purely zonal velocity field given as a function of the latitude $\theta$ by 
$v_\lambda(\theta) = v_0 \cos\theta$, where we define $v_0 = 2\pi a / (12 \ \mathrm{days})$ 
in terms of the Earth's radius $a = 6.37122 \times 10^3\ \mathrm{m}$. The height field 
then varies with the latitude as
```math
h(\theta) = \frac{1}{g} 
\Big(gh_0 - \Big(a \Omega v_0 + \frac{1}{2} v_0^2\Big)\sin^2\theta\Big),
```
where $gh_0 = 2.94 \times 10^4 \ \mathrm{m}^2/\mathrm{s}^2$, 
$g = 9.80616 \ \mathrm{m}/\mathrm{s}^2$, and 
$\Omega = 7.292 \times 10^{-5}\ \mathrm{s}^{-1}$. This problem corresponds to Case 2 of the
test suite described in the following paper:
- D. L. Williamson, J. B. Drake, J. J. Hack, R. Jakob, and P. N. Swarztrauber (1992). A  
  standard test set for numerical approximations to the shallow water equations in
  spherical geometry. Journal of Computational Physics, 102(1):211-224. 
  [DOI: 10.1016/S0021-9991(05)80016-6](https://doi.org/10.1016/S0021-9991(05)80016-6)
"""
@inline function initial_condition_geostrophic_balance(x, t, equations)
    RealT = eltype(x)
    a = sqrt(x[1]^2 + x[2]^2 + x[3]^2)  # radius of the sphere
    lat = asin(x[3] / a)

    # compute zonal and meridional components of the velocity
    v_0 = convert(RealT, 2π) * a / (12 * SECONDS_PER_DAY)
    vlon, vlat = v_0 * cos(lat), zero(eltype(x))

    # compute geopotential height 
    h = 1 / EARTH_GRAVITATIONAL_ACCELERATION *
        (2.94f4 - (a * EARTH_ROTATION_RATE * v_0 + 0.5f0 * v_0^2) * (sin(lat))^2)

    # Convert primitive variables from spherical coordinates to the chosen global 
    # coordinate system, which depends on the equation type
    return spherical2global(SVector(h, vlon, vlat, zero(RealT)), x, equations)
end

@doc raw"""
    initial_condition_rossby_haurwitz(x, t, equations)

Rossby-Haurwitz wave case for the spherical shallow water equations, where the zonal and 
meridional velocity components are given, respectively, as functions of the longitude 
$\lambda$ and latitude $\theta$ by
```math
\begin{aligned}
v_\lambda(\lambda,\theta) &= a \omega \cos \theta+a K \cos ^{R-1} \theta
\left(R \sin ^2 \theta-\cos ^2 \theta\right) \cos (R \lambda),\\
v_\theta(\lambda,\theta) &= -a K R \cos ^{R-1} \theta \sin \theta \sin (R \lambda),
\end{aligned}
```
where $\omega = K = 7.848 \times 10^{-6} \ \mathrm{s}^{-1}$ and $R = 4$ are given 
constants, and $a = 6.37122 \times 10^3\ \mathrm{m}$ is the Earth's radius. Taking 
$g = 9.80616 \ \mathrm{m}/\mathrm{s}^2$, $\Omega = 7.292 \times 10^{-5} \ \mathrm{s}^{-1}$, 
and $h_0 = 8000 \ \mathrm{m}$ and defining the functions 
```math
\begin{aligned}
A(\theta) &=  \frac{\omega}{2}(2 \Omega+\omega) \cos^2 \theta + 
\frac{1}{4} K^2 \cos^{2 R} \theta\Big((R+1) \cos^2\theta +\left(2 R^2-R-2\right) - 
\big(2 R^2 / \cos^2 \theta\big) \Big), \\
B(\theta) &= \frac{2(\Omega+\omega) K}{(R+1)(R+2)} \cos ^R \theta\big((R^2+2 R+2) - 
(R+1)^2 \cos^2 \theta\big), \\
C(\theta) &=  \frac{1}{4} K^2 \cos^{2 R} \theta\big((R+1) \cos^2 \theta-(R+2)\big),
\end{aligned}
```
the initial height field is given by
```math
h(\lambda,\theta) = h_0 + 
\frac{a^2}{g}\Big(A(\theta) + B(\theta)\cos(R\lambda) + C(\theta)\cos(2R\lambda) \Big).
```
This problem corresponds to Case 6 of the test suite described in the following paper:
- D. L. Williamson, J. B. Drake, J. J. Hack, R. Jakob, and P. N. Swarztrauber (1992). A  
  standard test set for numerical approximations to the shallow water equations in
  spherical geometry. Journal of Computational Physics, 102(1):211-224. 
  [DOI: 10.1016/S0021-9991(05)80016-6](https://doi.org/10.1016/S0021-9991(05)80016-6)
"""
@inline function initial_condition_rossby_haurwitz(x, t, equations)
    RealT = eltype(x)
    a = sqrt(x[1]^2 + x[2]^2 + x[3]^2)
    lon = atan(x[2], x[1])
    lat = asin(x[3] / a)

    h_0 = 8.0f3
    omega = 7.848f-6
    K = 7.848f-6
    R = 4.0f0

    A = 0.5f0 * omega * (2 * EARTH_ROTATION_RATE + omega) * (cos(lat))^2 +
        0.25f0 * K^2 * (cos(lat))^(2 * R) *
        ((R + 1) * (cos(lat))^2 + (2 * R^2 - R - 2) - 2 * R^2 / ((cos(lat))^2))
    B = 2 * (EARTH_ROTATION_RATE + omega) * K / ((R + 1) * (R + 2)) * (cos(lat))^R *
        ((R^2 + 2R + 2) - (R + 1)^2 * (cos(lat))^2)
    C = 0.25f0 * K^2 * (cos(lat))^(2 * R) * ((R + 1) * (cos(lat))^2 - (R + 2))

    # compute geopotential height
    h = h_0 +
        (1 / EARTH_GRAVITATIONAL_ACCELERATION) *
        (a^2 * A + a^2 * B * cos(R * lon) + a^2 * C * cos(2 * R * lon))

    # compute zonal and meridional components of the velocity
    vlon = a * omega * cos(lat) +
           a * K * (cos(lat))^(R - 1) * (R * (sin(lat))^2 -
                                         (cos(lat))^2) * cos(R * lon)
    vlat = -a * K * R * (cos(lat))^(R - 1) * sin(lat) * sin(R * lon)

    # Convert primitive variables from spherical coordinates to the chosen global 
    # coordinate system, which depends on the equation type
    return spherical2global(SVector(h, vlon, vlat, zero(RealT)), x, equations)
end

@doc raw"""
    initial_condition_isolated_mountain(x, t, equations)

Zonal flow over an isolated mountain with a profile given in terms of the latitude 
$\lambda$ and longitude $\theta$ as 
```math
b(\lambda,\theta) = 
b_0 (1 - \sqrt{\min(R^2, (\lambda-\lambda_0)^2 + (\theta-\theta_0)^2)}/R),
```
where $b_0 = 2000 \ \text{m}$, $\lambda_0 = -/\pi/2$, $\theta_0 = \pi/6$, and $R =\pi/9$. 
The initial velocity field is given by  $v_\lambda(\theta) = v_0 \cos\theta$, where 
$v_0 = 20 \ \mathrm{m/s}$, and the total height $H = h+b$ is given by 
```math
H(\theta) = H_0 - \frac{1}{g}\Big(a \Omega v_0 + \frac{1}{2} v_0^2\Big)\sin^2\theta,
```
where $H_0 = 5960 \ \mathrm{m}$, $g = 9.80616 \ \mathrm{m}/\mathrm{s}^2$, and 
$\Omega = 7.292 \times 10^{-5}\ \mathrm{s}^{-1}$. To use this test case with 
[`SplitCovariantShallowWaterEquations2D`](@ref), the keyword argument 
`auxiliary_field = bottom_topography_isolated_mountain` should be passed into the 
`SemidiscretizationHyperbolic` constructor. This problem corresponds to Case 5 of the
test suite described in the following paper:
- D. L. Williamson, J. B. Drake, J. J. Hack, R. Jakob, and P. N. Swarztrauber (1992). A  
  standard test set for numerical approximations to the shallow water equations in
  spherical geometry. Journal of Computational Physics, 102(1):211-224. 
  [DOI: 10.1016/S0021-9991(05)80016-6](https://doi.org/10.1016/S0021-9991(05)80016-6)
"""
@inline function initial_condition_isolated_mountain(x, t, equations)
    RealT = eltype(x)
    a = sqrt(x[1]^2 + x[2]^2 + x[3]^2)  # radius of the sphere
    lat = asin(x[3] / a)
    h_0 = 5960.0f0
    v_0 = 20.0f0

    # compute zonal and meridional components of the velocity
    vlon, vlat = v_0 * cos(lat), zero(eltype(x))

    # compute geopotential height 
    h = h_0 -
        1 / EARTH_GRAVITATIONAL_ACCELERATION *
        (a * EARTH_ROTATION_RATE * v_0 + 0.5f0 * v_0^2) * (sin(lat))^2

    # Convert primitive variables from spherical coordinates to the chosen global 
    # coordinate system, which depends on the equation type
    return spherical2global(SVector(h, vlon, vlat, zero(RealT)), x, equations)
end

# Bottom topography function to pass as auxiliary_field keyword argument in constructor for 
# SemidiscretizationHyperbolic, for use with initial_condition_isolated_mountain
@inline function bottom_topography_isolated_mountain(x)
    RealT = eltype(x)
    a = sqrt(x[1]^2 + x[2]^2 + x[3]^2)  # radius of the sphere
    lon, lat = atan(x[2], x[1]), asin(x[3] / a)

    # Position and height of mountain, noting that latitude is λ = -π/2 and not λ = 3π/2 
    # because atan(y,x) is in [-π, π]
    lon_0, lat_0 = convert(RealT, -π / 2), convert(RealT, π / 6)
    b_0 = 2000.0f0

    R = convert(RealT, π / 9)
    return b_0 * (1.0f0 - sqrt(min(R^2, (lon - lon_0)^2 + (lat - lat_0)^2)) / R)
end

@doc raw"""
    initial_condition_unsteady_solid_body_rotation(x, t, equations)

Unsteady solid body rotation for the spherical shallow water equations. This analytical 
solution was derived in the following paper:
- M. Läuter, D. Handorf, and K. Dethloff (2005). Unsteady analytical solutions of the 
  spherical shallow water equations. Journal of Computational Physics 210:535–553.
  [DOI: 10.1016/j.jcp.2005.04.022](https://doi.org/10.1016/j.jcp.2005.04.022)
"""
@inline function initial_condition_unsteady_solid_body_rotation(x, t, equations)
    RealT = eltype(x)
    a = sqrt(x[1]^2 + x[2]^2 + x[3]^2)  # radius of the sphere

    v_0 = convert(RealT, 2π) * a / (12 * SECONDS_PER_DAY)
    k1 = 133681.0f0
    alpha = convert(RealT, π / 4)

    b1 = SVector(cos(EARTH_ROTATION_RATE * t), sin(EARTH_ROTATION_RATE * t), 0.0f0)
    b2 = SVector(-sin(EARTH_ROTATION_RATE * t), cos(EARTH_ROTATION_RATE * t), 0.0f0)
    b3 = SVector(0.0f0, 0.0f0, 1.0f0)

    c = SVector(-sin(alpha), 0.0f0, cos(alpha))
    n = x / norm(x)

    Omega = SVector(0.0f0, 0.0f0, EARTH_ROTATION_RATE)
    phi_t = SVector(dot(c, b1), dot(c, b2), dot(c, b3))

    v = v_0 * cross(phi_t, n)
    h = (-0.5f0 * (v_0 * dot(phi_t, n) + dot(Omega, x))^2 +
         0.5f0 * dot(Omega, x)^2 + k1) / EARTH_GRAVITATIONAL_ACCELERATION

    # Convert primitive variables from Cartesian coordinates to the chosen global 
    # coordinate system, which depends on the equation type
    return cartesian2global(SVector(h, v[1], v[2], v[3]), x, equations)
end

# Bottom topography function to pass as auxiliary_field keyword argument in constructor for 
# SemidiscretizationHyperbolic, for use with initial_condition_unsteady_solid_body_rotation
@inline function bottom_topography_unsteady_solid_body_rotation(x)
    return 0.5f0 * dot(SVector(0.0f0, 0.0f0, EARTH_ROTATION_RATE), x)^2 /
           EARTH_GRAVITATIONAL_ACCELERATION
end
end # @muladd
