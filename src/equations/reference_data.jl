@muladd begin
#! format: noindent

# Physical constants in SI units, with values taken from the Williamson et al. test suite 
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
    axis_cross_x_0 = Trixi.cross(axis, x_0)
    x_0 = x_0 + sin(omega * t) * axis_cross_x_0 +
          (1 - cos(omega * t)) * Trixi.cross(axis, axis_cross_x_0)

    # compute Gaussian bump profile
    h = h_0 *
        exp(-b_0 * ((x[1] - x_0[1])^2 + (x[2] - x_0[2])^2 + (x[3] - x_0[3])^2) / (a^2))

    # get Cartesian velocity components
    vx, vy, vz = omega * Trixi.cross(axis, x)

    # Convert primitive variables from Cartesian coordinates to the chosen global 
    # coordinate system, which depends on the equation type
    return cartesian2global(SVector(h, vx, vy, vz, 0.0f0), x, equations)
end

@doc raw"""
    initial_condition_geostrophic_balance(x, t, equations)

Steady geostrophic balance for the spherical shallow water equations, corresponding to a 
purely zonal velocity field given as a function of the latitude $\theta$ by 
$v_\lambda(\theta) = v_0 \cos\theta$, where we define 
$v_0 = 2\pi a \cos(\theta) / 12 \ \mathrm{days}$ in terms of the Earth's radius 
$a = 6.37122 \times 10^3\ \mathrm{m}$. The height field then varies with the latitude as
```math
h(\theta) = 1/g \Big(gh_0 - \Big(a \omega v_0 + \frac{1}{2} v_0^2\Big)\sin^2\theta\Big),
```
where $gh_0 = 2.94 \times 10^4 \ \mathrm{m}^2/\mathrm{s}^2$, 
$g = 9.80616 \ \mathrm{m}/\mathrm{s}^2$, and $\omega = 7.292e-5 \mathrm{s}^{-1}$. This 
problem corresponds to Case 2 of the test suite described in the following paper:
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

Rossby-Haurwitz wave with wave number 4, corresponding to Case 5 of the test suite 
described in the following paper:
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
    K = 7.848f-6
    R = 4.0f0

    A = 0.5f0 * K * (2 * EARTH_ROTATION_RATE + K) * (cos(lat))^2 +
        0.25f0 * K^2 * (cos(lat))^(2 * R) *
        ((R + 1) * (cos(lat))^2 + (2 * R^2 - R - 2) - 2 * R^2 / ((cos(lat))^2))
    B = 2 * (EARTH_ROTATION_RATE + K) * K / ((R + 1) * (R + 2)) * (cos(lat))^R *
        ((R^2 + 2R + 2) - (R + 1)^2 * (cos(lat))^2)
    C = 0.25f0 * K^2 * (cos(lat))^(2 * R) * ((R + 1) * (cos(lat))^2 - (R + 2))

    # compute geopotential height
    h = h_0 +
        (1 / EARTH_GRAVITATIONAL_ACCELERATION) *
        (a^2 * A + a^2 * B * cos(R * lon) + a^2 * C * cos(2 * R * lon))

    # compute zonal and meridional components of the velocity
    vlon = a * K * cos(lat) +
           a * K * (cos(lat))^(R - 1) * (R * (sin(lat))^2 -
                                         (cos(lat))^2) * cos(R * lon)
    vlat = -a * K * R * (cos(lat))^(R - 1) * sin(lat) * sin(R * lon)

    # Convert primitive variables from spherical coordinates to the chosen global 
    # coordinate system, which depends on the equation type
    return spherical2global(SVector(h, vlon, vlat, zero(RealT)), x, equations)
end
end # muladd
