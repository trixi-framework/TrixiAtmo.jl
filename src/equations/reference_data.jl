@muladd begin
#! format: noindent

# Physical constants
const EARTH_RADIUS = 6.37122e6   # m
const EARTH_GRAVITATIONAL_ACCELERATION = 9.80616  # m/s²
const EARTH_ROTATION_RATE = 7.292e-5  # rad/s
const SECONDS_PER_DAY = 8.64e4

@doc raw"""
    initial_condition_gaussian(x, t, equations)

This Gaussian bell case is a smooth initial condition suitable for testing the convergence 
of discretizations of the linear advection equation on a spherical domain of radius $6.
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
@inline function initial_condition_gaussian(x, t, ::AbstractEquations)
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

    # Prescribe the rotated bell shape and Cartesian velocity components.
    # The last variable is the bottom topography, which we set to zero
    return SVector(h, vx, vy, vz, 0.0f0)
end

# Version for spherical coordinates (note: the velocity is not well defined at the poles)
@inline function initial_condition_gaussian(x, t,
                                            ::AbstractCovariantEquations{2, 3,
                                                                         GlobalSphericalCoordinates})
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

    # get zonal and meridional components of the velocity
    lon, lat = atan(x[2], x[1]), asin(x[3] / a)
    vlon = omega * a * (cos(lat) * cos(alpha) + sin(lat) * cos(lon) * sin(alpha))
    vlat = -omega * a * sin(lon) * sin(alpha)

    # Prescribe the rotated bell shape and spherical velocity components.
    # The last variable is the bottom topography, which we set to zero
    return SVector(h, vlon, vlat, 0.0f0, 0.0f0)
end
end # muladd
