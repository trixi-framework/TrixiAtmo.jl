@muladd begin
#! format: noindent

# Physical constants
const EARTH_RADIUS = 6.37122e6   # m
const EARTH_GRAVITATIONAL_ACCELERATION = 9.80616  # m/s²
const EARTH_ROTATION_RATE = 7.292e-5  # rad/s
const SECONDS_PER_DAY = 8.64e4

@doc raw"""
    initial_condition_gaussian(x, t)

This Gaussian bell case is a smooth initial condition suitable for testing the convergence 
of discretizations of the linear advection equation on a spherical domain of radius $a = 6.
37122 \times 10^3\ \mathrm{m}$, representing the surface of the Earth. The height field is 
given in terms of the Cartesian coordinates $x$, $y$, $z$ as 
```math
h(x,y,z) = h_0 \exp\Big(-b_0 \big((x-x_0)^2 + (y-y_0)^2 + (z-z_0)^2\big) \Big),
```
where $h_0 = 1 \times 10^3\ \mathrm{m}$ is the height of the bell, $b_0 = 5 / a$ is the 
width parameter, and the Cartesian coordinates of the bell's centre at longitude 
$\lambda_0 = 3\pi/2$ and latitude $\theta_0 = 0$ are given by 
```math
\begin{aligned}
x_0 &= a\cos\theta_0\cos\lambda_0,\\
y_0 &= a\cos\theta_0\sin\lambda_0,\\
z_0 &= a\sin\theta_0.
\end{aligned}
```
The velocity field corresponds to a solid body rotation with a period of 12 days, at an 
angle of $\alpha = \pi/4$ from the polar axis. The zonal and meridional components of the 
velocity field are then given in terms of the longitude $\lambda$ and latitude $\theta$ as 
```math
\begin{aligned}
u(\lambda,\theta) &= V (\cos\theta\cos\alpha + \sin\theta\cos\lambda\sin\alpha),\\
v(\lambda,\theta)  &= -V \sin\lambda\sin\alpha,
\end{aligned}
```
where we take $V = 2\pi a / 12\ \mathrm{days}$. This test case is adapted from Case 1 of 
the test suite described in the following paper:
- D. L. Williamson, J. B. Drake, J. J. Hack, R. Jakob, and P. N. Swarztrauber (1992). A  
  standard test set for numerical approximations to the shallow water equations in
  spherical geometry. Journal of Computational Physics, 102(1):211-224. 
  [DOI: 10.1016/S0021-9991(05)80016-6](https://doi.org/10.1016/S0021-9991(05)80016-6)

This function returns `SVector(h, vlon, vlat, b)`, where the first three entries are the 
height, zonal velocity, and meridional velocity. The fourth entry, representing variable 
bottom topography, is set to zero. The functions [`transform_to_contravariant`](@ref) and 
[`transform_to_cartesian`](@ref) are available for converting to the prognostic variables 
for Cartesian and covariant formulations.
"""
@inline function initial_condition_gaussian(x, t)
    RealT = eltype(x)

    a = sqrt(x[1]^2 + x[2]^2 + x[3]^2) #radius of the sphere
    V = convert(RealT, 2π) * a / (12 * SECONDS_PER_DAY)  # speed of rotation
    alpha = convert(RealT, π / 4)  # angle of rotation
    h_0 = 1000.0f0  # bump height in metres
    b_0 = 5.0f0 / (a^2)  # bump width
    lon_0, lat_0 = convert(RealT, 3π / 2), 0.0f0  # initial bump location

    # convert initial position to Cartesian coordinates
    x_0 = SVector(a * cos(lat_0) * cos(lon_0),
                  a * cos(lat_0) * sin(lon_0),
                  a * sin(lat_0))

    # compute Gaussian bump profile
    h = h_0 * exp(-b_0 * ((x[1] - x_0[1])^2 + (x[2] - x_0[2])^2 + (x[3] - x_0[3])^2))

    # get zonal and meridional components of the velocity
    lon, lat = atan(x[2], x[1]), asin(x[3] / a)
    vlon = V * (cos(lat) * cos(alpha) + sin(lat) * cos(lon) * sin(alpha))
    vlat = -V * sin(lon) * sin(alpha)

    # the last variable is the bottom topography, which we set to zero
    return SVector(h, vlon, vlat, 0.0)
end
end # muladd
