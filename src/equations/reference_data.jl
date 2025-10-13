@muladd begin
#! format: noindent

# Physical constants in SI units (reference values from the Williamson et al. test suite)
const EARTH_RADIUS = 6.371229 # 6.37122e6  # m
const EARTH_GRAVITATIONAL_ACCELERATION = 9.81 # 9.80616  # m/s²
const EARTH_ROTATION_RATE = 7.29212e-5 # 7.292e-5  # rad/s
const SECONDS_PER_DAY = 8.64e4

@doc raw"""
    initial_condition_gaussian(x, t, equations)

This Gaussian bell case is a smooth initial condition suitable for testing the convergence 
of discretizations of the linear advection equation on a spherical domain of radius $a = 6.
37122 \times 10^3\ \mathrm{m}$, representing the surface of the Earth. Denoting the 
Euclidean norm as $\lVert \cdot \rVert$, the initial height field is prescribed as a 
function of the position $\vec{x}$ relative to the centre of the Earth by
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
    return cartesian2global(SVector(h, vx, vy, vz, zero(RealT)), x, equations)
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
    return spherical2global(SVector(h, vlon, vlat, zero(RealT), zero(RealT)), x,
                            equations)
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
A(\theta) &= \frac{\omega}{2}(2 \Omega+\omega) \cos^2 \theta + 
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
    a = sqrt(x[1]^2 + x[2]^2 + x[3]^2)  # radius of the sphere
    lon, lat = atan(x[2], x[1]), asin(x[3] / a)

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
    return spherical2global(SVector(h, vlon, vlat, zero(RealT), zero(RealT)), x,
                            equations)
end

@doc raw"""
    initial_condition_isolated_mountain(x, t, equations)

Zonal flow over an isolated mountain with a profile given in terms of the latitude 
$\lambda$ and longitude $\theta$ as 
```math
h_s(\lambda,\theta) =
h_{s0} (1 - \sqrt{\min(R^2, (\lambda-\lambda_0)^2 + (\theta-\theta_0)^2)}/R),
```
where $h_{s0} = 2000 \ \text{m}$, $\lambda_0 = -\pi/2$, $\theta_0 = \pi/6$, and $R =\pi/9$. 
The initial velocity field is given by $v_\lambda(\theta) = v_0 \cos\theta$, where 
$v_0 = 20 \ \mathrm{m/s}$, and the total geopotential height $H = h+h_s$ is given by 
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
    return spherical2global(SVector(h, vlon, vlat, zero(RealT),
                                    bottom_topography_isolated_mountain(x)), x,
                            equations)
end

# Bottom topography function to pass as auxiliary_field keyword argument in constructor for 
# SemidiscretizationHyperbolic, used with initial_condition_isolated_mountain
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

Unsteady analytical solution to the spherical shallow water equations, corresponding to a
solid body rotation with a prescribed bottom topography. Assuming the domain to be a sphere 
of radius $a = 6.37122 \times 10^3\ \mathrm{m}$, letting $\vec{x}$ denote the position 
relative to the centre of the Earth, and letting $\vec{e}_x$, $\vec{e}_y$, and $\vec{e}_z$ 
denote the Cartesian basis vectors, we define the rotating frame 
```math
\vec{b}_x(t) = \cos(\Omega t)\vec{e}_x + \sin(\Omega t)\vec{e}_y, \quad
\vec{b}_y(t) = -\sin(\Omega t)\vec{e}_x + \cos(\Omega t)\vec{e}_y, \quad 
\vec{b}_z(t) = \vec{e}_z,
```
as a function of time $t$. We also define the associated coordinate transformation 
```math
\vec{\varphi}(\vec{x},t) = 
(\vec{x} \cdot \vec{b}_x(t)) \vec{e}_x + (\vec{x} \cdot \vec{b}_y(t)) \vec{e}_y + 
(\vec{x} \cdot \vec{b}_z(t)) \vec{e}_z
```
as well as a fixed axis $\vec{c} = -\sin(\alpha)\vec{e}_x + \cos(\alpha)\vec{e}_y$, where 
$\Omega = 7.292 \times 10^{-5}\ \mathrm{s}^{-1}$ is the Earth's rotation rate, and we take 
$\alpha = \pi/4$. For a bottom topography prescribed as
```math
h_s(\vec{x}) = \frac{1}{2g}(\vec{\Omega} \cdot \vec{x})^2
```
where $\vec{\Omega} = \Omega\vec{e}_z$ and $g = 9.80616 \ \mathrm{m}/\mathrm{s}^2$ are the 
Earth's axis of rotation and gravitational acceleration, respectively, the time-dependent 
velocity field is given as 
```math
\vec{v}(\vec{x},t) = v_0\, \vec{\varphi}(\vec{c},t) \times \vec{x}/\lVert \vec{x} \rVert,
```
and the total geopotential height $H = h+b$ is given by 
```math
H(\vec{x},t) = 
\frac{1}{2g}\left(\big(v_0 \,\vec{\Omega}\cdot \vec{x} - 
\vec{\varphi}(\vec{c},t) \cdot \vec{x}/\lVert \vec{x} \rVert \big)^2 +
(\vec{\Omega} \cdot \vec{x})^2 + 2k_1\right),
```
where $v_0 = 2\pi a / (12 \ \mathrm{days})$, $k_1 = 133681 \ \mathrm{m}^2/\mathrm{s}^2$,
and $\lVert \cdot \rVert$ denotes the Euclidean norm. To use this test case with
[`SplitCovariantShallowWaterEquations2D`](@ref), the keyword argument
`auxiliary_field = bottom_topography_unsteady_solid_body_rotation` should be passed into
the `SemidiscretizationHyperbolic` constructor. This analytical solution was derived in the
following paper:
- M. Läuter, D. Handorf, and K. Dethloff (2005). Unsteady analytical solutions of the 
  spherical shallow water equations. Journal of Computational Physics 210:535–553.
  [DOI: 10.1016/j.jcp.2005.04.022](https://doi.org/10.1016/j.jcp.2005.04.022)
"""
@inline function initial_condition_unsteady_solid_body_rotation(x, t, equations)
    RealT = eltype(x)
    a = sqrt(x[1]^2 + x[2]^2 + x[3]^2)  # radius of the sphere

    # parameters
    v_0 = convert(RealT, 2π) * a / (12 * SECONDS_PER_DAY)
    k1 = 133681.0f0
    alpha = convert(RealT, π / 4)

    # rotating frame
    b1 = SVector(cos(EARTH_ROTATION_RATE * t), sin(EARTH_ROTATION_RATE * t), 0.0f0)
    b2 = SVector(-sin(EARTH_ROTATION_RATE * t), cos(EARTH_ROTATION_RATE * t), 0.0f0)
    b3 = SVector(0.0f0, 0.0f0, 1.0f0)

    # axis and normal vectors 
    c = SVector(-sin(alpha), 0.0f0, cos(alpha))
    n = x / norm(x)

    # dot product of Earth's rotation axis and position vector
    Omega_dot_x = EARTH_ROTATION_RATE * x[3]

    # coordinate transformation
    phi_t = SVector(dot(c, b1), dot(c, b2), dot(c, b3))

    # compute velocity and total geopotential height
    v = v_0 * cross(phi_t, n)
    H = (-0.5f0 * (v_0 * dot(phi_t, n) + Omega_dot_x)^2 +
         0.5f0 * Omega_dot_x^2 + k1) / EARTH_GRAVITATIONAL_ACCELERATION

    # Convert primitive variables from Cartesian coordinates to the chosen global 
    # coordinate system, which depends on the equation type
    return cartesian2global(SVector(H, v[1], v[2], v[3],
                                    bottom_topography_unsteady_solid_body_rotation(x)),
                            x, equations)
end

# Bottom topography function to pass as auxiliary_field keyword argument in constructor for 
# SemidiscretizationHyperbolic, used with initial_condition_unsteady_solid_body_rotation
@inline function bottom_topography_unsteady_solid_body_rotation(x)
    return 0.5f0 * (EARTH_ROTATION_RATE * x[3])^2 / EARTH_GRAVITATIONAL_ACCELERATION
end

@doc raw"""
    initial_condition_barotropic_instability(x, t, equations)

Barotrotropic instability initiated by a perturbation applied to a mid-latitude jet. The  velocity field is a purely zonal flow, given as a function of the latitude $\theta$ as
```math
v_\lambda(\theta) = \begin{cases}
u_0 \exp(-4 / (\theta_1 - \theta_0)^{-2})\exp((\theta - \theta_0)^{-1}
(\theta - \theta_1)^{-1}), & \quad \theta_0 < \theta < \theta_1, \\
0 & \quad \text{otherwise},
\end{cases}
```
where $u_0 = 80 \ \mathrm{m}/\mathrm{s}$, $\theta_0 = \pi/7$, and 
$\theta_1 = \pi/2 - \theta_0$. The background geopotential height field is given by 
```math
h_0(\theta) = 10158 \ \mathrm{m} - 
\frac{a}{g} \int_{-\pi/2}^\theta v_\lambda(\theta')\big(2\Omega\sin\theta' + 
v_\lambda(\theta')\tan\theta' / a \big)\, \mathrm{d}\theta',
``` 
where $a = 6.37122 \times 10^3\ \mathrm{m}$ is the Earth's radius, 
$g = 9.80616 \ \mathrm{m}/\mathrm{s}^2$ is the Earth's gravitational acceleration, and
$\Omega = 7.292 \times 10^{-5}\ \mathrm{s}^{-1}$ is the Earth's rotation rate. The 
perturbation is then added to obtain the following geopotential height field:
```math
h(\lambda, \theta) = \begin{cases}
h_0(\theta) + \delta h \cos\theta \exp(-(\lambda/\alpha)^2) 
\exp(-((\theta_2 -\theta)/\beta)^2), & \quad -\pi < \lambda < \pi,\\
h_0(\theta), & \quad \text{otherwise},
\end{cases}
```
where $\lambda$ is the longitude coordinate, and we take $\alpha = 1/3$, $\beta = 1/15$, 
and $\delta h = 120 \ \mathrm{m}$. This problem was proposed in the following paper:
- J. Galewsky, R. K. Scott, and L. M. Polvani (2004). An initial-value problem for
  testing numerical models of the global shallow-water equations. Tellus A 56.5:429–440.
  [DOI: 10.3402/tellusa.v56i5.14436](https://doi.org/10.3402/tellusa.v56i5.14436)
"""
@inline function initial_condition_barotropic_instability(x, t, equations)
    RealT = eltype(x)
    a = sqrt(x[1]^2 + x[2]^2 + x[3]^2)  # radius of the sphere
    lon, lat = atan(x[2], x[1]), asin(x[3] / a)

    # compute zonal and meridional velocity components
    u_0 = 80.0f0
    lat_0 = convert(RealT, π / 7)
    lat_1 = convert(RealT, π / 2) - lat_0
    vlon = galewsky_velocity(lat, u_0, lat_0, lat_1)
    vlat = zero(eltype(x))

    # numerically integrate (here we use the QuadGK package) to get height
    galewsky_integral, _ = quadgk(latp -> galewsky_integrand(latp, u_0, lat_0, lat_1, a),
                                  convert(RealT, π / 2), lat)
    h = 10158.0f0 - a / EARTH_GRAVITATIONAL_ACCELERATION * galewsky_integral

    # add perturbation to initiate instability
    alpha, beta = convert(RealT, 1 / 3), convert(RealT, 1 / 15)
    lat_2 = convert(RealT, π / 4)
    if (-π < lon) && (lon < π)
        h = h +
            120.0f0 * cos(lat) *
            exp(-((lon / alpha)^2)) * exp(-((lat_2 - lat) / beta)^2)
    end

    # Convert primitive variables from spherical coordinates to the chosen global 
    # coordinate system, which depends on the equation type
    return spherical2global(SVector(h, vlon, vlat, zero(RealT), zero(RealT)), x,
                            equations)
end

@inline function galewsky_velocity(lat, u_0, lat_0, lat_1)
    if (lat_0 < lat) && (lat < lat_1)
        u = u_0 / exp(-4 / (lat_1 - lat_0)^2) *
            exp(1 / (lat - lat_0) * 1 / (lat - lat_1))
    else
        u = zero(lat)
    end
    return u
end

@inline function galewsky_integrand(lat, u_0, lat_0, lat_1, a)
    u = galewsky_velocity(lat, u_0, lat_0, lat_1)
    return u * (2 * EARTH_ROTATION_RATE * sin(lat) + u * tan(lat) / a)
end
end # @muladd
