@muladd begin
#! format: noindent

# Physical constants
const EARTH_RADIUS = 6.37122  # m
const EARTH_GRAVITATIONAL_ACCELERATION = 9.80616  # m/s²
const EARTH_ROTATION_RATE = 7.292e-5  # rad/s
const SECONDS_PER_DAY = 8.64e4

function spherical2cartesian(vlon, vlat, x)
    # Co-latitude
    colat = acos(x[3] / sqrt(x[1]^2 + x[2]^2 + x[3]^2))

    # Longitude
    if sign(x[2]) == 0.0
        signy = 1.0
    else
        signy = sign(x[2])
    end
    r_xy = sqrt(x[1]^2 + x[2]^2)
    if r_xy == 0.0
        lon = pi / 2
    else
        lon = signy * acos(x[1] / r_xy)
    end

    v1 = -cos(colat) * cos(lon) * vlat - sin(lon) * vlon
    v2 = -cos(colat) * sin(lon) * vlat + cos(lon) * vlon
    v3 = sin(colat) * vlat
    return SVector(v1, v2, v3)
end

function initial_condition_gaussian(x, t)
    RealT = eltype(x)

    # set parameters
    a = EARTH_RADIUS  # radius of the sphere in metres
    V = convert(RealT, 2π) * a / (12 * SECONDS_PER_DAY)  # speed of rotation
    alpha = convert(RealT, π / 4)  # angle of rotation
    h_0 = 1000.0f0  # bump height in metres
    b_0 = 5.0f0 / (a^2)  # bump width
    lon_0, lat_0 = convert(RealT, 3π / 2), 0.0f0  # initial bump location

    # convert initial position to Cartesian coordinates
    x_0 = SVector(a * cos(lat_0) * cos(lon_0), a * cos(lat_0) * sin(lon_0),
                  a * sin(lat_0))

    # compute Gaussian bump profile
    h = h_0 * exp(-b_0 * ((x[1] - x_0[1])^2 + (x[2] - x_0[2])^2 + (x[3] - x_0[3])^2))

    # get zonal and meridional components of the velocity
    lon, lat = atan(x[2], x[1]), asin(x[3] / a)
    vlon = V * (cos(lat) * cos(alpha) + sin(lat) * cos(lon) * sin(alpha))
    vlat = -V * sin(lon) * sin(alpha)

    return SVector(h, vlon, vlat)
end

function initial_condition_gaussian(x, t, equations)
    return initial_condition_gaussian(x, t)
end
end # muladd
