"""
    cartesian_to_spherical_coordinates(cartesian_coord)

# Input
    cartesian_coord: vector of Cartesian coordinates (x, y, z)

# Output
    vector of spherical coordinates, i.e.
    1. longitude (phi) in [-pi, pi]
       measured from x-axis (phi=0) with positive values in eastwards direction
    2. latitude (theta) in [-pi/2, pi/2]
       measuring elevation from x-y-plane (theta=0), i.e. north pole is at pi/2
    3. radius
"""
function cartesian_to_spherical_coordinates(cartesian_coord)
    r = norm(cartesian_coord)
    lon = atan(cartesian_coord[2], cartesian_coord[1])
    lat = asin(cartesian_coord[3] / r)
    return SVector(lon, lat, r)
end

"""
    spherical_to_cartesian_velocities(spherical_coord, spherical_velocities)

# Input
    spherical_coord: position in spherical coordinates lon, lat, r (see above)
    spherical_velocities: velocity components u (azimutal), v (meridional), w (radial)

# Output
    vector of Cartesian velocity components
"""
function spherical_to_cartesian_velocities(spherical_coord, spherical_velocities)
    phi, theta, _ = spherical_coord
    u, v, w = spherical_velocities

    v1 = cos(theta) * cos(phi) * w - sin(theta) * cos(phi) * v - sin(phi) * u
    v2 = cos(theta) * sin(phi) * w - sin(theta) * sin(phi) * v + cos(phi) * u
    v3 = sin(theta) * w + cos(theta) * v

    return SVector(v1, v2, v3)
end
