"""
    cartesian_to_spherical_coordinates(cartesian_coord)

# Input
    x: vector of Cartesian coordinates

# Output
    vector of spherical cordinates, i.e.
    1. longitude phi in [-pi, pi] 
       measured from x-axis (phi=0) with positive values in eastwards direction
    2. latitude theta in [-pi/2, pi/2]
       measuring elevation from x-y-plane (theta=0), i.e. north pole is at pi/2
    3. radius
"""
function cartesian_to_spherical_coordinates(cartesian_coord)
    r = norm(cartesian_coord)
    theta = atan(cartesian_coord[2], cartesian_coord[1])
    phi = asin(cartesian_coord[3] / r)
    return SVector(phi, theta, r)
end


"""
    spherical_to_cartesian_velocities(spherical_coord, spherical_velocities)

# Input
    coord: position in spherical coordinates r, theta, phi (see above)
    vel: velocity components w (radial), v (meridional), u (azimutal)

# Output
    vector of Cartesian velocity components
"""
function spherical_to_cartesian_velocities(spherical_coord, spherical_velocities)
    phi, theta, _ = spherical_coord
    u, v, w = spherical_velocities

    v1 = cos(theta)*cos(phi) * w -sin(theta)*cos(phi) * v -sin(phi) * u
    v2 = cos(theta)*sin(phi) * w -sin(theta)*sin(phi) * v +cos(phi) * u
    v3 = sin(theta)          * w +cos(theta)          * v 

    return SVector(v1, v2, v3)
end
