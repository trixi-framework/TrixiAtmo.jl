@muladd begin
#! format: noindent

function geopotential_cartZ(x, equations::CompressibleEulerAtmo{NDIMS}) where
         {NDIMS}
    g = equations.parameters.earth_gravitational_acceleration
    return g * x[NDIMS]
end

# TODO: constant term?
# phi = EARTH_GRAVITATIONAL_ACCELERATION * (EARTH_RADIUS - EARTH_RADIUS^2 / r)
function geopotential_spherical(x,
                                equations::CompressibleEulerAtmo{NDIMS}) where
         {NDIMS}
    g = equations.parameters.earth_gravitational_acceleration
    a = equations.parameters.earth_radius
    r = -norm(x)
    return a^2 * g / r
end
end # @muladd
