@muladd begin
#! format: noindent

function geopotential_cartZ(x, t, equations::CompressibleEulerAtmo{NDIMS}) where
         {NDIMS}
    g = equations.parameters.earth_gravitational_acceleration
    return g * x[NDIMS]
end

function geopotential_spherical(x, t,
                                equations::CompressibleEulerAtmo{NDIMS}) where
         {NDIMS}
    g = equations.parameters.earth_gravitational_acceleration
    a = equations.parameters.earth_radius
    r = -norm(x)
    return a^2 * g / r
end
end # @muladd
