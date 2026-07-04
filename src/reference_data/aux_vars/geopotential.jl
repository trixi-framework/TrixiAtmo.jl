@muladd begin
#! format: noindent

# Steady state for RHS correction below
function geopotential_cartZ(x, t, equations::CompressibleEulerAtmo{NDIMS}) where
         {NDIMS}
    g = equations.parameters.earth_gravitational_acceleration
    return g * x[NDIMS]
end
end # @muladd
