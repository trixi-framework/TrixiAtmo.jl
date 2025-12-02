@muladd begin
#! format: noindent

@inline function st_gravity_cartZ_td(u, ::CompressibleEulerAtmo{NDIMS}, ::TotalEnergy) where
    {NDIMS}
    return u[NDIMS+1]
end

@inline function st_gravity_cartZ_td(u, ::CompressibleEulerAtmo{NDIMS},
    ::PotentialTemperature) where
    {NDIMS}
    return zero(eltype(u))
end

function source_terms_gravity_cartZ_generator(equations::CompressibleEulerAtmo{NDIMS}) where
    {NDIMS}

    g = equations.parameters.earth_gravitational_acceleration
    function source_terms(u, x, t, equations)
        rho = u[1]
        u0 = zero(eltype(u))

        u_td = st_gravity_cartZ_td(u, equations, equations.td_equation)

        return SVector(ntuple(i->u0, NDIMS)..., -g * rho, -g * u_td)
    end
end
end # @muladd
