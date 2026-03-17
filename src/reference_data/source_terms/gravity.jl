# By convention refers to conservative variables as defined in the equations

@muladd begin
#! format: noindent

@inline function st_gravity_cartZ_td(u, ::CompressibleEulerAtmo{NDIMS}, ::TotalEnergy) where
    {NDIMS}
    return u[NDIMS]
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
        rho = density_total(u, equations)
        u0 = zero(eltype(u))

        u_td = st_gravity_cartZ_td(u, equations, equations.td_equation)

        return SVector(ntuple(i->u0, NDIMS-1)..., -g * rho, -g * u_td)
    end
end


@inline function st_gravity_spherical_td(u, x, equations::CompressibleEulerAtmo{3}, ::TotalEnergy)
    return dot(vars_moment(u, equations), x)
end

function source_terms_gravity_spherical_generator(equations::CompressibleEulerAtmo{3})

    g = equations.parameters.earth_gravitational_acceleration
    a = equations.parameters.radius_earth

    function source_terms(u, x, t, equations)
        r = norm(x)
        rho = u[1]
        u0 = zero(eltype(u))

        temp = -g * a^2 / r^3
        u_td = st_gravity_spherical_td(u, x, equations, equations.td_equation)

        return SVector(u0, (temp * rho_total) .* x, temp * u_td)
    end
end

function source_terms_coriolis_generator(equations::CompressibleEulerAtmo{3})

    omega = equations.parameters.angular_velocity

    function source_terms(u, x, t, equations)
        u0 = zero(eltype(u))

        # Coriolis term, -2Ω × ρv = -2 * angular_velocity * (0, 0, 1) × u[2:4]
        du2 =  2 * omega * u[3]
        du3 = -2 * omega * u[2]

        return SVector(u0, du2, du3, u0, u0)
    end
end
end # @muladd
