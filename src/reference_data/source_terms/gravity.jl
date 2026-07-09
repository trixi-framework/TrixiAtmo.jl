# By convention refers to conservative variables as defined in the equations

@muladd begin
#! format: noindent

@inline function st_gravity_cartZ_td(u, ::CompressibleEulerAtmo{NDIMS},
                                     ::EnergyTotal) where
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

        return SVector(ntuple(i -> u0, NDIMS - 1)..., -g * rho, -g * u_td)
    end
    return source_terms
end

@inline function st_gravity_spherical_td(u, x, equations::CompressibleEulerAtmo{3},
                                         ::EnergyTotal)
    return dot(vars_moment(u, equations), x)
end

function source_terms_gravity_spherical_generator(equations::CompressibleEulerAtmo{3,
                                                                                   NVARS,
                                                                                   NGAS,
                                                                                   NCONDENS,
                                                                                   NPRECIP,
                                                                                   NPASSIVE,
                                                                                   0}) where {
                                                                                              NVARS,
                                                                                              NGAS,
                                                                                              NCONDENS,
                                                                                              NPRECIP,
                                                                                              NPASSIVE
                                                                                              }
    g = equations.parameters.earth_gravitational_acceleration
    a = equations.parameters.earth_radius

    function source_terms(u, x, t, equations)
        r = norm(x)
        rho_total = density_total(u, equations)

        temp = -g * a^2 / r^3
        u_td = st_gravity_spherical_td(u, x, equations, equations.td_equation)

        return SVector(((temp * rho_total) * x)..., temp * u_td)
    end
    return source_terms
end

function source_terms_gravity_spherical_generator(equations::CompressibleEulerAtmo{3,
                                                                                   NVARS,
                                                                                   NGAS,
                                                                                   NCONDENS,
                                                                                   NPRECIP,
                                                                                   NPASSIVE,
                                                                                   NAUX}) where {
                                                                                                 NVARS,
                                                                                                 NGAS,
                                                                                                 NCONDENS,
                                                                                                 NPRECIP,
                                                                                                 NPASSIVE,
                                                                                                 NAUX
                                                                                                 }
    g = equations.parameters.earth_gravitational_acceleration
    a = equations.parameters.earth_radius

    function source_terms(u, aux, x, t, equations)
        r = norm(x)
        rho_total = density_total(u, equations)

        temp = -g * a^2 / r^3
        u_td = st_gravity_spherical_td(u, x, equations, equations.td_equation)

        return SVector(((temp * rho_total) * x)..., temp * u_td)
    end
    return source_terms
end

function source_terms_coriolis_generator(equations::CompressibleEulerAtmo{3, NVARS,
                                                                          NGAS,
                                                                          NCONDENS,
                                                                          NPRECIP,
                                                                          NPASSIVE, 0}) where {
                                                                                               NVARS,
                                                                                               NGAS,
                                                                                               NCONDENS,
                                                                                               NPRECIP,
                                                                                               NPASSIVE
                                                                                               }
    omega = equations.parameters.earth_rotation_rate

    function source_terms(u, aux, x, t, equations)
        # Coriolis term, -2Ω × ρv = -2 * angular_velocity * (0, 0, 1) × u[1:3]
        du1 = 2 * omega * u[2]
        du2 = -2 * omega * u[1]

        return SVector(du1, du2)
    end
    return source_terms
end

function source_terms_coriolis_generator(equations::CompressibleEulerAtmo{3, NVARS,
                                                                          NGAS,
                                                                          NCONDENS,
                                                                          NPRECIP,
                                                                          NPASSIVE,
                                                                          NAUX}) where {
                                                                                        NVARS,
                                                                                        NGAS,
                                                                                        NCONDENS,
                                                                                        NPRECIP,
                                                                                        NPASSIVE,
                                                                                        NAUX
                                                                                        }
    omega = equations.parameters.earth_rotation_rate

    function source_terms(u, aux, x, t, equations)
        # Coriolis term, -2Ω × ρv = -2 * angular_velocity * (0, 0, 1) × u[1:3]
        du1 = 2 * omega * u[2]
        du2 = -2 * omega * u[1]

        return SVector(du1, du2)
    end
    return source_terms
end
end # @muladd
