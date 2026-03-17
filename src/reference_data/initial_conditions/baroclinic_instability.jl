@muladd begin
#! format: noindent

#  Functions for setting up idealized initial conditions for the
#  Ullrich, Melvin, Staniforth and Jablonowski baroclinic instability.
#
#  Translated Fortran code orginally written by Paul Ullrich.
#
#  Options:
#    deep_atmosphere                 deep atmosphere
#    moisture                        include moisture
#    perturbation_stream_function    stream function pertrubation, exponential else
#    earth_scaling_factor            Earth scaling factor
#
#  Given a point specified by: 
#      lon    longitude (radians) 
#      lat    latitude (radians) 
#      p/z    pressure (Pa) / height (m)
#  zcoords    1 if z is specified, 0 if p is specified
#
#  the functions will return:
#        p    pressure if z is specified and zcoords = 1 (Pa)
#        u    zonal wind (m s^-1)
#        v    meridional wind (m s^-1)
#        t    temperature (K)
#   thetav    virtual potential temperature (K)
#     phis    surface geopotential (m^2 s^-2)
#       ps    surface pressure (Pa)
#      rho    density (kj m^-3)
#        q    water vapor mixing ratio (kg/kg)
#
function initial_condition_baroclinic_instability_generator(
        parameters::Parameters{RealType};
        deep_atmosphere              = True(),
        moisture                     = False(),
        perturbation_stream_function = False(),
        earth_scaling_factor         = RealType(1),
    ) where RealType

    # Test case parameters (currently fixed)
    T0E        = RealType(310)       # temperature at equatorial surface (K)
    T0P        = RealType(240)       # temperature at polar surface (K)
    B          = RealType(2)         # jet half-width parameter
    K          = RealType(3)         # jet width parameter
    lapse      = RealType(0.005)     # lapse rate parameter
  
    pertu0     = RealType(0.5)       # SF Perturbation wind velocity (m/s)
    pertr      = RealType(1/6)       # SF Perturbation radius (Earth radii)
    pertup     = RealType(1.0)       # Exp. perturbation wind velocity (m/s)
    pertexpr   = RealType(0.1)       # Exp. perturbation radius (Earth radii)
    pertlon    = RealType(pi/9)      # Perturbation longitude
    pertlat    = RealType(2*pi/9)    # Perturbation latitude
    pertz      = RealType(15000)     # Perturbation height cap

    moistqlat  = RealType(2*pi/9)    # Humidity latitudinal width
    moistqp    = RealType(34000)     # Humidity vertical pressure width
    moisttr    = RealType(0.1)       # Vertical cut-off pressure for humidity
    moistqs    = RealType(1e-12)     # Humidity above cut-off
    moistq0    = RealType(0.018)     # Maximum specific humidity
    Mvap       = RealType(0.608)     # Ratio of molar mass of dry air/water

    # Physical parameters (taken from TrixiAtmo parameters)
    a          = parameters.earth_radius
    g          = parameters.earth_gravitational_acceleration
    cp         = parameters.c_dry_air_const_pressure
    p0         = parameters.ref_pressure
    omega      = parameters.earth_rotation_rate
    dxepsilon  = parameters.tol_eps  # Small value for numerical derivatives
  
    # local constants
    X          = earth_scaling_factor
    Rd         = cp - cv             # Ideal gas const dry air (J kg^-1 K^1)
    aref       = a / X
    omegaref   = omega * X
    T0         = 0.5f0 * (T0E + T0P)
    constA     = 1 / lapse
    constB     = (T0 - T0P) / (T0 * T0P)
    constC     = 0.5f0 * (K + 2) * (T0E - T0P) / (T0E * T0P)
    constH     = Rd * T0 / g
    
    function initial_condtion(cart_coords, t, equations)
        spherical_coords = cartesian_to_spherical_coordinates(cart_coords)
        spherical_u, spherical_v, T, thetav, rho, p, q = baroclinic_wave_test(spherical_coords)
        spherical_velocities = SVector(spherical_u, spherical_v, 0)
        cart_u, cart_v, cart_w = spherical_to_cartesian_velocities(spherical_coords,
                                                                   spherical_velocities)
        return SVector(rho, cart_u, cart_v, cart_w, p)
    end

    function baroclinic_wave_test(spherical_coords)
        lon, lat, z = spherical_coords

        p, T = evaluate_pressure_temperature(lon, lat, z)

        scaledZ = z / (B * constH)
        inttau2 = constC * z * exp(- scaledZ^2)

        # Radius ratio
        if deep_atmosphere
            rratio = (z + aref) / aref;
        else
            rratio = 1
        end

        # Initialize velocity field
        inttermU = (rratio * cos(lat))^(K - 1) - (rratio * cos(lat))^(K + 1)
        bigU = g / aref * K * inttau2 * inttermU * T

        if !deep_atmosphere then
            rcoslat = aref * cos(lat)
        else
            rcoslat = (z + aref) * cos(lat)
        end

        omegarcoslat = omegaref * rcoslat
        u = - omegarcoslat + sqrt(omegarcoslat^2 + rcoslat * bigU)
        v = 0

        # Add perturbation to the velocity field
        if perturbation_stream_function
            u = u - 1 / (2 * dxepsilon) *
                ( evaluate_streamfunction(lon, lat + dxepsilon, z)
                - evaluate_streamfunction(lon, lat - dxepsilon, z))
            v = v + 1 / (2 * dxepsilon * cos(lat)) *
                ( evaluate_streamfunction(lon + dxepsilon, lat, z)
                - evaluate_streamfunction(lon - dxepsilon, lat, z))
        else
            u = u + evaluate_exponential(lon, lat, z)
        end

        # Initialize density
        rho = p / (Rd * T)

        # Initialize specific humidity
        if moisture
            eta = p/p0
            if eta > moisttr
                q = moistq0 * exp(- (lat/moistqlat)^4) * exp(- ((eta-1)*p0/moistqp)^2)
            else
                q = moistqs
            end
            # Convert virtual temperature to temperature
            T = T / (1 + Mvap * q)
        else
            q = 0
        end

        # Initialize virtual potential temperature
        thetav = T * (1 + RealType(0.61) * q) * (p0 / p)^(Rd / cp)

        # TODO: need to re-calculate p?

        return u, v, T, thetav, rho, p, q
    end

    function evaluate_pressure_temperature(lon, lat, z)

        scaledZ = z / (B * constH)
        tau1 = constA * lapse / T0 * exp(lapse * z / T0) +
               constB * (1 - 2 * scaledZ^2) * exp(- scaledZ^2)
        tau2 = constC * (1 - 2 * scaledZ^2) * exp(- scaledZ^2)

        inttau1 = constA * (exp(lapse * z / T0) - 1) +
                  constB * z * exp(- scaledZ^2)
        inttau2 = constC * z * exp(- scaledZ^2)

        if !deep_atmosphere
            rratio = 1
        else
            rratio = (z + aref) / aref
        end

        # interior term on temperature expression
        inttermT = (rratio * cos(lat))^K - K / (K + 2) * (rratio * cos(lat))^(K + 2)

        # temperature
        T = 1 / (rratio^2 * (tau1 - tau2 * inttermT))

        # hydrostatic pressure
        p = p0 * exp(- g / Rd * (inttau1 - inttau2 * inttermT))

        return p, T
    end

    function evaluate_exponential(lon, lat, z)
    
        # Great circle distance
        greatcircler = 1 / pertexpr * acos(sin(pertlat) * sin(lat) +
                       cos(pertlat) * cos(lat) * cos(lon - pertlon))

        # Vertical tapering of stream function
        if z < pertz
            perttaper = 1 - 3 * z^2 / pertz^2 + 2 * z^3 / pertz^3
        else
            perttaper = 0
        end

        # Zonal velocity perturbation
        if greatcircler < 1
            exponential = pertup * perttaper * exp(- greatcircler^2)
        else
            exponential = 0
        end
        return exponential
    end
  
    function evaluate_streamfunction(lon, lat, z)

        # Great circle distance
        greatcircler = 1 / pertr * acos(sin(pertlat) * sin(lat) +
                       cos(pertlat) * cos(lat) * cos(lon - pertlon))

        # Vertical tapering of stream function
        if z < pertz
            perttaper = 1 - 3 * z^2 / pertz^2 + 2 * z^3 / pertz^3
        else
            perttaper = 0
        endif

        # Horizontal tapering of stream function
        if greatcircler < 1
            cospert = cos(0.5f0 * pi * greatcircler)
        else
            cospert = 0
        end

        return - pertu0 * pertr * perttaper * cospert^4
    end

    return initial_condtion
end
end
end # @muladd
