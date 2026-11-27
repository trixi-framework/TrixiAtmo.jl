@testitem "Unit: Consistency check for EC flux with Potential Temperature: CEPTE" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquations2D(c_p = 1004.0,
                                                                 c_v = 717.0,
                                                                 gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 330.0)
    aux = SVector(1500.0)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        @test flux_ec(u, u, aux, aux, normal_direction, equations) ≈
              flux(u, aux, normal_direction, equations)
    end

    #check consistency between 1D and 2D EC fluxes
    u_1d = SVector(u[1], u[2], u[4], aux[1])
    u_2d = SVector(u[1], u[2], 0, u[4])
    normal_1d = SVector(-0.3)
    normal_2d = SVector(normal_1d[1], 0.0)
    equations_1d = CompressibleEulerPotentialTemperatureEquations1D(c_p = 1004.0,
                                                                    c_v = 717.0,
                                                                    gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    equations_2d = equations
    flux_1d = normal_1d[1] * flux_ec(u_1d, u_1d, 1, equations_1d)
    flux_2d = flux_ec(u_2d, u_2d, aux, aux, normal_2d, equations_2d)
    @test flux_1d ≈ flux(u_1d, normal_1d, equations_1d)
    @test flux_1d[1:3] ≈ flux_2d[[1, 2, 4]]

    # test when u_ll is not the same as u_rr
    u_rr_1d = SVector(2.1, 0.3, 280.5, 1700)
    u_rr_2d = SVector(u_rr_1d[1], u_rr_1d[2], 0.0, u_rr_1d[3])
    aux_rr = SVector(u_rr_1d[4])
    flux_1d = normal_1d[1] * flux_ec(u_1d, u_rr_1d, 1, equations_1d)
    flux_2d = flux_ec(u_2d, u_rr_2d, aux, aux_rr, normal_2d, equations_2d)
    @test flux_1d[1:3] ≈ flux_2d[[1, 2, 4]]

    # check consistency for 3D EC flux
    equations = CompressibleEulerPotentialTemperatureEquations3D(c_p = 1004.0,
                                                                 c_v = 717.0,
                                                                 gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 2.4, 330.0)
    aux = SVector(1500.0)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        @test flux_ec(u, u, aux, aux, normal_direction, equations) ≈
              flux(u, aux, normal_direction, equations)
    end

    # check consistency between 1D and 3D EC fluxes
    u_3d = SVector(u[1], u[2], 0.0, 0.0, u[5])
    normal_3d = SVector(normal_1d[1], 0.0, 0.0)
    equations_3d = equations
    flux_3d = flux_ec(u_3d, u_3d, aux, aux, normal_3d, equations_3d)
    flux_1d = normal_1d[1] * flux_ec(u_1d, u_1d, 1, equations_1d)
    @test flux_1d[1:3] ≈ flux_3d[[1, 2, 5]]
end

@testitem "Unit: Consistency check for FluxTurbo EC flux with Potential Temperature: CEPTE" setup=[Setup] tags=[:unit_fluxes] begin
    # The `FluxTurbo` specialization of `(flux_ec, flux_nonconservative_isothermal)`
    # must return the same left and right contributions as the fused evaluation of the
    # conservative and the nonconservative flux.
    equations = CompressibleEulerPotentialTemperatureEquations3D(c_p = 1004.0,
                                                                 c_v = 717.0,
                                                                 gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    volume_flux = (flux_ec, flux_nonconservative_isothermal)
    u_ll = SVector(1.1, -0.5, 2.34, 2.4, 330.0)
    u_rr = SVector(0.9, 0.3, -1.2, 1.7, 280.5)
    aux_ll = SVector(1500.0)
    aux_rr = SVector(1700.0)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    @test Trixi.nturbovars(volume_flux, equations) == Val(7)
    turbo_ll = Trixi.cons2turbo(volume_flux, u_ll..., aux_ll..., equations)
    turbo_rr = Trixi.cons2turbo(volume_flux, u_rr..., aux_rr..., equations)
    @test length(turbo_ll) == 7
    # rho, v1, v2, v3, theta, p, phi
    @test SVector(turbo_ll[1], turbo_ll[2], turbo_ll[3], turbo_ll[4], turbo_ll[6]) ≈
          cons2prim(u_ll, aux_ll, equations)
    @test turbo_ll[5] ≈ u_ll[5] / u_ll[1]
    @test turbo_ll[7] == aux_ll[1]

    for normal_direction in normal_directions
        flux_left, flux_right = Trixi.flux_turbo(volume_flux, turbo_ll..., turbo_rr...,
                                                 normal_direction..., equations)
        flux_cons = flux_ec(u_ll, u_rr, aux_ll, aux_rr, normal_direction, equations)
        flux_noncons = flux_nonconservative_isothermal(u_ll, u_rr, aux_ll, aux_rr,
                                                       normal_direction, equations)
        @test SVector(flux_left) ≈ flux_cons + 0.5 * flux_noncons
        @test SVector(flux_right) ≈ flux_cons - 0.5 * flux_noncons
    end
end

@testitem "Unit: Consistency check for TEC flux with Potential Temperature: CEPTE" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquations2D(c_p = 1004.0,
                                                                 c_v = 717.0,
                                                                 gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 330.0)
    aux = SVector(1500.0)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        @test flux_tec(u, u, aux, aux, normal_direction, equations) ≈
              flux(u, aux, normal_direction, equations)
    end

    #check consistency between 1D and 2D EC fluxes
    u_1d = SVector(u[1], u[2], u[4], aux[1])
    u_2d = SVector(u[1], u[2], 0, u[4])
    normal_1d = SVector(-0.3)
    normal_2d = SVector(normal_1d[1], 0.0)
    equations_1d = CompressibleEulerPotentialTemperatureEquations1D(c_p = 1004.0,
                                                                    c_v = 717.0,
                                                                    gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    equations_2d = equations
    flux_1d = normal_1d[1] * flux_tec(u_1d, u_1d, 1, equations_1d)
    flux_2d = flux_tec(u_2d, u_2d, aux, aux, normal_2d, equations_2d)
    @test flux_1d ≈ flux(u_1d, normal_1d, equations_1d)
    @test flux_1d[1:3] ≈ flux_2d[[1, 2, 4]]

    # test when u_ll is not the same as u_rr
    u_rr_1d = SVector(2.1, 0.3, 280.5, 1700)
    u_rr_2d = SVector(u_rr_1d[1], u_rr_1d[2], 0.0, u_rr_1d[3])
    aux_rr = SVector(u_rr_1d[4])
    flux_1d = normal_1d[1] * flux_tec(u_1d, u_rr_1d, 1, equations_1d)
    flux_2d = flux_tec(u_2d, u_rr_2d, aux, aux_rr, normal_2d, equations_2d)
    @test flux_1d[1:3] ≈ flux_2d[[1, 2, 4]]

    # check consistency for 3D EC flux
    equations = CompressibleEulerPotentialTemperatureEquations3D(c_p = 1004.0,
                                                                 c_v = 717.0,
                                                                 gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 2.4, 330.0)
    aux = SVector(1500.0)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        @test flux_tec(u, u, aux, aux, normal_direction, equations) ≈
              flux(u, aux, normal_direction, equations)
    end

    # check consistency between 1D and 3D EC fluxes
    u_3d = SVector(u[1], u[2], 0.0, 0.0, u[5])
    normal_3d = SVector(normal_1d[1], 0.0, 0.0)
    equations_3d = equations
    flux_3d = flux_tec(u_3d, u_3d, aux, aux, normal_3d, equations_3d)
    flux_1d = normal_1d[1] * flux_tec(u_1d, u_1d, 1, equations_1d)
    @test flux_1d[1:3] ≈ flux_3d[[1, 2, 5]]
end

@testitem "Unit: Consistency check for ETEC flux with Potential Temperature: CEPTE" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquations2D(c_p = 1004.0,
                                                                 c_v = 717.0,
                                                                 gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 330.0)
    aux = SVector(1500.0)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        @test flux_etec(u, u, aux, aux, normal_direction, equations) ≈
              flux(u, aux, normal_direction, equations)
    end

    #check consistency between 1D and 2D EC fluxes
    u_1d = SVector(u[1], u[2], u[4], aux[1])
    u_2d = SVector(u[1], u[2], 0, u[4])
    normal_1d = SVector(-0.3)
    normal_2d = SVector(normal_1d[1], 0.0)
    equations_1d = CompressibleEulerPotentialTemperatureEquations1D(c_p = 1004.0,
                                                                    c_v = 717.0,
                                                                    gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    equations_2d = equations
    flux_1d = normal_1d[1] * flux_etec(u_1d, u_1d, 1, equations_1d)
    flux_2d = flux_etec(u_2d, u_2d, aux, aux, normal_2d, equations_2d)
    @test flux_1d ≈ flux(u_1d, normal_1d, equations_1d)
    @test flux_1d[1:3] ≈ flux_2d[[1, 2, 4]]

    # test when u_ll is not the same as u_rr
    u_rr_1d = SVector(2.1, 0.3, 280.5, 1700)
    u_rr_2d = SVector(u_rr_1d[1], u_rr_1d[2], 0.0, u_rr_1d[3])
    aux_rr = SVector(u_rr_1d[4])
    flux_1d = normal_1d[1] * flux_etec(u_1d, u_rr_1d, 1, equations_1d)
    flux_2d = flux_etec(u_2d, u_rr_2d, aux, aux_rr, normal_2d, equations_2d)
    @test flux_1d[1:3] ≈ flux_2d[[1, 2, 4]]

    # check consistency for 3D EC flux
    equations = CompressibleEulerPotentialTemperatureEquations3D(c_p = 1004.0,
                                                                 c_v = 717.0,
                                                                 gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 2.4, 330.0)
    aux = SVector(1500.0)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        @test flux_etec(u, u, aux, aux, normal_direction, equations) ≈
              flux(u, aux, normal_direction, equations)
    end

    # check consistency between 1D and 3D EC fluxes
    u_3d = SVector(u[1], u[2], 0.0, 0.0, u[5])
    normal_3d = SVector(normal_1d[1], 0.0, 0.0)
    equations_3d = equations
    flux_3d = flux_etec(u_3d, u_3d, aux, aux, normal_3d, equations_3d)
    flux_1d = normal_1d[1] * flux_etec(u_1d, u_1d, 1, equations_1d)
    @test flux_1d[1:3] ≈ flux_3d[[1, 2, 5]]
end

@testitem "Unit: Consistency check for LMARS flux with Potential Temperature: CEPTE" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquations2D(c_p = 1004.0,
                                                                 c_v = 717.0,
                                                                 gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    flux_lmars = FluxLMARS(340)
    u = SVector(1.1, -0.5, 2.34, 330.0)
    aux = SVector(1700.0)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        @test flux_lmars(u, u, aux, aux, normal_direction, equations) ≈
              flux(u, aux, normal_direction, equations)
    end

    #check consistency between 1D and 2D EC fluxes
    u_1d = SVector(u[1], u[2], u[4], aux[1])
    u_2d = SVector(u[1], u[2], 0, u[4])
    normal_1d = SVector(-0.3)
    normal_2d = SVector(normal_1d[1], 0.0)
    equations_1d = CompressibleEulerPotentialTemperatureEquations1D(c_p = 1004.0,
                                                                    c_v = 717.0,
                                                                    gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    equations_2d = equations
    flux_1d = normal_1d[1] * flux_lmars(u_1d, u_1d, 1, equations_1d)
    flux_2d = flux_lmars(u_2d, u_2d, aux, aux, normal_2d, equations_2d)
    @test flux_1d ≈ flux(u_1d, normal_1d, equations_1d)
    @test flux_1d[1:3] ≈ flux_2d[[1, 2, 4]]

    # check consistency for 3D EC flux
    equations = CompressibleEulerPotentialTemperatureEquations3D(c_p = 1004.0,
                                                                 c_v = 717.0,
                                                                 gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 2.4, 330.0)
    aux = SVector(1700.0)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        @test flux_lmars(u, u, aux, aux, normal_direction, equations) ≈
              flux(u, aux, normal_direction, equations)
    end

    # check consistency between 1D and 3D EC fluxes
    u_3d = SVector(u[1], u[2], 0.0, 0.0, u[5])
    normal_3d = SVector(normal_1d[1], 0.0, 0.0)
    equations_3d = equations
    flux_3d = flux_lmars(u_3d, u_3d, aux, aux, normal_3d, equations_3d)
    flux_1d = normal_1d[1] * flux_lmars(u_1d, u_1d, 1, equations_1d)
    @test flux_1d[1:3] ≈ flux_3d[[1, 2, 5]]
end

@testitem "Unit: Consistency check for ES LMARS flux with Potential Temperature: CEPTE" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquations2D(c_p = 1004.0,
                                                                 c_v = 717.0,
                                                                 gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 330.0)
    u_rr = SVector(0.9, 0.3, -1.2, 280.5)
    aux = SVector(1700.0)
    aux_rr = SVector(1500.0)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        # Consistency with the physical flux
        for i in 1:2
            @test flux_lmars_es(u, u, aux, aux, normal_direction, equations)[i] ≈
                  flux(u, aux, normal_direction, equations)
        end

        flux_left, flux_right = flux_lmars_es(u, u_rr, aux, aux_rr, normal_direction,
                                              equations)
        @test flux_left == flux_right
        @test flux_lmars_es(u, u_rr, aux, aux_rr, 2.5 * normal_direction,
                            equations)[1] ≈ 2.5 * flux_left
    end

    # Check consistency for the 3D flux
    equations_3d = CompressibleEulerPotentialTemperatureEquations3D(c_p = 1004.0,
                                                                    c_v = 717.0,
                                                                    gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u_3d = SVector(1.1, -0.5, 2.34, 2.4, 330.0)
    u_rr_3d = SVector(0.9, 0.3, -1.2, 1.7, 280.5)
    aux_3d = SVector(1700.0)
    aux_rr_3d = SVector(1500.0)

    normal_directions_3d = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions_3d
        for i in 1:2
            @test flux_lmars_es(u_3d, u_3d, aux_3d, aux_3d, normal_direction,
                                equations_3d)[i] ≈
                  flux(u_3d, aux_3d, normal_direction, equations_3d)
        end
        flux_left, flux_right = flux_lmars_es(u_3d, u_rr_3d, aux_3d, aux_rr_3d,
                                              normal_direction, equations_3d)
        @test flux_left == flux_right
        @test flux_lmars_es(u_3d, u_rr_3d, aux_3d, aux_rr_3d, 2.5 * normal_direction,
                            equations_3d)[1] ≈ 2.5 * flux_left
    end

    # Check consistency between the 2D and 3D fluxes for different left and right states
    u_3d = SVector(u[1], u[2], u[3], 0.0, u[4])
    u_rr_3d = SVector(u_rr[1], u_rr[2], u_rr[3], 0.0, u_rr[4])
    for normal_direction in normal_directions
        normal_3d = SVector(normal_direction[1], normal_direction[2], 0.0)
        flux_2d = flux_lmars_es(u, u_rr, aux, aux_rr, normal_direction, equations)
        flux_3d = flux_lmars_es(u_3d, u_rr_3d, aux, aux_rr, normal_3d, equations_3d)
        for i in 1:2
            @test flux_2d[i] ≈ flux_3d[i][[1, 2, 3, 5]]
        end
    end
end

@testitem "Unit: Consistency check for gravity terms with energy: CEEE" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables states with the same temperature
    equations = CompressibleEulerEnergyEquations3D(c_p = 1004.0, c_v = 717.0,
                                                   gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    aux_ll = SVector(1700.0)
    aux_rr = SVector(2900.0)
    u_ll = prim2cons(SVector(1.1, 2.0, -1.5, 0.5, 9.0e4), aux_ll, equations)
    u_rr = prim2cons(SVector(0.9, -1.0, 3.0, 1.2, 0.9 / 1.1 * 9.0e4), aux_rr, equations)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    # Both gravity terms coincide for an isothermal state
    for normal_direction in normal_directions
        @test flux_nonconservative_waruszewski_etal(u_ll, u_rr, aux_ll, aux_rr,
                                                    normal_direction, equations) ≈
              flux_nonconservative_isothermal(u_ll, u_rr, aux_ll, aux_rr,
                                              normal_direction, equations)
    end

    # Check consistency between the 2D and 3D gravity terms
    equations_2d = CompressibleEulerEnergyEquations2D(c_p = 1004.0, c_v = 717.0,
                                                      gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    aux_ll_2d = SVector(1700.0)
    aux_rr_2d = SVector(2900.0)
    u_ll_2d = prim2cons(SVector(1.1, 2.0, -1.5, 9.0e4), aux_ll_2d, equations_2d)
    u_rr_2d = prim2cons(SVector(0.9, -1.0, 3.0, 7.5e4), aux_rr_2d, equations_2d)
    u_ll_3d = prim2cons(SVector(1.1, 2.0, -1.5, 0.0, 9.0e4), aux_ll, equations)
    u_rr_3d = prim2cons(SVector(0.9, -1.0, 3.0, 0.0, 7.5e4), aux_rr, equations)
    for normal_direction in (SVector(1.0, 0.0), SVector(0.0, 1.0), SVector(0.5, -0.5),
                             SVector(-1.2, 0.3))
        normal_direction_3d = SVector(normal_direction[1], normal_direction[2], 0.0)
        for flux_gravity in (flux_nonconservative_waruszewski_etal,
                             flux_nonconservative_isothermal)
            @test flux_gravity(u_ll_2d, u_rr_2d, aux_ll_2d, aux_rr_2d, normal_direction,
                               equations_2d) ≈
                  flux_gravity(u_ll_3d, u_rr_3d, aux_ll, aux_rr, normal_direction_3d,
                               equations)[[1, 2, 3, 5]]
        end
    end
end

@testitem "Unit: Consistency check for 3D shallow water fluxes: SWE" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = ShallowWaterEquations3D(gravity = 1.0)
    u = SVector(1.1, -0.5, 2.34, -3.5, 120.0)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        @test flux_wintermeyer_etal(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    for normal_direction in normal_directions
        @test flux_fjordholm_etal(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end
end

@testitem "Unit: Consistency check for split covariant shallow water fluxes: SWE" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = SplitCovariantShallowWaterEquations2D(EARTH_GRAVITATIONAL_ACCELERATION,
                                                      EARTH_ROTATION_RATE)
    u = SVector(1.1, -0.5, 2.34)
    aux_vars = SVector{26}(ones(26))
    orientation = 1
    @test flux_ec(u, u, aux_vars, aux_vars, orientation, equations) ≈
          flux(u, aux_vars, orientation, equations)
end

@testitem "Unit: Consistency check for EC flux with Rainy Euler" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleRainyEulerEquations2D()
    # Example state vector (ρ_d, ρ_m, ρ_r, ρu, ρv, ρe, ρq_v, ρq_c, T)
    u = SVector(1.0, 0.2, 0.1, 0.5, -0.4, 2.2, 0.1, 0.1, 300)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        equal = flux_ec_rain(u, u, normal_direction, equations) .≈
                flux(u, normal_direction, equations)
        # TODO
        expected = [true, true, true, true, true, false, true, true, true]
        @test equal == expected
    end
end

@testitem "Unit: check_axes for 2D manifolds in 3D" setup=[Setup] tags=[
    :unit_check_axes,
    :upstream
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "shallow_water", "cartesian",
                                 "elixir_unsteady_solid_body_rotation_ec_correction.jl"),
                        cells_per_dimension=(3, 3), maxiters=1)

    @testset "Cartesian form" begin
        mesh, equations, dg, cache = Trixi.mesh_equations_solver_cache(semi)
        u = Trixi.wrap_array(Trixi.compute_coefficients(0.0, semi), semi)

        @test Trixi.check_axes(u, mesh, equations, dg, cache) === nothing

        u_too_few = similar(u, size(u)[1:(end - 1)]..., size(u, ndims(u)) - 1)
        u_too_many = similar(u, size(u)[1:(end - 1)]..., size(u, ndims(u)) + 1)
        @test_throws DimensionMismatch Trixi.check_axes(u_too_few, mesh, equations, dg,
                                                        cache)
        @test_throws DimensionMismatch Trixi.check_axes(u_too_many, mesh, equations, dg,
                                                        cache)

        @test Trixi.ninterfaces(dg, cache) > 0
        for container in (cache.elements, cache.interfaces, cache.boundaries,
                          cache.mortars)
            @test Trixi.check_axes(container, equations, dg, cache) === nothing
        end

        # The element container of TrixiAtmo.jl is really checked, i.e., it does not use
        # the no-op fallback of Trixi.jl
        dg_wrong = DGSEM(polydeg = Trixi.polydeg(dg) + 1)
        @test_throws DimensionMismatch Trixi.check_axes(cache.elements, equations,
                                                        dg_wrong, cache)

        # Explicit bounds checks of TrixiAtmo.jl before assuming inbounds access
        @test_throws DimensionMismatch TrixiAtmo.calc_sources_2d_manifold_in_3d!(u_too_few,
                                                                                 u, 0.0,
                                                                                 semi.source_terms,
                                                                                 equations,
                                                                                 dg, cache)
    end

    @test_trixi_include(joinpath(EXAMPLES_DIR, "shallow_water", "covariant",
                                 "elixir_unsteady_solid_body_rotation_ec.jl"),
                        cells_per_dimension=(3, 3), maxiters=1)

    @testset "Covariant form" begin
        mesh, equations, dg, cache = Trixi.mesh_equations_solver_cache(semi)
        u = Trixi.wrap_array(Trixi.compute_coefficients(0.0, semi), semi)

        @test Trixi.check_axes(u, mesh, equations, dg, cache) === nothing

        u_too_few = similar(u, size(u)[1:(end - 1)]..., size(u, ndims(u)) - 1)
        u_too_many = similar(u, size(u)[1:(end - 1)]..., size(u, ndims(u)) + 1)
        @test_throws DimensionMismatch Trixi.check_axes(u_too_few, mesh, equations, dg,
                                                        cache)
        @test_throws DimensionMismatch Trixi.check_axes(u_too_many, mesh, equations, dg,
                                                        cache)

        @test Trixi.ninterfaces(dg, cache) > 0
        for container in (cache.elements, cache.interfaces, cache.boundaries,
                          cache.auxiliary_variables)
            @test Trixi.check_axes(container, equations, dg, cache) === nothing
        end

        # The containers of TrixiAtmo.jl are really checked, i.e., the element container
        # does not use the no-op fallback of Trixi.jl
        dg_wrong = DGSEM(polydeg = Trixi.polydeg(dg) + 1)
        @test_throws DimensionMismatch Trixi.check_axes(cache.elements, equations,
                                                        dg_wrong, cache)
        @test_throws DimensionMismatch Trixi.check_axes(cache.auxiliary_variables,
                                                        equations, dg_wrong, cache)

        # Explicit bounds checks of TrixiAtmo.jl before assuming inbounds access
        @test_throws DimensionMismatch Trixi.apply_jacobian!(nothing, u_too_few, mesh,
                                                             equations, dg, cache)
        @test_throws DimensionMismatch Trixi.calc_sources!(u_too_few, u, 0.0,
                                                           semi.source_terms, equations,
                                                           dg, cache)
    end
end

@testitem "Unit: Consistency check for ETEC flux with vector-invariant equations: CEVIEWG" setup=[Setup] tags=[:unit_fluxes] begin
    equations = CompressibleEulerVectorInvariantEquations3D(c_p = 1004.0,
                                                            c_v = 717.0,
                                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u_ll = SVector(1.1, -0.5, 2.34, 2.4, 330.0)
    u_rr = SVector(0.9, 0.3, -1.2, 1.7, 280.5)
    aux_ll = SVector(1500.0)
    aux_rr = SVector(1700.0)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        # Consistency with the conservative part of the physical flux
        rho, v1, v2, v3, rho_theta = u_ll
        v_normal = v1 * normal_direction[1] + v2 * normal_direction[2] +
                   v3 * normal_direction[3]
        kinetic_energy = 0.5 * (v1^2 + v2^2 + v3^2)
        flux_physical = SVector(rho * v_normal, kinetic_energy * normal_direction[1],
                                kinetic_energy * normal_direction[2],
                                kinetic_energy * normal_direction[3],
                                rho_theta * v_normal)
        for i in 1:2
            @test flux_etec(u_ll, u_ll, aux_ll, aux_ll, normal_direction, equations)[i] ≈
                  flux_physical
        end

        # The nonconservative contribution is antisymmetric
        flux_left, flux_right = flux_etec(u_ll, u_rr, aux_ll, aux_rr, normal_direction,
                                          equations)
        flux_left_swapped, flux_right_swapped = flux_etec(u_rr, u_ll, aux_rr, aux_ll,
                                                          normal_direction, equations)
        @test flux_left ≈ flux_right_swapped
        @test flux_right ≈ flux_left_swapped

        # The `FluxTurbo` specialization returns the same contributions
        turbo_ll = Trixi.cons2turbo(flux_etec, u_ll..., aux_ll..., equations)
        turbo_rr = Trixi.cons2turbo(flux_etec, u_rr..., aux_rr..., equations)
        flux_left_turbo, flux_right_turbo = Trixi.flux_turbo(flux_etec, turbo_ll...,
                                                             turbo_rr...,
                                                             normal_direction...,
                                                             equations)
        @test SVector(flux_left_turbo) ≈ flux_left
        @test SVector(flux_right_turbo) ≈ flux_right
    end
    @test Trixi.nturbovars(flux_etec, equations) == Val(8)

    # The entropy variables are the derivatives of the entropy
    epsilon = 1.0e-6
    for v in (1, 5)
        e_v = SVector(ntuple(i -> i == v ? epsilon : 0.0, 5))
        @test isapprox(cons2entropy(u_ll, aux_ll, equations)[v],
                       (entropy(u_ll + e_v, aux_ll, equations) -
                        entropy(u_ll - e_v, aux_ll, equations)) / (2 * epsilon);
                       rtol = 1.0e-6)
    end
end

@testitem "Unit: Consistency check for ES LMARS flux with vector-invariant equations: CEVIEWG" setup=[Setup] tags=[:unit_fluxes] begin
    equations = CompressibleEulerVectorInvariantEquations3D(c_p = 1004.0,
                                                            c_v = 717.0,
                                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u_ll = SVector(1.1, -0.5, 2.34, 2.4, 330.0)
    u_rr = SVector(0.9, 0.3, -1.2, 1.7, 280.5)
    aux_ll = SVector(1500.0)
    aux_rr = SVector(1500.0)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        # Consistency with the conservative part of the physical flux
        rho, v1, v2, v3, rho_theta = u_ll
        v_normal = v1 * normal_direction[1] + v2 * normal_direction[2] +
                   v3 * normal_direction[3]
        kinetic_energy = 0.5 * (v1^2 + v2^2 + v3^2)
        flux_physical = SVector(rho * v_normal, kinetic_energy * normal_direction[1],
                                kinetic_energy * normal_direction[2],
                                kinetic_energy * normal_direction[3],
                                rho_theta * v_normal)
        for i in 1:2
            @test flux_lmars_es(u_ll, u_ll, aux_ll, aux_ll, normal_direction,
                                equations)[i] ≈ flux_physical
        end

        # The flux is homogeneous of degree one in the normal direction
        for i in 1:2
            @test flux_lmars_es(u_ll, u_rr, aux_ll, aux_rr, 2.5 * normal_direction,
                                equations)[i] ≈
                  2.5 * flux_lmars_es(u_ll, u_rr, aux_ll, aux_rr, normal_direction,
                                equations)[i]
        end
    end
end

@testitem "Unit: Conversions with auxiliary variables for 3D equations" setup=[Setup] tags=[:unit_fluxes] begin
    # `prim2cons` is the inverse of `cons2prim`
    prim = SVector(1.1, 2.0, -1.5, 0.5, 9.0e4)
    aux = SVector(1700.0)
    for equations in (CompressibleEulerPotentialTemperatureEquations3D(c_p = 1004.0,
                                                                       c_v = 717.0,
                                                                       gravity = EARTH_GRAVITATIONAL_ACCELERATION),
                      CompressibleEulerEnergyEquations3D(c_p = 1004.0,
                                                         c_v = 717.0,
                                                         gravity = EARTH_GRAVITATIONAL_ACCELERATION),
                      CompressibleEulerVectorInvariantEquations3D(c_p = 1004.0,
                                                                  c_v = 717.0,
                                                                  gravity = EARTH_GRAVITATIONAL_ACCELERATION))
        u = prim2cons(prim, aux, equations)
        @test cons2prim(u, aux, equations) ≈ prim
    end
end

@testitem "Unit: Consistency check for fluxes with Moist Potential Temperature: CEMPTE3D" setup=[Setup] tags=[:unit_fluxes] begin
    equations = CompressibleEulerMoistPotentialTemperatureEquations3D(c_p = 1004.0,
                                                                      c_v = 717.0,
                                                                      gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u_ll = SVector(1.1, -0.5, 2.34, 0.7, 310.0, 0.012, 0.003, 0.002)
    u_rr = SVector(1.05, 0.4, -1.3, -0.2, 300.0, 0.01, 0.0, 0.004)
    # Geopotential and radial upward unit vector at two nearby points on a sphere
    function aux_sphere(x)
        r = sqrt(sum(abs2, x))
        return SVector(EARTH_GRAVITATIONAL_ACCELERATION * (r - 3.0e5), x[1] / r, x[2] / r,
                       x[3] / r)
    end
    x_ll = SVector(2.0e5, -1.0e5, 2.0e5)
    x_rr = SVector(2.1e5, -0.9e5, 2.05e5)
    aux_ll = aux_sphere(x_ll)
    aux_rr = aux_sphere(x_rr)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions, (u, aux) in ((u_ll, aux_ll), (u_rr, aux_rr))
        for i in 1:2
            @test flux_lmars_combined(u, u, aux, aux, normal_direction, equations)[i] ≈
                  flux(u, aux, normal_direction, equations)
            @test flux_kennedy_gruber_souza_etal(u, u, aux, aux, normal_direction,
                                                 equations)[i] ≈
                  flux(u, aux, normal_direction, equations)
        end
    end

    for normal_direction in normal_directions
        # The conservative part of the volume flux is symmetric and the nonconservative
        # part is antisymmetric
        flux_left, flux_right = flux_kennedy_gruber_souza_etal(u_ll, u_rr, aux_ll, aux_rr,
                                                               normal_direction,
                                                               equations)
        flux_left_swapped, flux_right_swapped = flux_kennedy_gruber_souza_etal(u_rr, u_ll,
                                                                               aux_rr,
                                                                               aux_ll,
                                                                               normal_direction,
                                                                               equations)
        @test flux_left + flux_right ≈ flux_left_swapped + flux_right_swapped
        @test flux_left - flux_right ≈ -(flux_left_swapped - flux_right_swapped)

        # The `FluxTurbo` specialization returns the same contributions
        turbo_ll = Trixi.cons2turbo(flux_kennedy_gruber_souza_etal, u_ll..., aux_ll...,
                                    equations)
        turbo_rr = Trixi.cons2turbo(flux_kennedy_gruber_souza_etal, u_rr..., aux_rr...,
                                    equations)
        flux_left_turbo, flux_right_turbo = Trixi.flux_turbo(flux_kennedy_gruber_souza_etal,
                                                             turbo_ll..., turbo_rr...,
                                                             normal_direction...,
                                                             equations)
        @test SVector(flux_left_turbo) ≈ flux_left
        @test SVector(flux_right_turbo) ≈ flux_right

        # The surface flux and the wave speed estimates are homogeneous of degree one in
        # the normal direction
        for i in 1:2
            @test flux_lmars_combined(u_ll, u_rr, aux_ll, aux_rr, 2.5 * normal_direction,
                                      equations)[i] ≈
                  2.5 * flux_lmars_combined(u_ll, u_rr, aux_ll, aux_rr,
                                      normal_direction, equations)[i]
        end
        @test max_abs_speed_naive(u_ll, u_rr, aux_ll, aux_rr, 2.5 * normal_direction,
                                  equations) ≈
              2.5 * max_abs_speed_naive(u_ll, u_rr, aux_ll, aux_rr, normal_direction,
                                  equations)
        @test max_abs_speed(u_ll, u_rr, aux_ll, aux_rr, 2.5 * normal_direction,
                            equations) ≈
              2.5 * max_abs_speed(u_ll, u_rr, aux_ll, aux_rr, normal_direction, equations)
        @test max_abs_speed(u_ll, u_rr, aux_ll, aux_rr, normal_direction, equations) <=
              max_abs_speed_naive(u_ll, u_rr, aux_ll, aux_rr, normal_direction, equations)
    end

    # Without water, the physical and the volume flux are the ones of the dry equations,
    # whose auxiliary variable is the geopotential. The surface flux is not compared since
    # `flux_lmars_combined` of the dry equations uses a constant speed of sound.
    equations_dry = CompressibleEulerPotentialTemperatureEquations3D(c_p = 1004.0,
                                                                     c_v = 717.0,
                                                                     gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    dry(u) = SVector(u[1], u[2], u[3], u[4], u[5], 0.0, 0.0, 0.0)
    u_ll_dry = SVector(u_ll[1], u_ll[2], u_ll[3], u_ll[4], u_ll[5])
    u_rr_dry = SVector(u_rr[1], u_rr[2], u_rr[3], u_rr[4], u_rr[5])
    aux_ll_dry = SVector(aux_ll[1])
    aux_rr_dry = SVector(aux_rr[1])
    for normal_direction in normal_directions
        @test flux(dry(u_ll), aux_ll, normal_direction, equations)[1:5] ≈
              flux(u_ll_dry, aux_ll_dry, normal_direction, equations_dry)
        for i in 1:2
            @test flux_kennedy_gruber_souza_etal(dry(u_ll), dry(u_rr), aux_ll, aux_rr,
                                                 normal_direction, equations)[i][1:5] ≈
                  flux_kennedy_gruber_souza_etal(u_ll_dry, u_rr_dry, aux_ll_dry, aux_rr_dry,
                                                 normal_direction, equations_dry)[i]
        end
    end

    # The rain falls in the direction of the auxiliary variables: there is no
    # sedimentation through a face whose normal is horizontal, and the full terminal
    # velocity through a face whose normal is vertical
    v_ll = SVector(u_ll[2], u_ll[3], u_ll[4]) / u_ll[1]
    rain = u_ll[8] *
           TrixiAtmo.terminal_velocity_rain(TrixiAtmo.density_dry(u_ll, aux_ll, equations),
                                            u_ll[8],
                                            equations.microphysics)
    vertical = SVector(aux_ll[2], aux_ll[3], aux_ll[4])
    horizontal = SVector(x_ll[2], -x_ll[1], 0.0) / 1.0e5
    @test flux(u_ll, aux_ll, horizontal, equations)[8] ≈ u_ll[8] * sum(v_ll .* horizontal)
    @test flux(u_ll, aux_ll, vertical, equations)[8] ≈
          u_ll[8] * sum(v_ll .* vertical) - rain

    # The rain leaves the domain through the lower boundary but does not enter it through
    # the upper boundary
    flux_bottom = boundary_condition_slip_wall(u_ll, aux_ll, -2 * vertical, x_ll, 0.0,
                                               flux_lmars_combined, equations)
    flux_top = boundary_condition_slip_wall(u_ll, aux_ll, 2 * vertical, x_ll, 0.0,
                                            flux_lmars_combined, equations)
    @test flux_bottom[8] ≈ 2 * rain
    @test flux_bottom[1] ≈ 2 * rain
    @test isapprox(flux_top[8], 0, atol = 1e-12)
    @test isapprox(flux_top[1], 0, atol = 1e-12)

    # `prim2cons` is the inverse of `cons2prim`
    for (u, aux) in ((u_ll, aux_ll), (u_rr, aux_rr))
        @test prim2cons(cons2prim(u, aux, equations), aux, equations) ≈ u
    end
end

@testitem "Unit: Consistency check for Kessler source terms: CEMPTE3D" setup=[Setup] tags=[:unit_source_terms] begin
    equations = CompressibleEulerMoistPotentialTemperatureEquations3D(c_p = 1004.5,
                                                                      c_v = 717.5,
                                                                      gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    @unpack R, R_v, c_p, p_0 = equations
    @unpack latent_heat, condensation_time = equations.microphysics
    aux = SVector(1000.0, 0.0, 0.0, 1.0)

    # Right-hand sides of the mixing ratios `q_v`, `q_c`, and `q_r` and isobaric heating
    # rate at the temperature `T` and the pressure `p`
    function rates(T, p, q_v, q_c, q_r)
        rho_d = p / ((R + R_v * q_v) * T)
        exner = (p / p_0)^(R / c_p)
        u = SVector(rho_d * (1 + q_v + q_c + q_r), 1.0, -2.0, 0.5,
                    rho_d * T / exner * (1 + R_v / R * q_v),
                    rho_d * q_v, rho_d * q_c, rho_d * q_r)
        @test pressure(u, aux, equations) ≈ p
        du = source_terms_microphysics(u, aux, nothing, 0.0, equations)

        # The microphysics only moves water between the species and does not change the
        # total density and the momentum
        @test du[1] == du[2] == du[3] == du[4] == 0
        @test isapprox(du[6] + du[7] + du[8], 0, atol = 1e-15)

        # Heating rate from the right-hand sides of
        # `rho * theta_v = rho_d * theta * (1 + q_v R_v / R)` and of the vapor, which equals
        # the latent heat released by the condensed vapor
        heating = exner / rho_d * (du[5] - R_v / R * T / exner * du[6]) /
                  (1 + R_v / R * q_v)
        @test isapprox(heating, -latent_heat / c_p * du[6] / rho_d, rtol = 1e-12,
                       atol = 1e-15)
        return du[6:8] / rho_d, heating
    end

    for (T, p) in ((285.0, 85000.0), (300.0, 95000.0), (260.0, 60000.0))
        q_vs = TrixiAtmo.saturation_mixing_ratio_tetens(T, p)
        zero_rates = SVector(0.0, 0.0, 0.0)

        # No processes in subsaturated air without condensate, in saturated air with cloud
        # water below the autoconversion threshold of 1 g/kg, and in saturated air with rain
        @test all(iszero, rates(T, p, 0.5 * q_vs, 0.0, 0.0)[1])
        @test isapprox(rates(T, p, q_vs, 5.0e-4, 0.0)[1], zero_rates, atol = 1e-15)
        @test isapprox(rates(T, p, q_vs, 0.0, 2.0e-3)[1], zero_rates, atol = 1e-15)

        # Autoconversion of the cloud water above the threshold
        rates_q, _ = rates(T, p, q_vs, 2.0e-3, 0.0)
        @test isapprox(rates_q[3], 1.0e-3 * (2.0e-3 - 1.0e-3), rtol = 1e-12)
        @test isapprox(rates_q[2], -rates_q[3], rtol = 1e-10)

        # Subsaturated air evaporates the cloud water within `condensation_time` and the
        # rain with latent cooling
        rates_q, _ = rates(T, p, 0.5 * q_vs, 1.0e-4, 0.0)
        @test rates_q[2] ≈ -1.0e-4 / condensation_time
        rates_q, heating = rates(T, p, 0.8 * q_vs, 0.0, 2.0e-3)
        @test rates_q[3] < 0 && rates_q[2] == 0 && heating < 0

        # Condensing for `condensation_time` adjusts supersaturated air to saturation up to
        # second order in the supersaturation, since the condensation accounts for the
        # latent heating with the derivative of `saturation_mixing_ratio_tetens`
        residuals = map((1.0e-3, 1.0e-4)) do supersaturation
            q_v = (1 + supersaturation) * q_vs
            rates_q, heating = rates(T, p, q_v, 0.0, 0.0)
            q_v_new = q_v + condensation_time * rates_q[1]
            T_new = T + condensation_time * heating
            return ((q_v_new - TrixiAtmo.saturation_mixing_ratio_tetens(T_new, p)) /
                    (q_v - q_vs))
        end
        @test isapprox(residuals[1] / residuals[2], 10, rtol = 1e-2)
    end
end
