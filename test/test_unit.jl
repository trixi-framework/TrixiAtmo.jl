@testitem "Unit: Consistency check for EC flux with Potential Temperature: CEPTE" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquations2D(c_p = 1004.0, c_v = 717.0)
    u = SVector(1.1, -0.5, 2.34, 330.0)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        @test flux_ec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    #check consistency between 1D and 2D EC fluxes
    u_1d = SVector(u[1], u[2], u[4])
    u_2d = SVector(u[1], u[2], 0, u[4])
    normal_1d = SVector(-0.3)
    normal_2d = SVector(normal_1d[1], 0.0)
    equations_1d = CompressibleEulerPotentialTemperatureEquations1D(c_p = 1004.0,
                                                                    c_v = 717.0)
    equations_2d = equations
    flux_1d = normal_1d[1] * flux_ec(u_1d, u_1d, 1, equations_1d)
    flux_2d = flux_ec(u_2d, u_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux(u_1d, normal_1d, equations_1d)
    @test flux_1d ≈ flux_2d[[1, 2, 4]]

    # test when u_ll is not the same as u_rr
    u_rr_1d = SVector(2.1, 0.3, 280.5)
    u_rr_2d = SVector(u_rr_1d[1], u_rr_1d[2], 0.0, u_rr_1d[3])
    flux_1d = normal_1d[1] * flux_ec(u_1d, u_rr_1d, 1, equations_1d)
    flux_2d = flux_ec(u_2d, u_rr_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux_2d[[1, 2, 4]]

    # check consistency for 3D EC flux
    equations = CompressibleEulerPotentialTemperatureEquations3D(c_p = 1004.0, c_v = 717.0)
    u = SVector(1.1, -0.5, 2.34, 2.4, 330.0)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        @test flux_ec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    # check consistency between 1D and 3D EC fluxes
    u_3d = SVector(u[1], u[2], 0.0, 0.0, u[5])
    normal_3d = SVector(normal_1d[1], 0.0, 0.0)
    equations_3d = equations
    flux_3d = flux_ec(u_3d, u_3d, normal_3d, equations_3d)
    flux_1d = normal_1d[1] * flux_ec(u_1d, u_1d, 1, equations_1d)
    @test flux_1d ≈ flux_3d[[1, 2, 5]]
end

@testitem "Unit: Consistency check for TEC flux with Potential Temperature: CEPTE" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquations2D(c_p = 1004.0, c_v = 717.0)
    u = SVector(1.1, -0.5, 2.34, 330.0)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        @test flux_tec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    #check consistency between 1D and 2D EC fluxes
    u_1d = SVector(u[1], u[2], u[4])
    u_2d = SVector(u[1], u[2], 0, u[4])
    normal_1d = SVector(-0.3)
    normal_2d = SVector(normal_1d[1], 0.0)
    equations_1d = CompressibleEulerPotentialTemperatureEquations1D(c_p = 1004.0,
                                                                    c_v = 717.0)
    equations_2d = equations
    flux_1d = normal_1d[1] * flux_tec(u_1d, u_1d, 1, equations_1d)
    flux_2d = flux_tec(u_2d, u_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux(u_1d, normal_1d, equations_1d)
    @test flux_1d ≈ flux_2d[[1, 2, 4]]

    # test when u_ll is not the same as u_rr
    u_rr_1d = SVector(2.1, 0.3, 280.5)
    u_rr_2d = SVector(u_rr_1d[1], u_rr_1d[2], 0.0, u_rr_1d[3])
    flux_1d = normal_1d[1] * flux_tec(u_1d, u_rr_1d, 1, equations_1d)
    flux_2d = flux_tec(u_2d, u_rr_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux_2d[[1, 2, 4]]

    # check consistency for 3D EC flux
    equations = CompressibleEulerPotentialTemperatureEquations3D(c_p = 1004.0, c_v = 717.0)
    u = SVector(1.1, -0.5, 2.34, 2.4, 330.0)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        @test flux_tec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    # check consistency between 1D and 3D EC fluxes
    u_3d = SVector(u[1], u[2], 0.0, 0.0, u[5])
    normal_3d = SVector(normal_1d[1], 0.0, 0.0)
    equations_3d = equations
    flux_3d = flux_tec(u_3d, u_3d, normal_3d, equations_3d)
    flux_1d = normal_1d[1] * flux_tec(u_1d, u_1d, 1, equations_1d)
    @test flux_1d ≈ flux_3d[[1, 2, 5]]
end

@testitem "Unit: Consistency check for ETEC flux with Potential Temperature: CEPTE" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquations2D(c_p = 1004.0, c_v = 717.0)
    u = SVector(1.1, -0.5, 2.34, 330.0)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        @test flux_etec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    #check consistency between 1D and 2D EC fluxes
    u_1d = SVector(u[1], u[2], u[4])
    u_2d = SVector(u[1], u[2], 0, u[4])
    normal_1d = SVector(-0.3)
    normal_2d = SVector(normal_1d[1], 0.0)
    equations_1d = CompressibleEulerPotentialTemperatureEquations1D(c_p = 1004.0,
                                                                    c_v = 717.0)
    equations_2d = equations
    flux_1d = normal_1d[1] * flux_etec(u_1d, u_1d, 1, equations_1d)
    flux_2d = flux_etec(u_2d, u_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux(u_1d, normal_1d, equations_1d)
    @test flux_1d ≈ flux_2d[[1, 2, 4]]

    # test when u_ll is not the same as u_rr
    u_rr_1d = SVector(2.1, 0.3, 280.5)
    u_rr_2d = SVector(u_rr_1d[1], u_rr_1d[2], 0.0, u_rr_1d[3])
    flux_1d = normal_1d[1] * flux_etec(u_1d, u_rr_1d, 1, equations_1d)
    flux_2d = flux_etec(u_2d, u_rr_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux_2d[[1, 2, 4]]

    # check consistency for 3D EC flux
    equations = CompressibleEulerPotentialTemperatureEquations3D(c_p = 1004.0, c_v = 717.0)
    u = SVector(1.1, -0.5, 2.34, 2.4, 330.0)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        @test flux_etec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    # check consistency between 1D and 3D EC fluxes
    u_3d = SVector(u[1], u[2], 0.0, 0.0, u[5])
    normal_3d = SVector(normal_1d[1], 0.0, 0.0)
    equations_3d = equations
    flux_3d = flux_etec(u_3d, u_3d, normal_3d, equations_3d)
    flux_1d = normal_1d[1] * flux_etec(u_1d, u_1d, 1, equations_1d)
    @test flux_1d ≈ flux_3d[[1, 2, 5]]
end

@testitem "Unit: Consistency check for LMARS flux with Potential Temperature: CEPTE" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquations2D(c_p = 1004.0, c_v = 717.0)
    flux_lmars = FluxLMARS(340)
    u = SVector(1.1, -0.5, 2.34, 330.0)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        @test flux_lmars(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    #check consistency between 1D and 2D EC fluxes
    u_1d = SVector(u[1], u[2], u[4])
    u_2d = SVector(u[1], u[2], 0, u[4])
    normal_1d = SVector(-0.3)
    normal_2d = SVector(normal_1d[1], 0.0)
    equations_1d = CompressibleEulerPotentialTemperatureEquations1D(c_p = 1004.0,
                                                                    c_v = 717.0)
    equations_2d = equations
    flux_1d = normal_1d[1] * flux_lmars(u_1d, u_1d, 1, equations_1d)
    flux_2d = flux_lmars(u_2d, u_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux(u_1d, normal_1d, equations_1d)
    @test flux_1d ≈ flux_2d[[1, 2, 4]]

    # check consistency for 3D EC flux
    equations = CompressibleEulerPotentialTemperatureEquations3D(c_p = 1004.0, c_v = 717.0)
    u = SVector(1.1, -0.5, 2.34, 2.4, 330.0)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        @test flux_lmars(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    # check consistency between 1D and 3D EC fluxes
    u_3d = SVector(u[1], -u[2], 0.0, 0.0, u[5])
    u_1d = SVector(u[1], -u[2], u[5])
    normal_3d = SVector(normal_1d[1], 0.0, 0.0)
    equations_3d = equations
    flux_3d = flux_lmars(u_3d, u_3d, normal_3d, equations_3d)
    flux_1d = normal_1d[1] * flux_lmars(u_1d, u_1d, 1, equations_1d)
    @test flux_1d ≈ flux_3d[[1, 2, 5]]
end

@testitem "Unit: Consistency check for EC flux with Potential Temperature with gravity: CEPTEWG" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity2D(c_p = 1004.0,
                                                                            c_v = 717.0,
                                                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 330.0, 1500)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        @test flux_ec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    #check consistency between 1D and 2D EC fluxes
    u_1d = SVector(u[1], u[2], u[4], u[5])
    u_2d = SVector(u[1], u[2], 0, u[4], u[5])
    normal_1d = SVector(-0.3)
    normal_2d = SVector(normal_1d[1], 0.0)
    equations_1d = CompressibleEulerPotentialTemperatureEquationsWithGravity1D(c_p = 1004.0,
                                                                               c_v = 717.0,
                                                                               gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    equations_2d = equations
    flux_1d = normal_1d[1] * flux_ec(u_1d, u_1d, 1, equations_1d)
    flux_2d = flux_ec(u_2d, u_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux(u_1d, normal_1d, equations_1d)
    @test flux_1d ≈ flux_2d[[1, 2, 4, 5]]

    # test when u_ll is not the same as u_rr
    u_rr_1d = SVector(2.1, 0.3, 280.5, 1700)
    u_rr_2d = SVector(u_rr_1d[1], u_rr_1d[2], 0.0, u_rr_1d[3], u_rr_1d[4])
    flux_1d = normal_1d[1] * flux_ec(u_1d, u_rr_1d, 1, equations_1d)
    flux_2d = flux_ec(u_2d, u_rr_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux_2d[[1, 2, 4, 5]]

    # check consistency for 3D EC flux
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = 1004.0,
                                                                            c_v = 717.0,
                                                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 2.4, 330.0, 1500)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        @test flux_ec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    # check consistency between 1D and 3D EC fluxes
    u_3d = SVector(u[1], u[2], 0.0, 0.0, u[5], u[6])
    normal_3d = SVector(normal_1d[1], 0.0, 0.0)
    equations_3d = equations
    flux_3d = flux_ec(u_3d, u_3d, normal_3d, equations_3d)
    flux_1d = normal_1d[1] * flux_ec(u_1d, u_1d, 1, equations_1d)
    @test flux_1d ≈ flux_3d[[1, 2, 5, 6]]
end

@testitem "Unit: Consistency check for FluxTurbo EC flux with Potential Temperature with gravity: CEPTEWG" setup=[Setup] tags=[:unit_fluxes] begin
    # The `FluxTurbo` specialization of `(flux_ec, flux_nonconservative_isothermal)`
    # must return the same left and right contributions as the fused evaluation of the
    # conservative and the nonconservative flux.
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = 1004.0,
                                                                            c_v = 717.0,
                                                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    volume_flux = (flux_ec, flux_nonconservative_isothermal)
    u_ll = SVector(1.1, -0.5, 2.34, 2.4, 330.0, 1500.0)
    u_rr = SVector(0.9, 0.3, -1.2, 1.7, 280.5, 1700.0)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    @test Trixi.nturbovars(volume_flux, equations) == Val(7)
    turbo_ll = Trixi.cons2turbo(volume_flux, u_ll..., equations)
    turbo_rr = Trixi.cons2turbo(volume_flux, u_rr..., equations)
    @test length(turbo_ll) == 7
    # rho, v1, v2, v3, theta, p, phi
    @test SVector(turbo_ll[1], turbo_ll[2], turbo_ll[3], turbo_ll[4], turbo_ll[6],
                  turbo_ll[7]) ≈ cons2prim(u_ll, equations)
    @test turbo_ll[5] ≈ u_ll[5] / u_ll[1]

    for normal_direction in normal_directions
        flux_left, flux_right = Trixi.flux_turbo(volume_flux, turbo_ll..., turbo_rr...,
                                                 normal_direction..., equations)
        flux_cons = flux_ec(u_ll, u_rr, normal_direction, equations)
        flux_noncons = flux_nonconservative_isothermal(u_ll, u_rr,
                                                       normal_direction, equations)
        @test SVector(flux_left) ≈ flux_cons + 0.5 * flux_noncons
        @test SVector(flux_right) ≈ flux_cons - 0.5 * flux_noncons
    end
end

@testitem "Unit: Consistency check for TEC flux with Potential Temperature with gravity: CEPTEWG" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity2D(c_p = 1004.0,
                                                                            c_v = 717.0,
                                                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 330.0, 1500)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        @test flux_tec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    #check consistency between 1D and 2D EC fluxes
    u_1d = SVector(u[1], u[2], u[4], u[5])
    u_2d = SVector(u[1], u[2], 0, u[4], u[5])
    normal_1d = SVector(-0.3)
    normal_2d = SVector(normal_1d[1], 0.0)
    equations_1d = CompressibleEulerPotentialTemperatureEquationsWithGravity1D(c_p = 1004.0,
                                                                               c_v = 717.0,
                                                                               gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    equations_2d = equations
    flux_1d = normal_1d[1] * flux_tec(u_1d, u_1d, 1, equations_1d)
    flux_2d = flux_tec(u_2d, u_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux(u_1d, normal_1d, equations_1d)
    @test flux_1d ≈ flux_2d[[1, 2, 4, 5]]

    # test when u_ll is not the same as u_rr
    u_rr_1d = SVector(2.1, 0.3, 280.5, 1700)
    u_rr_2d = SVector(u_rr_1d[1], u_rr_1d[2], 0.0, u_rr_1d[3], u_rr_1d[4])
    flux_1d = normal_1d[1] * flux_tec(u_1d, u_rr_1d, 1, equations_1d)
    flux_2d = flux_tec(u_2d, u_rr_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux_2d[[1, 2, 4, 5]]

    # check consistency for 3D EC flux
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = 1004.0,
                                                                            c_v = 717.0,
                                                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 2.4, 330.0, 1500)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        @test flux_tec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    # check consistency between 1D and 3D EC fluxes
    u_3d = SVector(u[1], u[2], 0.0, 0.0, u[5], u[6])
    normal_3d = SVector(normal_1d[1], 0.0, 0.0)
    equations_3d = equations
    flux_3d = flux_tec(u_3d, u_3d, normal_3d, equations_3d)
    flux_1d = normal_1d[1] * flux_tec(u_1d, u_1d, 1, equations_1d)
    @test flux_1d ≈ flux_3d[[1, 2, 5, 6]]
end

@testitem "Unit: Consistency check for ETEC flux with Potential Temperature with gravity: CEPTEWG" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity2D(c_p = 1004.0,
                                                                            c_v = 717.0,
                                                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 330.0, 1500)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        @test flux_etec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    #check consistency between 1D and 2D EC fluxes
    u_1d = SVector(u[1], u[2], u[4], u[5])
    u_2d = SVector(u[1], u[2], 0, u[4], u[5])
    normal_1d = SVector(-0.3)
    normal_2d = SVector(normal_1d[1], 0.0)
    equations_1d = CompressibleEulerPotentialTemperatureEquationsWithGravity1D(c_p = 1004.0,
                                                                               c_v = 717.0,
                                                                               gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    equations_2d = equations
    flux_1d = normal_1d[1] * flux_etec(u_1d, u_1d, 1, equations_1d)
    flux_2d = flux_etec(u_2d, u_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux(u_1d, normal_1d, equations_1d)
    @test flux_1d ≈ flux_2d[[1, 2, 4, 5]]

    # test when u_ll is not the same as u_rr
    u_rr_1d = SVector(2.1, 0.3, 280.5, 1700)
    u_rr_2d = SVector(u_rr_1d[1], u_rr_1d[2], 0.0, u_rr_1d[3], u_rr_1d[4])
    flux_1d = normal_1d[1] * flux_etec(u_1d, u_rr_1d, 1, equations_1d)
    flux_2d = flux_etec(u_2d, u_rr_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux_2d[[1, 2, 4, 5]]

    # check consistency for 3D EC flux
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = 1004.0,
                                                                            c_v = 717.0,
                                                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 2.4, 330.0, 1500)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        @test flux_etec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    # check consistency between 1D and 3D EC fluxes
    u_3d = SVector(u[1], u[2], 0.0, 0.0, u[5], u[6])
    normal_3d = SVector(normal_1d[1], 0.0, 0.0)
    equations_3d = equations
    flux_3d = flux_etec(u_3d, u_3d, normal_3d, equations_3d)
    flux_1d = normal_1d[1] * flux_etec(u_1d, u_1d, 1, equations_1d)
    @test flux_1d ≈ flux_3d[[1, 2, 5, 6]]
end

@testitem "Unit: Consistency check for LMARS flux with Potential Temperature with gravity: CEPTEWG" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity2D(c_p = 1004.0,
                                                                            c_v = 717.0,
                                                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    flux_lmars = FluxLMARS(340)
    u = SVector(1.1, -0.5, 2.34, 330.0, 1700)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        @test flux_lmars(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    #check consistency between 1D and 2D EC fluxes
    u_1d = SVector(u[1], u[2], u[4], u[5])
    u_2d = SVector(u[1], u[2], 0, u[4], u[5])
    normal_1d = SVector(-0.3)
    normal_2d = SVector(normal_1d[1], 0.0)
    equations_1d = CompressibleEulerPotentialTemperatureEquationsWithGravity1D(c_p = 1004.0,
                                                                               c_v = 717.0,
                                                                               gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    equations_2d = equations
    flux_1d = normal_1d[1] * flux_lmars(u_1d, u_1d, 1, equations_1d)
    flux_2d = flux_lmars(u_2d, u_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux(u_1d, normal_1d, equations_1d)
    @test flux_1d ≈ flux_2d[[1, 2, 4, 5]]

    # check consistency for 3D EC flux
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = 1004.0,
                                                                            c_v = 717.0,
                                                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 2.4, 330.0, 1700)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        @test flux_lmars(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    # check consistency between 1D and 3D EC fluxes
    u_3d = SVector(u[1], u[2], 0.0, 0.0, u[5], u[6])
    normal_3d = SVector(normal_1d[1], 0.0, 0.0)
    equations_3d = equations
    flux_3d = flux_lmars(u_3d, u_3d, normal_3d, equations_3d)
    flux_1d = normal_1d[1] * flux_lmars(u_1d, u_1d, 1, equations_1d)
    @test flux_1d ≈ flux_3d[[1, 2, 5, 6]]
end

@testitem "Unit: Consistency check for ES LMARS flux with Potential Temperature with gravity: CEPTEWG" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity2D(c_p = 1004.0,
                                                                            c_v = 717.0,
                                                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 330.0, 1700)
    u_rr = SVector(0.9, 0.3, -1.2, 280.5, 1500)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        # Consistency with the physical flux
        for i in 1:2
            @test flux_lmars_es(u, u, normal_direction, equations)[i] ≈
                  flux(u, normal_direction, equations)
        end

        flux_left, flux_right = flux_lmars_es(u, u_rr, normal_direction, equations)
        @test flux_left == flux_right
        @test flux_lmars_es(u, u_rr, 2.5 * normal_direction, equations)[1] ≈ 2.5 * flux_left
    end

    # Check consistency for the 3D flux
    equations_3d = CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = 1004.0,
                                                                               c_v = 717.0,
                                                                               gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u_3d = SVector(1.1, -0.5, 2.34, 2.4, 330.0, 1700)
    u_rr_3d = SVector(0.9, 0.3, -1.2, 1.7, 280.5, 1500)

    normal_directions_3d = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions_3d
        for i in 1:2
            @test flux_lmars_es(u_3d, u_3d, normal_direction, equations_3d)[i] ≈
                  flux(u_3d, normal_direction, equations_3d)
        end
        flux_left, flux_right = flux_lmars_es(u_3d, u_rr_3d, normal_direction, equations_3d)
        @test flux_left == flux_right
        @test flux_lmars_es(u_3d, u_rr_3d, 2.5 * normal_direction, equations_3d)[1] ≈
              2.5 * flux_left
    end

    # Check consistency between the 2D and 3D fluxes for different left and right states
    u_3d = SVector(u[1], u[2], u[3], 0.0, u[4], u[5])
    u_rr_3d = SVector(u_rr[1], u_rr[2], u_rr[3], 0.0, u_rr[4], u_rr[5])
    for normal_direction in normal_directions
        normal_3d = SVector(normal_direction[1], normal_direction[2], 0.0)
        flux_2d = flux_lmars_es(u, u_rr, normal_direction, equations)
        flux_3d = flux_lmars_es(u_3d, u_rr_3d, normal_3d, equations_3d)
        for i in 1:2
            @test flux_2d[i] ≈ flux_3d[i][[1, 2, 3, 5, 6]]
        end
    end
end

@testitem "Unit: Consistency check for gravity terms with energy with gravity: CEEEWG" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables states with the same temperature
    equations = CompressibleEulerEnergyEquationsWithGravity3D(c_p = 1004.0, c_v = 717.0,
                                                              gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u_ll = prim2cons(SVector(1.1, 2.0, -1.5, 0.5, 9.0e4, 1700.0), equations)
    u_rr = prim2cons(SVector(0.9, -1.0, 3.0, 1.2, 0.9 / 1.1 * 9.0e4, 2900.0), equations)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    # Both gravity terms coincide for an isothermal state
    for normal_direction in normal_directions
        @test flux_nonconservative_waruszewski_etal(u_ll, u_rr, normal_direction,
                                                    equations) ≈
              flux_nonconservative_isothermal(u_ll, u_rr, normal_direction, equations)
    end

    # Check consistency between the 2D and 3D gravity terms
    equations_2d = CompressibleEulerEnergyEquationsWithGravity2D(c_p = 1004.0, c_v = 717.0,
                                                                 gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u_ll_2d = prim2cons(SVector(1.1, 2.0, -1.5, 9.0e4, 1700.0), equations_2d)
    u_rr_2d = prim2cons(SVector(0.9, -1.0, 3.0, 7.5e4, 2900.0), equations_2d)
    u_ll_3d = prim2cons(SVector(1.1, 2.0, -1.5, 0.0, 9.0e4, 1700.0), equations)
    u_rr_3d = prim2cons(SVector(0.9, -1.0, 3.0, 0.0, 7.5e4, 2900.0), equations)
    for normal_direction in (SVector(1.0, 0.0), SVector(0.0, 1.0), SVector(0.5, -0.5),
                             SVector(-1.2, 0.3))
        normal_direction_3d = SVector(normal_direction[1], normal_direction[2], 0.0)
        for flux_gravity in (flux_nonconservative_waruszewski_etal,
                             flux_nonconservative_isothermal)
            @test flux_gravity(u_ll_2d, u_rr_2d, normal_direction, equations_2d) ≈
                  flux_gravity(u_ll_3d, u_rr_3d, normal_direction_3d,
                               equations)[[1, 2, 3, 5, 6]]
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
    equations = CompressibleEulerVectorInvariantEquationsWithGravity3D(c_p = 1004.0,
                                                                       c_v = 717.0,
                                                                       gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u_ll = SVector(1.1, -0.5, 2.34, 2.4, 330.0, 1500.0)
    u_rr = SVector(0.9, 0.3, -1.2, 1.7, 280.5, 1700.0)

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
                                rho_theta * v_normal, 0.0)
        for i in 1:2
            @test flux_etec(u_ll, u_ll, normal_direction, equations)[i] ≈ flux_physical
        end

        # The nonconservative contribution is antisymmetric
        flux_left, flux_right = flux_etec(u_ll, u_rr, normal_direction, equations)
        flux_left_swapped, flux_right_swapped = flux_etec(u_rr, u_ll, normal_direction,
                                                          equations)
        @test flux_left ≈ flux_right_swapped
        @test flux_right ≈ flux_left_swapped

        # The `FluxTurbo` specialization returns the same contributions
        turbo_ll = Trixi.cons2turbo(flux_etec, u_ll..., equations)
        turbo_rr = Trixi.cons2turbo(flux_etec, u_rr..., equations)
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
        e_v = SVector(ntuple(i -> i == v ? epsilon : 0.0, 6))
        @test isapprox(cons2entropy(u_ll, equations)[v],
                       (entropy(u_ll + e_v, equations) -
                        entropy(u_ll - e_v, equations)) / (2 * epsilon); rtol = 1.0e-6)
    end
end

@testitem "Unit: Consistency check for ES LMARS flux with vector-invariant equations: CEVIEWG" setup=[Setup] tags=[:unit_fluxes] begin
    equations = CompressibleEulerVectorInvariantEquationsWithGravity3D(c_p = 1004.0,
                                                                       c_v = 717.0,
                                                                       gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u_ll = SVector(1.1, -0.5, 2.34, 2.4, 330.0, 1500.0)
    u_rr = SVector(0.9, 0.3, -1.2, 1.7, 280.5, 1500.0)

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
                                rho_theta * v_normal, 0.0)
        for i in 1:2
            @test flux_lmars_es(u_ll, u_ll, normal_direction, equations)[i] ≈ flux_physical
        end

        # The flux is homogeneous of degree one in the normal direction
        for i in 1:2
            @test flux_lmars_es(u_ll, u_rr, 2.5 * normal_direction, equations)[i] ≈
                  2.5 * flux_lmars_es(u_ll, u_rr, normal_direction, equations)[i]
        end
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
    # which carry the geopotential in the state. The surface flux is not compared since
    # `flux_lmars_combined` of the dry equations uses a constant speed of sound.
    equations_dry = CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = 1004.0,
                                                                                c_v = 717.0,
                                                                                gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    dry(u) = SVector(u[1], u[2], u[3], u[4], u[5], 0.0, 0.0, 0.0)
    dry_with_geopotential(u, aux) = SVector(u[1], u[2], u[3], u[4], u[5], aux[1])
    u_ll_dry = dry_with_geopotential(u_ll, aux_ll)
    u_rr_dry = dry_with_geopotential(u_rr, aux_rr)
    for normal_direction in normal_directions
        @test flux(dry(u_ll), aux_ll, normal_direction, equations)[1:5] ≈
              flux(u_ll_dry, normal_direction, equations_dry)[1:5]
        for i in 1:2
            @test flux_kennedy_gruber_souza_etal(dry(u_ll), dry(u_rr), aux_ll, aux_rr,
                                                 normal_direction, equations)[i][1:5] ≈
                  flux_kennedy_gruber_souza_etal(u_ll_dry, u_rr_dry, normal_direction,
                                                 equations_dry)[i][1:5]
        end
    end

    # The rain falls in the direction of the auxiliary variables: there is no
    # sedimentation through a face whose normal is horizontal, and the full terminal
    # velocity through a face whose normal is vertical
    v_ll = SVector(u_ll[2], u_ll[3], u_ll[4]) / u_ll[1]
    rain = u_ll[8] * TrixiAtmo.terminal_velocity_rain(u_ll, equations)
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
end

@testitem "Unit: Consistency check for Kessler source terms: CEMPTE3D" setup=[Setup] tags=[:unit_source_terms] begin
    # Reference values of the unmodified `interface/kessler.f90` of DCMIP2016-v1.0, compiled
    # with `gfortran -fdefault-real-8` and called on a column with two levels 1e12 m apart,
    # such that the sedimentation vanishes. The inputs are
    # `(T, p, rho_d, q_v, q_c, q_r, dt)` and the outputs `(T, q_v, q_c, q_r)` after one
    # call. The states are chosen such that one step of the reference isolates a single
    # process, which makes the rates recoverable exactly:
    # - rain evaporation: no cloud, subsaturated air, and no active limit, i.e.,
    #   `q_r` decreases by `dt * E`;
    # - autoconversion and accretion: supersaturated air without evaporation and the
    #   semi-implicit update `A + C_r = (q_r_ref - q_r) / dt * (1 + 2.2 dt q_r^0.875)`;
    # - condensation: supersaturated air without condensate, where the reference condenses
    #   the amount `condensation_time * C` in one step.
    cases_evap = [
        (245.50998664661992, 73997.41793135805, 1.0499832912908547, 0.00011903200024704985,
         0.0, 0.006393154400063015, 1.0,
         245.5036540052133, 0.00012157265597939809, 0.0, 0.006390613744330666),
        (296.35049717342577, 31587.777324256378, 0.3550477579366977, 0.028626266804346493,
         0.0, 0.003731926498856784, 1.0,
         296.31113128261774, 0.028642060399738673, 0.0, 0.003716132903464603),
        (250.7753430577538, 61818.16354439236, 0.8581864416592901, 0.0005267777206176081,
         0.0, 0.006069668253512719, 1.0,
         250.7690603696752, 0.0005292983350747266, 0.0, 0.0060671476390556),
        (275.6460746683124, 82240.213872993, 1.0347295886276209, 0.002904519781261292, 0.0,
         0.009921625517690312, 1.0,
         275.61172446315226, 0.002918301083571528, 0.0, 0.009907844215380076),
        (286.67191053564613, 97714.05982903691, 1.172832052842216, 0.00785934566403949, 0.0,
         0.0016833001691726457, 1.0,
         286.6661654032873, 0.007861650611141866, 0.0, 0.0016809952220702695),
        (257.37306392152624, 35104.93090225085, 0.47323164035742876, 0.0026538242501123093,
         0.0, 0.0031831304051924674, 1.0,
         257.370362613704, 0.002654908014810588, 0.0, 0.0031820466404941885)
    ]
    cases_acc = [
        (288.71421777755427, 59557.576438472744, 0.6974786440154381, 0.018979358547250663,
         0.0016283580189628585, 0.00584538192872492, 10.0,
         288.9054240925639, 0.0189026465736688, 0.001380044487208449, 0.006170407434061194),
        (264.16208683569846, 44567.65919018755, 0.5836767601416812, 0.0044477511228992244,
         0.0016070674996579966, 0.002704639928955728, 10.0,
         264.2697176949871, 0.004404569622152644, 0.0014667958137769754,
         0.0028880931155833294),
        (248.06498973676273, 51099.97867003122, 0.716618293055866, 0.0009823104325069128,
         0.0005059320530804828, 0.0022687157312843635, 10.0,
         248.07778476812564, 0.0009771770659240987, 0.0004622296595094133,
         0.002317551491438247),
        (239.13887147827748, 36204.20583137108, 0.5270026531816496, 0.0005933243900429469,
         0.0016422105065254248, 0.003270197190982811, 10.0,
         239.14622998890806, 0.0005903721555779585, 0.001428946789064268,
         0.003486413142908956),
        (302.9645922419518, 31769.295210032924, 0.32148422642530144, 0.08489497152304637,
         0.0012371135773060191, 0.0013274845241556088, 10.0,
         303.3104413715413, 0.08475621685225505, 0.00129613145315346,
         0.0014072213190994866),
        (308.4970964719618, 102531.08608987181, 1.0956234592201892, 0.035426673322217925,
         0.0029139658834745077, 0.004995226953357617, 10.0,
         308.61415982317095, 0.0353797075057128, 0.002433198063936951,
         0.0055229605894002946)
    ]
    cases_adj = [
        (249.5156179563099, 86694.81311489329, 1.2093597477695357, 0.0006561386014287879,
         0.0, 0.0, 10.0,
         249.51611783312103, 0.000655938050852163, 2.005505766248733e-7, 0.0),
        (254.23609963123334, 72019.924953978, 0.9851399109581318, 0.0011977363056368782,
         0.0, 0.0, 10.0,
         254.24102853125575, 0.0011957588309478976, 1.9774746889806047e-6, 0.0),
        (300.69145676048333, 46537.86978536745, 0.4992140912135074, 0.049894533548792544,
         0.0, 0.0, 10.0,
         300.7404531400263, 0.04987487620131989, 1.965734747265266e-5, 0.0),
        (272.22984078615616, 83817.9821516052, 1.0654402035149335, 0.004296997070198336,
         0.0, 0.0, 10.0,
         272.24600424373887, 0.004290512291016155, 6.484779182180898e-6, 0.0),
        (265.42754684726765, 26987.417817790632, 0.34977338941256425, 0.007993464987640687,
         0.0, 0.0, 10.0,
         265.4575089225474, 0.007981444203038441, 1.2020784602245963e-5, 0.0),
        (289.8157786024733, 36791.915105501896, 0.420402066170489, 0.03244030564979132, 0.0,
         0.0, 10.0,
         289.8202293979613, 0.03243851999064151, 1.785659149807772e-6, 0.0)
    ]
    equations = CompressibleEulerMoistPotentialTemperatureEquations3D(c_p = 1004.5,
                                                                      c_v = 717.5,
                                                                      gravity = 9.80616,
                                                                      R_v = 461.5)
    aux = SVector(1000.0, 0.0, 0.0, 1.0)
    function state(T, p, rho_d, q_v, q_c, q_r)
        exner = (p / equations.p_0)^(equations.R / equations.c_p)
        rho_theta_v = rho_d * T / exner * (1 + equations.R_v / equations.R * q_v)
        u = SVector(rho_d * (1 + q_v + q_c + q_r), 1.0, -2.0, 0.5, rho_theta_v,
                    rho_d * q_v, rho_d * q_c, rho_d * q_r)
        # The states of the reference satisfy the equation of state
        @test pressure(u, equations) ≈ p
        return u, exner
    end
    # Right-hand side of the dry potential temperature `rho_d * theta` from the ones of
    # `rho * theta_v = rho_d * theta * (1 + q_v R_v / R)` and of the vapor
    function source_rho_theta(du, T, p, q_v)
        exner = (p / equations.p_0)^(equations.R / equations.c_p)
        return (du[5] - equations.R_v / equations.R * T / exner * du[6]) /
               (1 + equations.R_v / equations.R * q_v)
    end

    for (T, p, rho_d, q_v, q_c, q_r, dt, T_ref, q_v_ref, q_c_ref, q_r_ref) in cases_evap
        u, exner = state(T, p, rho_d, q_v, q_c, q_r)
        du = source_terms_kessler(u, aux, nothing, 0.0, equations)
        evaporation = (q_r - q_r_ref) / dt
        @test isapprox(du[8] / rho_d, -evaporation, rtol = 1e-10)
        @test isapprox(du[6] / rho_d, evaporation, rtol = 1e-10)
        # Isobaric latent cooling
        @test isapprox(source_rho_theta(du, T, p, q_v) / rho_d * exner, (T_ref - T) / dt,
                       rtol = 1e-8)
    end

    for (T, p, rho_d, q_v, q_c, q_r, dt, T_ref, q_v_ref, q_c_ref, q_r_ref) in cases_acc
        u, _ = state(T, p, rho_d, q_v, q_c, q_r)
        du = source_terms_kessler(u, aux, nothing, 0.0, equations)
        @test isapprox(du[8] / rho_d, (q_r_ref - q_r) / dt * (1 + dt * 2.2 * q_r^0.875),
                       rtol = 1e-10)
    end

    for (T, p, rho_d, q_v, q_c, q_r, dt, T_ref, q_v_ref, q_c_ref, q_r_ref) in cases_adj
        u, exner = state(T, p, rho_d, q_v, q_c, q_r)
        du = source_terms_kessler(u, aux, nothing, 0.0, equations)
        condensation = du[7] / rho_d
        @test isapprox(equations.condensation_time * condensation, q_c_ref, rtol = 1e-10)
        @test isapprox(equations.condensation_time * source_rho_theta(du, T, p, q_v) /
                       rho_d * exner, T_ref - T, rtol = 1e-8)
    end

    # The microphysics only moves water between the species and does not change the
    # total density and the momentum
    for (T, p, rho_d, q_v, q_c, q_r) in vcat(cases_evap, cases_acc, cases_adj)
        du = source_terms_kessler(state(T, p, rho_d, q_v, q_c, q_r)[1], aux, nothing,
                                  0.0, equations)
        @test isapprox(du[6] + du[7] + du[8], 0, atol = 1e-15)
        @test du[1] == du[2] == du[3] == du[4] == 0
    end
end
