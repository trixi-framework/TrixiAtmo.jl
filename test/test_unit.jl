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
                                 "elixir_unsteady_solid_body_rotation_EC_correction.jl"),
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
                                 "elixir_unsteady_solid_body_rotation_EC.jl"),
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
