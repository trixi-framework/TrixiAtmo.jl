module TestUnit

using TrixiAtmo

include("test_trixiatmo.jl")

# Start with a clean environment: remove Trixi.jl output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

@timed_testset "Consistency check for EC flux with Potential Temperature: CEPTE" begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquations2D()
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
    equations_1d = CompressibleEulerPotentialTemperatureEquations1D()
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
    equations = CompressibleEulerPotentialTemperatureEquations3D()
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

@timed_testset "Consistency check for TEC flux with Potential Temperature: CEPTE" begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquations2D()
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
    equations_1d = CompressibleEulerPotentialTemperatureEquations1D()
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
    equations = CompressibleEulerPotentialTemperatureEquations3D()
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

@timed_testset "Consistency check for ETEC flux with Potential Temperature: CEPTE" begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquations2D()
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
    equations_1d = CompressibleEulerPotentialTemperatureEquations1D()
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
    equations = CompressibleEulerPotentialTemperatureEquations3D()
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

@timed_testset "Consistency check for LMARS flux with Potential Temperature: CEPTE" begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquations2D()
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
    equations_1d = CompressibleEulerPotentialTemperatureEquations1D()
    equations_2d = equations
    flux_1d = normal_1d[1] * flux_lmars(u_1d, u_1d, 1, equations_1d)
    flux_2d = flux_lmars(u_2d, u_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux(u_1d, normal_1d, equations_1d)
    @test flux_1d ≈ flux_2d[[1, 2, 4]]

    # check consistency for 3D EC flux
    equations = CompressibleEulerPotentialTemperatureEquations3D()
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
    u_3d = SVector(u[1], u[2], 0.0, 0.0, u[5])
    normal_3d = SVector(normal_1d[1], 0.0, 0.0)
    equations_3d = equations
    flux_3d = flux_lmars(u_3d, u_3d, normal_3d, equations_3d)
    flux_1d = normal_1d[1] * flux_lmars(u_1d, u_1d, 1, equations_1d)
    @test flux_1d ≈ flux_3d[[1, 2, 5]]
end

@timed_testset "Consistency check for EC flux with Potential Temperature with gravity: CEPTE" begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity2D()
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
    equations_1d = CompressibleEulerPotentialTemperatureEquationsWithGravity1D()
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
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity3D()
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

@timed_testset "Consistency check for TEC flux with Potential Temperature with gravity: CEPTE" begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity2D()
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
    equations_1d = CompressibleEulerPotentialTemperatureEquationsWithGravity1D()
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
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity3D()
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

@timed_testset "Consistency check for ETEC flux with Potential Temperature with gravity: CEPTE" begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity2D()
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
    equations_1d = CompressibleEulerPotentialTemperatureEquationsWithGravity1D()
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
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity3D()
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

@timed_testset "Consistency check for LMARS flux with Potential Temperature with gravity: CEPTE" begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity2D()
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
    equations_1d = CompressibleEulerPotentialTemperatureEquationsWithGravity1D()
    equations_2d = equations
    flux_1d = normal_1d[1] * flux_lmars(u_1d, u_1d, 1, equations_1d)
    flux_2d = flux_lmars(u_2d, u_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux(u_1d, normal_1d, equations_1d)
    @test flux_1d ≈ flux_2d[[1, 2, 4, 5]]

    # check consistency for 3D EC flux
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity3D()
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
end
