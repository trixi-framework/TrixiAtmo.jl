module TestTypeStable

using TrixiAtmo
using Trixi
using Test

include("test_trixiatmo.jl")

# Start with a clean environment: remove Trixi.jl output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

@timed_testset "Compressible Euler Potential Temperature 1D" begin
    for RealT in (Float32, Float64)
        equations = @inferred CompressibleEulerPotentialTemperatureEquations1D(c_p = RealT(1004),
                                                                               c_v = RealT(717))

        x = SVector(zero(RealT))
        t = zero(RealT)
        u = u_ll = u_rr = u_inner = cons = SVector(one(RealT), one(RealT), one(RealT))
        orientation = 1

        surface_flux_function = flux_lax_friedrichs

        @test eltype(@inferred flux(u, orientation, equations)) == RealT
        @test eltype(@inferred flux_ec(u_ll, u_rr, orientation, equations)) ==
              RealT
        @test eltype(@inferred flux_tec(u_ll, u_rr, orientation, equations)) ==
              RealT
        @test eltype(@inferred flux_etec(u_ll, u_rr, orientation, equations)) ==
              RealT

        @test eltype(@inferred cons2prim(u, equations)) == RealT
        @test eltype(@inferred prim2cons(u, equations)) == RealT
        @test eltype(@inferred cons2entropy(u, equations)) == RealT
        @test typeof(@inferred pressure(u, equations)) == RealT
        @test typeof(@inferred entropy(cons, equations)) == RealT
        @test typeof(@inferred energy_kinetic(cons, equations)) == RealT
        @test typeof(@inferred energy_total(cons, equations)) == RealT
    end
end

@timed_testset "Compressible Euler Potential Temperature With Gravity 1D" begin
    for RealT in (Float32, Float64)
        equations = @inferred CompressibleEulerPotentialTemperatureEquationsWithGravity1D(c_p = RealT(1004),
                                                                                          c_v = RealT(717),
                                                                                          gravity = RealT(EARTH_GRAVITATIONAL_ACCELERATION))

        x = SVector(zero(RealT))
        t = zero(RealT)
        u = u_ll = u_rr = u_inner = cons = SVector(one(RealT), one(RealT), one(RealT),
                                                   one(RealT))
        orientation = 1
        directions = [1, 2]

        surface_flux_function = (flux_lax_friedrichs, flux_zero)

        for direction in directions
            @test eltype(@inferred boundary_condition_slip_wall(u_inner, orientation,
                                                                direction,
                                                                x, t,
                                                                surface_flux_function,
                                                                equations)) ==
                  SVector{4, RealT}
        end

        @test eltype(@inferred flux(u, orientation, equations)) == RealT
        @test eltype(@inferred flux_ec(u_ll, u_rr, orientation, equations)) ==
              RealT
        @test eltype(@inferred flux_tec(u_ll, u_rr, orientation, equations)) ==
              RealT
        @test eltype(@inferred flux_etec(u_ll, u_rr, orientation, equations)) ==
              RealT
        @test eltype(@inferred flux_nonconservative_waruszewski_etal(u_ll, u_rr,
                                                                     orientation,
                                                                     equations)) ==
              RealT
        @test eltype(@inferred flux_nonconservative_artiano_etal(u_ll, u_rr,
                                                                 orientation,
                                                                 equations)) ==
              RealT
        @test eltype(@inferred flux_nonconservative_souza_etal(u_ll, u_rr, orientation,
                                                               equations)) ==
              RealT
        flux_lmars = FluxLMARS(RealT(340))
        @test eltype(@inferred flux_lmars(u_ll, u_rr, orientation, equations)) ==
              RealT

        @test eltype(@inferred cons2prim(u, equations)) == RealT
        @test eltype(@inferred prim2cons(u, equations)) == RealT
        @test eltype(@inferred cons2entropy(u, equations)) == RealT
        @test typeof(@inferred pressure(u, equations)) == RealT
        @test typeof(@inferred entropy(cons, equations)) == RealT
        @test typeof(@inferred energy_kinetic(cons, equations)) == RealT
        @test typeof(@inferred energy_total(cons, equations)) == RealT
        @test eltype(@inferred Trixi.max_abs_speeds(u, equations)) == RealT
        @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, orientation, equations)) ==
              RealT
    end
end

@timed_testset "Compressible Euler Potential Temperature 2D" begin
    for RealT in (Float32, Float64)
        equations = @inferred CompressibleEulerPotentialTemperatureEquations2D(c_p = RealT(1004),
                                                                               c_v = RealT(717))

        x = SVector(zero(RealT))
        t = zero(RealT)
        u = u_ll = u_rr = u_inner = cons = SVector(one(RealT), one(RealT), one(RealT),
                                                   one(RealT))

        normal_direction = SVector(one(RealT), one(RealT))
        surface_flux_function = flux_lax_friedrichs
        orientations = [1, 2]
        directions = [1, 2]

        for direction in directions, orientation in orientations
            @test eltype(@inferred boundary_condition_slip_wall(u_inner, orientation,
                                                                direction,
                                                                x, t,
                                                                surface_flux_function,
                                                                equations)) == RealT
            @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, orientation,
                                                       equations)) == RealT
        end

        @test eltype(@inferred flux(u, normal_direction, equations)) == RealT
        @test eltype(@inferred flux_ec(u_ll, u_rr, normal_direction, equations)) ==
              RealT
        @test eltype(@inferred flux_tec(u_ll, u_rr, normal_direction, equations)) ==
              RealT
        @test eltype(@inferred flux_etec(u_ll, u_rr, normal_direction, equations)) ==
              RealT

        @test eltype(@inferred cons2prim(u, equations)) == RealT
        @test eltype(@inferred prim2cons(u, equations)) == RealT
        @test eltype(@inferred cons2entropy(u, equations)) == RealT
        @test typeof(@inferred pressure(u, equations)) == RealT
        @test typeof(@inferred entropy(cons, equations)) == RealT
        @test typeof(@inferred energy_kinetic(cons, equations)) == RealT
        @test typeof(@inferred energy_total(cons, equations)) == RealT
        @test eltype(@inferred Trixi.max_abs_speeds(u, equations)) == RealT
        @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, normal_direction,
                                                   equations)) == RealT
        @test typeof(@inferred max_abs_speed(u_ll, u_rr, normal_direction,
                                             equations)) == RealT
    end
end

@timed_testset "Compressible Euler Potential Temperature With Gravity 2D" begin
    for RealT in (Float32, Float64)
        equations = @inferred CompressibleEulerPotentialTemperatureEquationsWithGravity2D(c_p = RealT(1004),
                                                                                          c_v = RealT(717),
                                                                                          gravity = RealT(EARTH_GRAVITATIONAL_ACCELERATION))

        x = SVector(zero(RealT))
        t = zero(RealT)
        u = u_ll = u_rr = u_inner = cons = SVector(one(RealT), one(RealT), one(RealT),
                                                   one(RealT), one(RealT))

        normal_direction = SVector(one(RealT), one(RealT))
        surface_flux_function = (flux_lax_friedrichs, flux_zero)
        orientations = [1, 2]
        directions = [1, 2]

        for direction in directions, orientation in orientations
            @test eltype(@inferred boundary_condition_slip_wall(u_inner, orientation,
                                                                direction,
                                                                x, t,
                                                                surface_flux_function,
                                                                equations)) ==
                  SVector{5, RealT}
            @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, orientation,
                                                       equations)) == RealT
            @test typeof(@inferred max_abs_speed(u_ll, u_rr, orientation,
                                                 equations)) == RealT
        end
        @test eltype(@inferred flux(u, normal_direction, equations)) == RealT
        @test eltype(@inferred flux_ec(u_ll, u_rr, normal_direction, equations)) ==
              RealT
        @test eltype(@inferred flux_tec(u_ll, u_rr, normal_direction, equations)) ==
              RealT
        @test eltype(@inferred flux_etec(u_ll, u_rr, normal_direction, equations)) ==
              RealT

        @test eltype(@inferred cons2prim(u, equations)) == RealT
        @test eltype(@inferred prim2cons(u, equations)) == RealT
        @test eltype(@inferred cons2entropy(u, equations)) == RealT
        @test typeof(@inferred pressure(u, equations)) == RealT
        @test typeof(@inferred entropy(cons, equations)) == RealT
        @test typeof(@inferred energy_kinetic(cons, equations)) == RealT
        @test typeof(@inferred energy_total(cons, equations)) == RealT
        @test eltype(@inferred Trixi.max_abs_speeds(u, equations)) == RealT
        @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, normal_direction,
                                                   equations)) == RealT
        @test typeof(@inferred max_abs_speed(u_ll, u_rr, normal_direction,
                                             equations)) == RealT
    end
end

@timed_testset "Compressible Euler Potential Temperature 3D" begin
    for RealT in (Float32, Float64)
        equations = @inferred CompressibleEulerPotentialTemperatureEquations3D(c_p = RealT(1004),
                                                                               c_v = RealT(717))

        x = SVector(zero(RealT))
        t = zero(RealT)
        u = u_ll = u_rr = u_inner = cons = SVector(one(RealT), one(RealT), one(RealT),
                                                   one(RealT), one(RealT))

        surface_flux_function = flux_lax_friedrichs
        orientations = [1, 2, 3]

        for orientation in orientations
            @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, orientation,
                                                       equations)) == RealT
            @test typeof(@inferred max_abs_speed(u_ll, u_rr, orientation,
                                                 equations)) == RealT
        end
        normal_direction = SVector(one(RealT), one(RealT), one(RealT))
        @test eltype(@inferred flux(u, normal_direction, equations)) == RealT
        @test eltype(@inferred flux_ec(u_ll, u_rr, normal_direction, equations)) ==
              RealT
        @test eltype(@inferred flux_tec(u_ll, u_rr, normal_direction, equations)) ==
              RealT
        @test eltype(@inferred flux_etec(u_ll, u_rr, normal_direction, equations)) ==
              RealT

        @test eltype(@inferred cons2prim(u, equations)) == RealT
        @test eltype(@inferred prim2cons(u, equations)) == RealT
        @test eltype(@inferred cons2entropy(u, equations)) == RealT
        @test typeof(@inferred pressure(u, equations)) == RealT
        @test typeof(@inferred entropy(cons, equations)) == RealT
        @test typeof(@inferred energy_kinetic(cons, equations)) == RealT
        @test typeof(@inferred energy_total(cons, equations)) == RealT
        @test eltype(@inferred Trixi.max_abs_speeds(u, equations)) == RealT
        @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, normal_direction,
                                                   equations)) == RealT
        @test typeof(@inferred max_abs_speed(u_ll, u_rr, normal_direction,
                                             equations)) == RealT
    end
end
@timed_testset "Compressible Euler Potential Temperature With Gravity 3D" begin
    for RealT in (Float32, Float64)
        equations = @inferred CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = RealT(1004),
                                                                                          c_v = RealT(717),
                                                                                          gravity = RealT(EARTH_GRAVITATIONAL_ACCELERATION))

        x = SVector(zero(RealT))
        t = zero(RealT)
        u = u_ll = u_rr = u_inner = cons = SVector(one(RealT), one(RealT), one(RealT),
                                                   one(RealT), one(RealT), one(RealT))

        normal_direction = SVector(one(RealT), one(RealT), one(RealT))

        orientations = [1, 2, 3]

        for orientation in orientations
            @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, orientation,
                                                       equations)) == RealT
            @test typeof(@inferred max_abs_speed(u_ll, u_rr, orientation,
                                                 equations)) == RealT
        end
        @test eltype(@inferred flux(u, normal_direction, equations)) == RealT
        @test eltype(@inferred flux_ec(u_ll, u_rr, normal_direction, equations)) ==
              RealT
        @test eltype(@inferred flux_tec(u_ll, u_rr, normal_direction, equations)) ==
              RealT
        @test eltype(@inferred flux_etec(u_ll, u_rr, normal_direction, equations)) ==
              RealT
        @test eltype(@inferred cons2prim(u, equations)) == RealT
        @test eltype(@inferred prim2cons(u, equations)) == RealT
        @test eltype(@inferred cons2entropy(u, equations)) == RealT
        @test typeof(@inferred pressure(u, equations)) == RealT
        @test typeof(@inferred entropy(cons, equations)) == RealT
        @test typeof(@inferred energy_kinetic(cons, equations)) == RealT
        @test typeof(@inferred energy_total(cons, equations)) == RealT
        @test eltype(@inferred Trixi.max_abs_speeds(u, equations)) == RealT
        @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, normal_direction,
                                                   equations)) == RealT
        @test typeof(@inferred max_abs_speed(u_ll, u_rr, normal_direction,
                                             equations)) == RealT
    end
end

@timed_testset "Compressible Euler Moist Euler 2D" begin
    for RealT in (Float32, Float64)
        equations = @inferred CompressibleMoistEulerEquations2D(c_pd = RealT(2.1),
                                                                c_vd = RealT(2),
                                                                c_pv = RealT(2.1),
                                                                c_vv = RealT(2),
                                                                gravity = RealT(1))

        x = SVector(zero(RealT))
        t = zero(RealT)
        u = u_ll = u_rr = u_inner = cons = SVector(one(RealT), one(RealT), one(RealT),
                                                   one(RealT), one(RealT), one(RealT))

        normal_direction = SVector(one(RealT), one(RealT))
        surface_flux_function = FluxLMARS(RealT(340))
        directions = [1, 2]

        for direction in directions
            @test eltype(@inferred boundary_condition_slip_wall(u_inner,
                                                                normal_direction,
                                                                direction,
                                                                x, t,
                                                                surface_flux_function,
                                                                equations)) == RealT
        end
        @test eltype(@inferred flux(u, normal_direction, equations)) == RealT
        @test eltype(@inferred flux_chandrashekar(u_ll, u_rr, normal_direction,
                                                  equations)) == RealT

        @test eltype(@inferred TrixiAtmo.cons2temp(u, equations)) == RealT
        @test eltype(@inferred TrixiAtmo.cons2drypot(u, equations)) == RealT
        @test eltype(@inferred TrixiAtmo.cons2moistpot(u, equations)) == RealT
        @test eltype(@inferred cons2prim(u, equations)) == RealT
        @test eltype(@inferred TrixiAtmo.moist_pottemp_thermodynamic(u, equations)) ==
              RealT
        @test eltype(@inferred TrixiAtmo.dry_pottemp_thermodynamic(u, equations)) ==
              RealT
        @test eltype(@inferred prim2cons(u, equations)) == RealT
        @test eltype(@inferred TrixiAtmo.density(u, equations)) == RealT
        @test eltype(@inferred TrixiAtmo.density_dry(u, equations)) == RealT
        @test eltype(@inferred TrixiAtmo.density_vapor(u, equations)) == RealT
        @test eltype(@inferred TrixiAtmo.temperature(u, equations)) == RealT
        @test eltype(@inferred TrixiAtmo.density_liquid(u, equations)) == RealT
        @test eltype(@inferred TrixiAtmo.ratio_liquid(u, equations)) == RealT
        @test eltype(@inferred TrixiAtmo.ratio_vapor(u, equations)) == RealT
        @test eltype(@inferred TrixiAtmo.density_pressure(u, equations)) == RealT
        @test typeof(@inferred pressure(u, equations)) == RealT
        @test eltype(@inferred energy_internal(u, equations)) == RealT
        @test typeof(@inferred energy_kinetic(cons, equations)) == RealT
        @test typeof(@inferred energy_total(cons, equations)) == RealT
        @test eltype(@inferred Trixi.max_abs_speeds(u, equations)) == RealT
        @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, normal_direction,
                                                   equations)) == RealT
    end
end

@timed_testset "Compressible Rainy Euler" begin
    using TrixiAtmo: boundary_condition_simple_slip_wall, cons2eq_pot_temp,
                     cons2nonlinearsystemsol, cons2speeds, densities, velocities,
                     energy_density, speed_of_sound, terminal_velocity_rain,
                     saturation_vapour_pressure, saturation_vapour_pressure_derivative,
                     saturation_residual, saturation_residual_jacobian, theta_d,
                     AtmosphereLayersRainyBubble
    for RealT in (Float32, Float64)
        equations = @inferred CompressibleRainyEulerEquations2D(RealT)

        x = SVector(zero(RealT))
        t = zero(RealT)
        u = u_ll = u_rr = SVector(ntuple(i -> one(RealT), 9))

        normal_direction = SVector(one(RealT), one(RealT))
        surface_flux_function = flux_lax_friedrichs
        orientations = [1, 2]
        directions = [1, 2]

        for direction in directions, orientation in orientations
            @test eltype(@inferred boundary_condition_simple_slip_wall(u, orientation,
                                                                       direction,
                                                                       x, t,
                                                                       surface_flux_function,
                                                                       equations)) ==
                  RealT
        end
        for direction in directions
            @test eltype(@inferred boundary_condition_slip_wall(u, normal_direction,
                                                                direction,
                                                                x, t,
                                                                surface_flux_function,
                                                                equations)) == RealT
        end
        for orientation in orientations
            @test eltype(@inferred max_abs_speed(u_ll, u_rr, orientation, equations)) ==
                  RealT
            @test eltype(@inferred max_abs_speed_naive(u_ll, u_rr, orientation,
                                                       equations)) == RealT
            @test eltype(@inferred flux(u, orientation, equations)) == RealT
            @test eltype(@inferred flux_ec_rain(u_ll, u_rr, orientation, equations)) ==
                  RealT
        end

        @test eltype(@inferred flux(u, normal_direction, equations)) == RealT
        @test eltype(@inferred flux_ec_rain(u_ll, u_rr, normal_direction, equations)) ==
              RealT
        @test eltype(@inferred flux_LMARS(u_ll, u_rr, normal_direction, equations)) ==
              RealT
        @test eltype(@inferred cons2prim(u, equations)) == RealT
        @test eltype(@inferred cons2entropy(u, equations)) == RealT
        @test eltype(@inferred cons2eq_pot_temp(u, equations)) == RealT
        @test eltype(@inferred cons2nonlinearsystemsol(u, equations)) == RealT
        @test eltype(@inferred cons2speeds(u, equations)) == RealT
        @test eltype(@inferred densities(u, equations)) == RealT
        @test eltype(@inferred velocities(u, u[1], equations)) == RealT
        @test eltype(@inferred energy_density(u, equations)) == RealT
        @test eltype(@inferred pressure(u, equations)) == RealT
        @test eltype(@inferred speed_of_sound(u, equations)) == RealT
        @test eltype(@inferred terminal_velocity_rain(u[2], u[3], equations)) == RealT
        @test eltype(@inferred saturation_vapour_pressure(u[9], equations)) == RealT
        @test eltype(@inferred saturation_vapour_pressure_derivative(u[9], equations)) ==
              RealT
        @test eltype(@inferred Trixi.max_abs_speeds(u, equations)) == RealT
        @test eltype(@inferred boundary_condition_slip_wall(u, normal_direction,
                                                            x, t,
                                                            surface_flux_function,
                                                            equations)) == RealT
        @test eltype(@inferred saturation_residual(u, u, equations)) == RealT
        @test eltype(@inferred saturation_residual_jacobian(u, u, equations)) == RealT
        @test eltype(@inferred max_abs_speed_naive(u_ll, u_rr, normal_direction,
                                                   equations)) == RealT
        @test eltype(@inferred theta_d(x[1], equations)) == RealT

        @inferred AtmosphereLayersRainyBubble(equations; total_height = x[1])
    end
end

@timed_testset "Shallow Water 3D" begin
    for RealT in (Float32, Float64)
        equations = @inferred ShallowWaterEquations3D(gravity = RealT(1),
                                                      rotation_rate = RealT(1),
                                                      H0 = RealT(1))

        x = SVector(zero(RealT))
        t = zero(RealT)
        u = u_ll = u_rr = u_inner = cons = SVector(one(RealT), one(RealT), one(RealT),
                                                   one(RealT), one(RealT))

        normal_direction = SVector(one(RealT), one(RealT), one(RealT))
        surface_flux_function = FluxLMARS(RealT(340))
        directions = [1, 2]
        @test eltype(@inferred flux(u, normal_direction, equations)) == RealT
        @test eltype(@inferred flux_wintermeyer_etal(u_ll, u_rr, normal_direction,
                                                     equations)) == RealT
        @test eltype(@inferred flux_fjordholm_etal(u_ll, u_rr, normal_direction,
                                                   equations)) == RealT
        @test eltype(@inferred flux_nonconservative_wintermeyer_etal(u_ll, u_rr,
                                                                     normal_direction,
                                                                     equations)) ==
              RealT
        @test eltype(@inferred flux_nonconservative_fjordholm_etal(u_ll, u_rr,
                                                                   normal_direction,
                                                                   equations)) == RealT
    end
end

@timed_testset "Compressible Euler Energy With Gravity 2D" begin
    for RealT in (Float32, Float64)
        equations = @inferred CompressibleEulerEnergyEquationsWithGravity2D(c_p = RealT(1004),
                                                                            c_v = RealT(717),
                                                                            gravity = RealT(9.81))

        x = SVector(zero(RealT))
        t = zero(RealT)
        u = u_ll = u_rr = u_inner = cons = SVector(one(RealT), one(RealT), one(RealT),
                                                   one(RealT), one(RealT))

        normal_direction = SVector(one(RealT), one(RealT))
        surface_flux_function = (flux_lax_friedrichs, flux_zero)
        orientations = [1, 2]
        directions = [1, 2]

        for direction in directions, orientation in orientations
            @test eltype(@inferred boundary_condition_slip_wall(u_inner, orientation,
                                                                direction,
                                                                x, t,
                                                                surface_flux_function,
                                                                equations)) ==
                  SVector{5, RealT}
            @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, orientation,
                                                       equations)) == RealT
            @test typeof(@inferred max_abs_speed(u_ll, u_rr, orientation,
                                                 equations)) == RealT
        end
        @test eltype(@inferred flux(u, normal_direction, equations)) == RealT
        @test eltype(@inferred flux_ranocha(u_ll, u_rr, normal_direction, equations)) ==
              RealT
        @test eltype(@inferred flux_kennedy_gruber(u_ll, u_rr, normal_direction,
                                                   equations)) ==
              RealT
        @test eltype(@inferred flux_shima_etal(u_ll, u_rr, normal_direction, equations)) ==
              RealT
        @test eltype(@inferred flux_nonconservative_artiano_etal(u_ll, u_rr,
                                                                 normal_direction,
                                                                 equations)) ==
              RealT
        @test eltype(@inferred flux_nonconservative_waruszewski_etal(u_ll, u_rr,
                                                                     normal_direction,
                                                                     equations)) ==
              RealT
        @test eltype(@inferred flux_nonconservative_souza_etal(u_ll, u_rr,
                                                               normal_direction,
                                                               equations)) ==
              RealT
        @test eltype(@inferred boundary_condition_slip_wall(u_inner,
                                                            normal_direction,
                                                            x, t,
                                                            surface_flux_function,
                                                            equations)) ==
              SVector{5, RealT}
        @test eltype(varnames(cons2prim, equations)) == String

        @test eltype(@inferred cons2prim(u, equations)) == RealT
        @test eltype(@inferred prim2cons(u, equations)) == RealT
        @test eltype(@inferred cons2entropy(u, equations)) == RealT
        @test eltype(@inferred entropy2cons(cons2entropy(u, equations), equations)) ==
              RealT
        @test typeof(@inferred pressure(u, equations)) == RealT
        @test typeof(@inferred entropy(cons, equations)) == RealT
        @test typeof(@inferred energy_kinetic(cons, equations)) == RealT
        @test typeof(@inferred energy_total(cons, equations)) == RealT
        @test typeof(@inferred energy_internal(cons, equations)) == RealT
        @test eltype(@inferred Trixi.max_abs_speeds(u, equations)) == RealT
        @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, normal_direction,
                                                   equations)) == RealT
        @test typeof(@inferred max_abs_speed(u_ll, u_rr, normal_direction,
                                             equations)) == RealT
    end
end

@timed_testset "Compressible Euler Energy With Gravity 3D" begin
    for RealT in (Float32, Float64)
        equations = @inferred CompressibleEulerEnergyEquationsWithGravity3D(c_p = RealT(1004),
                                                                            c_v = RealT(717),
                                                                            gravity = RealT(9.81))

        x = SVector(zero(RealT))
        t = zero(RealT)
        u = u_ll = u_rr = u_inner = cons = SVector(one(RealT), one(RealT), one(RealT),
                                                   one(RealT), RealT(2), one(RealT))

        normal_direction = SVector(one(RealT), one(RealT), one(RealT))

        orientations = [1, 2, 3]

        for orientation in orientations
            @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, orientation,
                                                       equations)) == RealT
            @test typeof(@inferred max_abs_speed(u_ll, u_rr, orientation,
                                                 equations)) == RealT
        end
        @test eltype(@inferred flux(u, normal_direction, equations)) == RealT
        @test eltype(@inferred flux_ranocha(u_ll, u_rr, normal_direction, equations)) ==
              RealT
        @test eltype(@inferred flux_kennedy_gruber(u_ll, u_rr, normal_direction,
                                                   equations)) ==
              RealT
        @test eltype(@inferred flux_shima_etal(u_ll, u_rr, normal_direction, equations)) ==
              RealT
        @test eltype(@inferred flux_chandrashekar(u_ll, u_rr, normal_direction,
                                                  equations)) == RealT
        @test eltype(@inferred flux_nonconservative_artiano_etal(u_ll, u_rr,
                                                                 normal_direction,
                                                                 equations)) ==
              RealT
        @test eltype(@inferred flux_nonconservative_waruszewski_etal(u_ll, u_rr,
                                                                     normal_direction,
                                                                     equations)) ==
              RealT
        @test eltype(@inferred flux_nonconservative_souza_etal(u_ll, u_rr,
                                                               normal_direction,
                                                               equations)) ==
              RealT
        @test eltype(@inferred cons2prim(u, equations)) == RealT
        @test eltype(@inferred prim2cons(u, equations)) == RealT
        @test eltype(@inferred cons2entropy(u, equations)) == RealT
        @test eltype(@inferred entropy2cons(cons2entropy(u, equations), equations)) ==
              RealT
        @test typeof(@inferred pressure(u, equations)) == RealT
        @test typeof(@inferred entropy(cons, equations)) == RealT
        @test typeof(@inferred energy_kinetic(cons, equations)) == RealT
        @test typeof(@inferred energy_total(cons, equations)) == RealT
        @test typeof(@inferred energy_internal(cons, equations)) == RealT
        @test eltype(@inferred Trixi.max_abs_speeds(u, equations)) == RealT
        @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, normal_direction,
                                                   equations)) == RealT
        @test typeof(@inferred max_abs_speed(u_ll, u_rr, normal_direction,
                                             equations)) == RealT
    end
end
end
