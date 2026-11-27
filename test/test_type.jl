@testitem "Type stability: Compressible Euler Potential Temperature 1D" setup=[Setup] tags=[:type_stable_tests] begin
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

        # Adapt the real type, e.g., to move a semidiscretization to a GPU
        equations_float64 = @inferred similar(equations, Float64)
        @test typeof(equations_float64.c_p) == Float64
        @test typeof(@inferred similar(equations_float64, RealT)) == typeof(equations)
    end
end

@testitem "Type stability: Compressible Euler Potential Temperature With Gravity 1D" setup=[Setup] tags=[:type_stable_tests] begin
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
        @test eltype(@inferred flux_nonconservative_isothermal(u_ll, u_rr,
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

        # Adapt the real type, e.g., to move a semidiscretization to a GPU
        equations_float64 = @inferred similar(equations, Float64)
        @test typeof(equations_float64.c_p) == Float64
        @test typeof(equations_float64.g) == Float64
        @test typeof(@inferred similar(equations_float64, RealT)) == typeof(equations)
    end
end

@testitem "Type stability: Compressible Euler Potential Temperature 2D" setup=[Setup] tags=[:type_stable_tests] begin
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

        # Adapt the real type, e.g., to move a semidiscretization to a GPU
        equations_float64 = @inferred similar(equations, Float64)
        @test typeof(equations_float64.c_p) == Float64
        @test typeof(@inferred similar(equations_float64, RealT)) == typeof(equations)
    end
end

@testitem "Type stability: Compressible Euler Potential Temperature With Gravity 2D" setup=[Setup] tags=[:type_stable_tests] begin
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
        for i in 1:2
            @test eltype(@inferred flux_lmars_es(u_ll, u_rr, normal_direction, equations)[i]) ==
                  RealT
        end

        @test eltype(@inferred cons2prim(u, equations)) == RealT
        @test eltype(@inferred prim2cons(u, equations)) == RealT
        @test eltype(@inferred cons2entropy(u, equations)) == RealT
        @test typeof(@inferred pressure(u, equations)) == RealT
        @test typeof(@inferred entropy(cons, equations)) == RealT
        @test typeof(@inferred energy_kinetic(cons, equations)) == RealT
        @test typeof(@inferred energy_internal_kinetic(cons, equations)) == RealT
        @test eltype(@inferred Trixi.max_abs_speeds(u, equations)) == RealT
        @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, normal_direction,
                                                   equations)) == RealT
        @test typeof(@inferred max_abs_speed(u_ll, u_rr, normal_direction,
                                             equations)) == RealT

        # Adapt the real type, e.g., to move a semidiscretization to a GPU
        equations_float64 = @inferred similar(equations, Float64)
        @test typeof(equations_float64.c_p) == Float64
        @test typeof(equations_float64.g) == Float64
        @test typeof(@inferred similar(equations_float64, RealT)) == typeof(equations)
    end
end

@testitem "Type stability: Compressible Euler Potential Temperature 3D" setup=[Setup] tags=[:type_stable_tests] begin
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

        # Adapt the real type, e.g., to move a semidiscretization to a GPU
        equations_float64 = @inferred similar(equations, Float64)
        @test typeof(equations_float64.c_p) == Float64
        @test typeof(@inferred similar(equations_float64, RealT)) == typeof(equations)
    end
end
@testitem "Type stability: Compressible Euler Potential Temperature With Gravity 3D" setup=[Setup] tags=[:type_stable_tests] begin
    for RealT in (Float32, Float64)
        equations = @inferred CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = RealT(1004),
                                                                                          c_v = RealT(717),
                                                                                          gravity = RealT(EARTH_GRAVITATIONAL_ACCELERATION))

        x = SVector(zero(RealT))
        t = zero(RealT)
        u = u_ll = u_rr = u_inner = cons = SVector(one(RealT), one(RealT), one(RealT),
                                                   one(RealT), one(RealT))
        aux = aux_ll = aux_rr = SVector(one(RealT))

        normal_direction = SVector(one(RealT), one(RealT), one(RealT))

        orientations = [1, 2, 3]

        for orientation in orientations
            @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, aux_ll, aux_rr,
                                                       orientation,
                                                       equations)) == RealT
            @test typeof(@inferred max_abs_speed(u_ll, u_rr, aux_ll, aux_rr, orientation,
                                                 equations)) == RealT
        end
        @test eltype(@inferred flux(u, aux, normal_direction, equations)) == RealT
        @test eltype(@inferred flux_ec(u_ll, u_rr, aux_ll, aux_rr, normal_direction,
                                       equations)) ==
              RealT
        @test eltype(@inferred flux_tec(u_ll, u_rr, aux_ll, aux_rr, normal_direction,
                                        equations)) ==
              RealT
        @test eltype(@inferred flux_etec(u_ll, u_rr, aux_ll, aux_rr, normal_direction,
                                         equations)) ==
              RealT
        for i in 1:2
            @test eltype(@inferred flux_kennedy_gruber(u_ll, u_rr, aux_ll, aux_rr,
                                                       normal_direction,
                                                       equations)[i]) ==
                  RealT
            @test eltype(@inferred flux_lmars_combined(u_ll, u_rr, aux_ll, aux_rr,
                                                       normal_direction,
                                                       equations)[i]) ==
                  RealT
            @test eltype(@inferred flux_kennedy_gruber_souza_etal(u_ll, u_rr, aux_ll,
                                                                  aux_rr,
                                                                  normal_direction,
                                                                  equations)[i]) ==
                  RealT
            @test eltype(@inferred flux_lmars_es(u_ll, u_rr, aux_ll, aux_rr,
                                                 normal_direction, equations)[i]) ==
                  RealT
        end
        for volume_flux in ((flux_kennedy_gruber, flux_nonconservative_souza_etal),
                            (flux_ec, flux_nonconservative_isothermal))
            turbo = @inferred Trixi.cons2turbo(volume_flux, u..., aux..., equations)
            @test eltype(turbo) == RealT
            for i in 1:2
                @test eltype(@inferred Trixi.flux_turbo(volume_flux, turbo..., turbo...,
                                                        normal_direction..., equations)[i]) ==
                      RealT
            end
        end
        x_3d = SVector(RealT(EARTH_RADIUS), zero(RealT), zero(RealT))
        @test eltype(@inferred aux_field_gravity_earth(x_3d, equations)) == RealT
        @test eltype(@inferred aux_field_gravity_small_earth(x_3d, equations)) == RealT
        @test eltype(@inferred source_terms_coriolis(u, aux, x_3d, t, equations)) == RealT
        @test eltype(@inferred cons2prim(u, aux, equations)) == RealT
        @test eltype(@inferred prim2cons(u, aux, equations)) == RealT
        @test eltype(@inferred cons2entropy(u, aux, equations)) == RealT
        @test typeof(@inferred pressure(u, aux, equations)) == RealT
        @test typeof(@inferred entropy(cons, aux, equations)) == RealT
        @test typeof(@inferred energy_kinetic(cons, aux, equations)) == RealT
        @test typeof(@inferred energy_internal_kinetic(cons, aux, equations)) == RealT
        @test eltype(@inferred Trixi.max_abs_speeds(u, aux, equations)) == RealT
        @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, aux_ll, aux_rr,
                                                   normal_direction,
                                                   equations)) == RealT
        @test typeof(@inferred max_abs_speed(u_ll, u_rr, aux_ll, aux_rr, normal_direction,
                                             equations)) == RealT

        # Adapt the real type, e.g., to move a semidiscretization to a GPU
        equations_float64 = @inferred similar(equations, Float64)
        @test typeof(equations_float64.c_p) == Float64
        @test typeof(equations_float64.g) == Float64
        @test typeof(@inferred similar(equations_float64, RealT)) == typeof(equations)
    end
end

@testitem "Type stability: Compressible Euler Moist Potential Temperature 3D" setup=[Setup] tags=[:type_stable_tests] begin
    for RealT in (Float32, Float64)
        equations = @inferred CompressibleEulerMoistPotentialTemperatureEquations3D(c_p = RealT(1004),
                                                                                    c_v = RealT(717),
                                                                                    gravity = RealT(EARTH_GRAVITATIONAL_ACCELERATION))

        x = SVector(RealT(10_000), RealT(2_000), RealT(1_000))
        t = zero(RealT)
        u = u_ll = u_inner = SVector(RealT(1.1), RealT(-0.5), RealT(2.3), RealT(0.7),
                                     RealT(310), RealT(0.012), RealT(0.003),
                                     RealT(0.002))
        u_rr = SVector(RealT(1.05), RealT(0.4), RealT(-1.3), RealT(-0.2), RealT(300),
                       RealT(0.01), zero(RealT), RealT(0.004))
        aux = aux_ll = aux_rr = aux_inner = SVector(RealT(1700), RealT(0.6), RealT(0),
                                                    RealT(0.8))

        normal_direction = SVector(one(RealT), one(RealT), one(RealT))

        @test eltype(@inferred boundary_condition_slip_wall(u_inner, aux_inner,
                                                            normal_direction, x, t,
                                                            flux_lmars_combined,
                                                            equations)) == RealT
        @test eltype(@inferred flux(u, aux, normal_direction, equations)) == RealT
        for i in 1:2
            @test eltype(@inferred flux_lmars_combined(u_ll, u_rr, aux_ll, aux_rr,
                                                       normal_direction, equations)[i]) ==
                  RealT
            @test eltype(@inferred flux_kennedy_gruber_souza_etal(u_ll, u_rr, aux_ll,
                                                                  aux_rr,
                                                                  normal_direction,
                                                                  equations)[i]) ==
                  RealT
        end

        @test eltype(@inferred cons2prim(u, aux, equations)) == RealT
        @test eltype(@inferred prim2cons(u, aux, equations)) == RealT
        @test typeof(@inferred pressure(u, aux, equations)) == RealT
        @test typeof(@inferred TrixiAtmo.temperature(u, aux, equations)) == RealT
        @test typeof(@inferred TrixiAtmo.terminal_velocity_rain(TrixiAtmo.density_dry(u,
                                                                                      aux,
                                                                                      equations),
                                                                u[8],
                                                                equations.microphysics)) ==
              RealT
        @test typeof(@inferred TrixiAtmo.density_dry(u, aux, equations)) == RealT
        @test typeof(@inferred TrixiAtmo.saturation_mixing_ratio_tetens(RealT(280),
                                                                        RealT(90_000))) ==
              RealT
        @test eltype(@inferred Trixi.max_abs_speeds(u, aux, equations)) == RealT
        @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, aux_ll, aux_rr,
                                                   normal_direction, equations)) == RealT
        @test typeof(@inferred max_abs_speed(u_ll, u_rr, aux_ll, aux_rr,
                                             normal_direction, equations)) == RealT
        @test eltype(@inferred source_terms_microphysics(u, aux, x, t, equations)) == RealT
        @test eltype(@inferred source_terms_coriolis(u, aux, x, t, equations)) == RealT
        @test eltype(@inferred aux_field_gravity_earth(x, equations)) == RealT
        x_sphere = SVector(RealT(EARTH_RADIUS + 1_000), RealT(2_000), RealT(1_000))
        @test eltype(@inferred initial_condition_moist_baroclinic_instability(x_sphere, t,
                                                                              equations)) ==
              RealT
        turbo_ll = @inferred Trixi.cons2turbo(flux_kennedy_gruber_souza_etal, u_ll...,
                                              aux_ll..., equations)
        turbo_rr = Trixi.cons2turbo(flux_kennedy_gruber_souza_etal, u_rr..., aux_rr...,
                                    equations)
        @test eltype(turbo_ll) == RealT
        for flux_turbo_ in @inferred Trixi.flux_turbo(flux_kennedy_gruber_souza_etal,
                                                      turbo_ll..., turbo_rr...,
                                                      normal_direction..., equations)
            @test eltype(flux_turbo_) == RealT
        end

        # Adapt the real type, e.g., to move a semidiscretization to a GPU
        equations_float64 = @inferred similar(equations, Float64)
        @test typeof(equations_float64.c_p) == Float64
        @test typeof(equations_float64.microphysics.condensation_time) == Float64
        @test typeof(@inferred KesslerMicrophysics(latent_heat = RealT(2.5e6),
                                                   rho_surface = RealT(1.2),
                                                   condensation_time = one(RealT))) ==
              KesslerMicrophysics{RealT}
        @test typeof(@inferred similar(equations_float64, RealT)) == typeof(equations)
    end
end

@testitem "Type stability: Compressible Euler Moist Euler 2D" setup=[Setup] tags=[:type_stable_tests] begin
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

@testitem "Type stability: Compressible Rainy Euler" setup=[Setup] tags=[:type_stable_tests] begin
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

@testitem "Type stability: Shallow Water 3D" setup=[Setup] tags=[:type_stable_tests] begin
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

@testitem "Type stability: Compressible Euler Energy With Gravity 2D" setup=[Setup] tags=[:type_stable_tests] begin
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
        @test eltype(@inferred flux_nonconservative_isothermal(u_ll, u_rr,
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
        @test typeof(@inferred energy_internal_kinetic(cons, equations)) == RealT
        @test typeof(@inferred energy_internal(cons, equations)) == RealT
        @test eltype(@inferred Trixi.max_abs_speeds(u, equations)) == RealT
        @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, normal_direction,
                                                   equations)) == RealT
        @test typeof(@inferred max_abs_speed(u_ll, u_rr, normal_direction,
                                             equations)) == RealT
    end
end

@testitem "Type stability: Compressible Euler Energy With Gravity 3D" setup=[Setup] tags=[:type_stable_tests] begin
    for RealT in (Float32, Float64)
        equations = @inferred CompressibleEulerEnergyEquationsWithGravity3D(c_p = RealT(1004),
                                                                            c_v = RealT(717),
                                                                            gravity = RealT(9.81))

        x = SVector(zero(RealT))
        t = zero(RealT)
        u = u_ll = u_rr = u_inner = cons = SVector(one(RealT), one(RealT), one(RealT),
                                                   one(RealT), RealT(2))
        aux = aux_ll = aux_rr = SVector(one(RealT))

        normal_direction = SVector(one(RealT), one(RealT), one(RealT))

        orientations = [1, 2, 3]

        for orientation in orientations
            @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, aux_ll, aux_rr,
                                                       orientation,
                                                       equations)) == RealT
            @test typeof(@inferred max_abs_speed(u_ll, u_rr, aux_ll, aux_rr, orientation,
                                                 equations)) == RealT
        end
        @test eltype(@inferred flux(u, aux, normal_direction, equations)) == RealT
        @test eltype(@inferred flux_ranocha(u_ll, u_rr, aux_ll, aux_rr, normal_direction,
                                            equations)) ==
              RealT
        @test eltype(@inferred flux_kennedy_gruber(u_ll, u_rr, aux_ll, aux_rr,
                                                   normal_direction,
                                                   equations)) ==
              RealT
        @test eltype(@inferred flux_shima_etal(u_ll, u_rr, aux_ll, aux_rr, normal_direction,
                                               equations)) ==
              RealT
        @test eltype(@inferred flux_chandrashekar(u_ll, u_rr, aux_ll, aux_rr,
                                                  normal_direction,
                                                  equations)) == RealT
        @test eltype(@inferred flux_nonconservative_artiano_etal(u_ll, u_rr, aux_ll, aux_rr,
                                                                 normal_direction,
                                                                 equations)) ==
              RealT
        @test eltype(@inferred flux_nonconservative_waruszewski_etal(u_ll, u_rr, aux_ll,
                                                                     aux_rr,
                                                                     normal_direction,
                                                                     equations)) ==
              RealT
        @test eltype(@inferred flux_nonconservative_isothermal(u_ll, u_rr, aux_ll, aux_rr,
                                                               normal_direction,
                                                               equations)) ==
              RealT
        @test eltype(@inferred flux_nonconservative_souza_etal(u_ll, u_rr, aux_ll, aux_rr,
                                                               normal_direction,
                                                               equations)) ==
              RealT
        @test eltype(@inferred cons2prim(u, aux, equations)) == RealT
        @test eltype(@inferred prim2cons(u, aux, equations)) == RealT
        @test eltype(@inferred cons2entropy(u, aux, equations)) == RealT
        @test eltype(@inferred entropy2cons(cons2entropy(u, aux, equations), aux,
                                            equations)) ==
              RealT
        @test typeof(@inferred pressure(u, aux, equations)) == RealT
        @test typeof(@inferred entropy(cons, aux, equations)) == RealT
        @test typeof(@inferred energy_kinetic(cons, aux, equations)) == RealT
        @test typeof(@inferred energy_internal_kinetic(cons, aux, equations)) == RealT
        @test typeof(@inferred energy_internal(cons, aux, equations)) == RealT
        @test eltype(@inferred Trixi.max_abs_speeds(u, aux, equations)) == RealT
        @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, aux_ll, aux_rr,
                                                   normal_direction,
                                                   equations)) == RealT
        @test typeof(@inferred max_abs_speed(u_ll, u_rr, aux_ll, aux_rr, normal_direction,
                                             equations)) == RealT
    end
end

@testitem "Type stability: Compressible Euler Internal Energy With Gravity 2D" setup=[Setup] tags=[:type_stable_tests] begin
    for RealT in (Float32, Float64)
        equations = @inferred CompressibleEulerInternalEnergyEquationsWithGravity2D(c_p = RealT(1004),
                                                                                    c_v = RealT(717),
                                                                                    gravity = RealT(9.81))
        x = SVector(zero(RealT))
        t = zero(RealT)
        u = u_ll = u_rr = u_inner = cons = SVector(one(RealT), one(RealT), one(RealT),
                                                   RealT(2), one(RealT))
        surface_flux_functions = (flux_conservative_etec,
                                  flux_nonconservative_etec)
        normal_direction = SVector(one(RealT), one(RealT))

        @test eltype(@inferred flux_nonconservative_etec(u_ll, u_rr,
                                                         normal_direction,
                                                         equations)) ==
              RealT
        @test eltype(@inferred flux_conservative_etec(u_ll, u_rr,
                                                      normal_direction,
                                                      equations)) ==
              RealT
        @test eltype(@inferred flux_nonconservative_es(u_ll, u_rr,
                                                       normal_direction,
                                                       equations)) ==
              RealT
        @test eltype(@inferred flux_conservative_es(u_ll, u_rr,
                                                    normal_direction,
                                                    equations)) ==
              RealT
        @test eltype(@inferred boundary_condition_slip_wall(u_inner,
                                                            normal_direction,
                                                            x, t,
                                                            surface_flux_functions,
                                                            equations)) ==
              SVector{5, RealT}
        @test @inferred Trixi.have_nonconservative_terms(equations) === Trixi.True()
        @test eltype(varnames(cons2cons, equations)) == String
        @test eltype(varnames(cons2prim, equations)) == String
        @test eltype(@inferred cons2prim(u, equations)) == RealT
        @test eltype(@inferred cons2cons(u, equations)) == RealT
        @test eltype(@inferred prim2cons(u, equations)) == RealT
        @test eltype(@inferred cons2entropy(u, equations)) == RealT
        @test typeof(@inferred pressure(u, equations)) == RealT
        @test typeof(@inferred entropy(cons, equations)) == RealT
        @test eltype(@inferred Trixi.max_abs_speeds(u, equations)) == RealT
        @test typeof(@inferred max_abs_speed(u_ll, u_rr, normal_direction,
                                             equations)) == RealT
    end
end

@testitem "Type stability: Compressible Euler Vector Invariant With Gravity 3D" setup=[Setup] tags=[:type_stable_tests] begin
    for RealT in (Float32, Float64)
        equations = @inferred CompressibleEulerVectorInvariantEquationsWithGravity3D(c_p = RealT(1004),
                                                                                     c_v = RealT(717),
                                                                                     gravity = RealT(EARTH_GRAVITATIONAL_ACCELERATION))

        u = u_ll = u_rr = cons = SVector(one(RealT), one(RealT), one(RealT),
                                         one(RealT), one(RealT))
        aux = aux_ll = aux_rr = SVector(one(RealT))

        normal_direction = SVector(one(RealT), one(RealT), one(RealT))

        for i in 1:2
            @test eltype(@inferred flux_etec(u_ll, u_rr, aux_ll, aux_rr, normal_direction,
                                             equations)[i]) ==
                  RealT
            @test eltype(@inferred flux_lmars_es(u_ll, u_rr, aux_ll, aux_rr,
                                                 normal_direction, equations)[i]) ==
                  RealT
        end
        turbo = @inferred Trixi.cons2turbo(flux_etec, u..., aux..., equations)
        @test eltype(turbo) == RealT
        for i in 1:2
            @test eltype(@inferred Trixi.flux_turbo(flux_etec, turbo..., turbo...,
                                                    normal_direction..., equations)[i]) ==
                  RealT
        end
        @test eltype(@inferred cons2prim(u, aux, equations)) == RealT
        @test eltype(@inferred prim2cons(u, aux, equations)) == RealT
        @test eltype(@inferred cons2entropy(u, aux, equations)) == RealT
        @test typeof(@inferred pressure(u, aux, equations)) == RealT
        @test typeof(@inferred entropy(cons, aux, equations)) == RealT
        @test eltype(@inferred Trixi.max_abs_speeds(u, aux, equations)) == RealT
    end
end
