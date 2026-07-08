using Test
using BenchmarkTools

###############################################################################
# Set up unified equations

RealType = Float64
c_p = 1004.0
c_v = 717.0
g = EARTH_GRAVITATIONAL_ACCELERATION

parameters = Parameters{RealType}(;
                                  earth_gravitational_acceleration = g,
                                  c_dry_air_const_pressure = c_p,
                                  c_dry_air_const_volume = c_v)

td_single = IdealGas(; parameters)
td_potT = PotentialTemperature(td_single)

equations_uni = [
    CompressibleEulerAtmo(; n_dims = 1, n_vars_aux = 1,
                          parameters = parameters,
                          thermodynamic_state = td_single,
                          thermodynamic_equation = td_potT),
    CompressibleEulerAtmo(; n_dims = 2, n_vars_aux = 1,
                          parameters = parameters,
                          thermodynamic_state = td_single,
                          thermodynamic_equation = td_potT),
    CompressibleEulerAtmo(; n_dims = 3, n_vars_aux = 1,
                          parameters = parameters,
                          thermodynamic_state = td_single,
                          thermodynamic_equation = td_potT)]

###############################################################################
# Set up original equation

equations_orig = [
    CompressibleEulerPotentialTemperatureEquationsWithGravity1D(c_p = c_p, c_v = c_v,
                                                                gravity = g),
    CompressibleEulerPotentialTemperatureEquationsWithGravity2D(c_p = c_p, c_v = c_v,
                                                                gravity = g),
    CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = c_p, c_v = c_v,
                                                                gravity = g)]

###############################################################################
# Test

u_ll = SVector(1.1, -0.5, 2.34, 2.4, 330.0, 1000)
u_rr = SVector(1.2, -0.45, 1.89, 2.56, 321.0, 1100)

aux_ll = SVector(u_ll[6])
aux_rr = SVector(u_rr[6])

normal_directions = [SVector(0.5, -0.5, 0.2),
    SVector(-1.2, 0.3, 1.4),
    SVector(1.0, 0.0, 0.0),
    SVector(0.0, 1.0, 0.0),
    SVector(0.0, 0.0, 1.0)]

fluxes = [flux_ec, flux_tec, flux_etec, FluxLMARS(340),
    flux_nonconservative_waruszewski_etal,
    flux_nonconservative_artiano_etal,
    flux_nonconservative_souza_etal
]

bench = Array{BenchmarkTools.Trial, 3}(undef, 2, 2, length(fluxes))

for dim in 2:3
    println("Dimension $dim")
    u_orig_ll = SVector{dim + 3}(u_ll[1], u_ll[SVector{dim}(2:(dim + 1))]...,
                                 u_ll[5], u_ll[6])
    u_orig_rr = SVector{dim + 3}(u_rr[1], u_rr[SVector{dim}(2:(dim + 1))]...,
                                 u_rr[5], u_rr[6])

    u_uni_ll = SVector{dim + 2}(u_ll[SVector{dim}(2:(dim + 1))]..., u_ll[5],
                                u_ll[1])
    u_uni_rr = SVector{dim + 2}(u_rr[SVector{dim}(2:(dim + 1))]..., u_rr[5],
                                u_rr[1])

    for n in 1:(dim + 2)
        println("  Normal $n")
        normal_direction = SVector{dim}(normal_directions[n][SVector{dim}(1:dim)]...)

        for (iflux, _flux) in enumerate(fluxes)
            println("    Flux ", _flux)

            flux_orig = _flux(u_orig_ll, u_orig_rr, normal_direction, equations_orig[dim])
            flux_rearrange = SVector{dim + 2}(flux_orig[SVector{dim}(2:(dim + 1))]...,
                                              flux_orig[dim + 2],
                                              flux_orig[1])
            @test flux_rearrange ≈
                  _flux(u_uni_ll, u_uni_rr, aux_ll, aux_rr, normal_direction,
                        equations_uni[dim])

            if n == 1
                bench[1, dim - 1, iflux] = @benchmark ($_flux)($u_orig_ll, $u_orig_rr,
                                                               $normal_direction,
                                                               $(equations_orig[dim]))
                bench[2, dim - 1, iflux] = @benchmark ($_flux)($u_uni_ll, $u_uni_rr,
                                                               $aux_ll, $aux_rr,
                                                               $normal_direction,
                                                               $(equations_uni[dim]))
            end
        end
    end
end

for dim in 2:3
    println("\nDimension $dim")
    for (iflux, _flux) in enumerate(fluxes)
        println("    Flux ", _flux)
        display(bench[1, dim - 1, iflux])
        display(bench[2, dim - 1, iflux])
    end
end
