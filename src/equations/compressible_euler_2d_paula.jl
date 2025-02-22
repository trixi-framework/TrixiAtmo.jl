using Trixi 
import Trixi: varnames

#function getpottemp(u, equations::CompressibleEulerEquations2D)
#    @unpack g, c_p, c_v = interia_wave_setup

#    p_0 = 100_000.0  # reference pressure
#    R = c_p - c_v    # gas constant (dry air)
    
#    rho, rho_v1, rho_v2, rho_e = u
#    E = rho_e / rho 
#    v1 = rho_v1 / rho 
#    v2 = rho_v2 / rho 

#    T = (E - 0.5 * (v1^2 + v2^2)) / c_v
#    p = R * T * rho 
#    exner = (p / p_0) ^ (R/c_p)
    
#    pot = T / exner 
#    return pot
#end

@inline function cons2pot(u, equations::CompressibleEulerEquations2D)
    rho, rho_v1, rho_v2, rho_e = u

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    E = rho_e / rho 


    p_0 = 100_000.0  # reference pressure
    c_p = 1004.0
    c_v = 717.0
    R = c_p - c_v    # gas constant (dry air)

    T = (E - 0.5 * (v1^2 + v2^2)) / c_v
    p = R * T * rho 
    exner = (p / p_0) ^ (R/c_p)
    

    potential_temperature = p_0 /(R * rho) * exner ^ (c_v /R)
    #potential_temperature = T / exner 

    #test = potential_temperature / rho

    #potential_temperature_int = 300.0 #constant of integration 
    #bvfrequency = 0.01 #Brunt-Väisälä frequency
    #g = 9.81
    
    #potential_temperature_mean = potential_temperature_int * exp(bvfrequency^2 / g * x[2])

    #potential_temperature_perturbation  = potential_temperature - potential_temperature_mean

    return SVector(rho, v1, v2, potential_temperature) 
end

varnames(::typeof(cons2pot), ::CompressibleEulerEquations2D) = ("rho", "v1", "v2", "pot_temp")

#Rayleigh dumping sponge source term 
@inline function source_terms_rayleigh_sponge(u, x, t, equations:: CompressibleEulerEquations2D)
    rho, rho_v1, rho_v2, rho_e = u

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    z = x[2]

    # relaxed background velocity
    vr1, vr2 = (10.0, 0.0)
    # damping threshold
    z_s = 15000.0
    # boundary top
    z_top = 21000.0
    # positive even power with default value 2
    gamma = 2.0
    # relaxation coefficient > 0
    alpha = 0.5

    tau_s = zero(eltype(u))
    if z > z_s
        tau_s = alpha * sin(0.5 * (z - z_s) * inv(z_top - z_s))^(gamma)
    end

    return SVector(zero(eltype(u)),
                   -tau_s * rho * (v1 - vr1),
                   -tau_s * rho * (v2 - vr2),
                   zero(eltype(u)))
end

@inline function source_terms_rayleigh_sponge_left(u, x, t, equations:: CompressibleEulerEquations2D)
    rho, rho_v1, rho_v2, rho_e = u

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    y = x[1]

    # relaxed background velocity
    vr1, vr2 = (10.0, 0.0)
    # damping threshold
    y_s = -20000.0
    # boundary left
    y_top = -25000.0
    # positive even power with default value 2
    gamma = 2.0
    # relaxation coefficient > 0
    alpha = 0.5

    tau_s = zero(eltype(u))
    if y < y_s
        tau_s = alpha * sin(0.5 * (y - y_s) * inv(y_top - y_s))^(gamma)
    end

    return SVector(zero(eltype(u)),
                   -tau_s * rho * (v1 - vr1),
                   -tau_s * rho * (v2 - vr2),
                   zero(eltype(u)))
end


@inline function source_terms_rayleigh_sponge_right(u, x, t, equations:: CompressibleEulerEquations2D)
    rho, rho_v1, rho_v2, rho_e = u

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    y = x[1]

    # relaxed background velocity
    vr1, vr2 = (10.0, 0.0)
    # damping threshold
    y_s = 20000.0
    # boundary right
    y_top = 25000.0
    # positive even power with default value 2
    gamma = 2.0
    # relaxation coefficient > 0
    alpha = 0.5

    tau_s = zero(eltype(u))
    if y > y_s
        tau_s = alpha * sin(0.5 * (y - y_s) * inv(y_top - y_s))^(gamma)
    end

    return SVector(zero(eltype(u)),
                   -tau_s * rho * (v1 - vr1),
                   -tau_s * rho * (v2 - vr2),
                   zero(eltype(u)))
end



@inline function source_terms_gravity(u, equations:: CompressibleEulerEquations2D)
    g = 9.81
    rho, rho_v1, rho_v2, rho_e = u
    return SVector(zero(eltype(u)), zero(eltype(u)), -g * rho, -g * rho_v2)
end