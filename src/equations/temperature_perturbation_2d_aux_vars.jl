@muladd begin

struct PerturbationEulerEquations2DAuxVars{RealT <: Real} <:
       AbstractVariableCoefficientEquations{2,4} 
    gamma::RealT               # ratio of specific heats
    inv_gamma_minus_one::RealT # = inv(gamma - 1); can be used to write slow divisions as fast multiplications

    function PerturbationEulerEquations2DAuxVars(gamma)
        γ, inv_gamma_minus_one = promote(gamma, inv(gamma - 1))
        return new{typeof(γ)}(γ, inv_gamma_minus_one)
    end
end 


@inline n_aux_node_vars(::PerturbationEulerEquations2DAuxVars{}) = 4

varnames(::typeof(cons2cons), ::PerturbationEulerEquations2DAuxVars) = ("rho_prime", "rho_v1", "rho_v2", "rhoe_prime")
varnames(::typeof(cons2prim), ::PerturbationEulerEquations2DAuxVars) = ("rho_prime", "v1", "v2", "p_prime")
varnames(::typeof(cons2aux), ::PerturbationEulerEquations2DAuxVars) = ("rho_mean", "v1_mean", "v2_mean", "e_mean")



@inline function Trixi.cons2prim(u, aux, equations::PerturbationEulerEquations2DAuxVars)  
    rho_prime, rho_v1, rho_v2, rhoe_prime = u
    rho_mean, v1_mean, v2_mean, e_mean = aux
    
    rho = rho_prime + rho_mean
    
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho

    e_prime = rhoe_prime / rho
    e = e_prime + e_mean

    p = (equations.gamma - 1) * rho * (e - 0.5f0  * (v1^2 + v2^2))
    p_mean = (equations.gamma - 1) * rho_mean * (e_mean - 0.5f0  * (v1_mean^2 + v2_mean^2))
    p_prime = p - p_mean #p_prime as defined in giraldo under eq (6)

    return SVector(rho_prime, v1, v2, p_prime)
end

@inline function cons2primtotal(u, aux, equations::PerturbationEulerEquations2DAuxVars)
    rho_prime, rho_v1, rho_v2, rhoe_prime = u
    rho_mean, v1_mean, v2_mean, e_mean = aux
    
    rho = rho_prime + rho_mean
    
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho

    e_prime = rhoe_prime / rho
    e = e_prime + e_mean

    p = (equations.gamma - 1) * rho * (e - 0.5f0  * (v1^2 + v2^2))
    return SVector(rho, v1, v2, p)
end 

@inline function cons2constotal(u, aux, equations::PerturbationEulerEquations2DAuxVars)
    rho_prime, rho_v1, rho_v2, rhoe_prime = u
    rho_mean, v1_mean, v2_mean, e_mean = aux

    rho = rho_prime + rho_mean

    rho_e = rhoe_prime + (rho * e_mean)
    return SVector(rho, rho_v1, rho_v2, rho_e)
end 


# TODO -? 
@inline function cons2entropy(u, aux, equations::PerturbationEulerEquations2DAuxVars)
    rho, rho_v1, rho_v2, rho_e = cons2constotal(u, aux, equations)

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v_square = v1^2 + v2^2
    p = (equations.gamma - 1) * (rho_e - 0.5f0 * rho * v_square)
    s = log(p) - equations.gamma * log(rho)
    rho_p = rho / p

    w1 = (equations.gamma - s) * equations.inv_gamma_minus_one -
         0.5f0 * rho_p * v_square
    w2 = rho_p * v1
    w3 = rho_p * v2
    w4 = -rho_p

    return SVector(w1, w2, w3, w4)
end

@inline function cons2temppert(u, aux, equations::PerturbationEulerEquations2DAuxVars)
    rho, rho_v1, rho_v2, rho_e = cons2constotal(u, aux, equations)
    e = rho_e / rho #    e = c_v * theta * exner + 0.5 * (v1^2 + v2^2)
    v1 = rho_v1 / rho 
    v2 = rho_v2 / rho 

        # constants 
    g = 9.81 
    c_p = 1004.0 
    c_v = 717.0

    theta_c = 0.01 
    h_c = 10_000.0
    a_c = 5_000.0
    x_c = 100_000.0 
    theta_0 = 300.0 # constant of integration 
    bvfrequency = 0.01 # Brunt-Väisälä frequency
    p_0 = 100_000.0  # reference pressure
    R = c_p - c_v

    rho_mean, v1_mean, v2_mean, e_mean = aux
    eint_mean = e_mean - 0.5 * (v1_mean^2 + v2_mean^2)
    p_mean = (equations.gamma - 1) * rho_mean * eint_mean
    exner_mean = (p_mean / p_0)^(R / c_p)
    theta_mean = eint_mean / (c_v * exner_mean)

    eint = e - 0.5 * (v1^2 + v2^2)
    p = (equations.gamma - 1) * rho * eint
    exner = (p / p_0)^(R / c_p)
    theta = eint / (c_v * exner)
    theta_prime = theta - theta_mean
    return SVector(rho, v1, v2, theta_prime)
end 
varnames(::typeof(cons2temppert), ::PerturbationEulerEquations2DAuxVars) = ("rho", "v1", "v2", "theta_prime")

@inline function cons2all(u, aux, equations::PerturbationEulerEquations2DAuxVars) #for solution variables like cons2temppert and background state 
    rho, rho_v1, rho_v2, rho_e = cons2constotal(u, aux, equations)
    rho_mean, v1_mean, v2_mean, e_mean = aux

    e = rho_e / rho #    e = c_v * theta * exner + 0.5 * (v1^2 + v2^2)
    v1 = rho_v1 / rho 
    v2 = rho_v2 / rho 

        # constants 
    g = 9.81 
    c_p = 1004.0 
    c_v = 717.0

    theta_c = 0.01 
    h_c = 10_000.0
    a_c = 5_000.0
    x_c = 100_000.0 
    theta_0 = 300.0 # constant of integration 
    bvfrequency = 0.01 # Brunt-Väisälä frequency
    p_0 = 100_000.0  # reference pressure
    R = c_p - c_v


    # exner pressure 

    T = (e - 0.5 * (v1^2 + v2^2)) / c_v
    p = R * T * rho 
    exner = (p / p_0) ^ (R/c_p)

    theta = (e - 0.5 * (v1^2 + v2^2)) / (c_v * exner)

    theta_mean = p_0 / (R * rho_mean) * exner ^ (c_v / R)
    theta_prime = theta - theta_mean
    return SVector(rho, v1, v2, theta_prime, rho_mean, v1_mean, v2_mean, e_mean, theta_mean, theta)
end 

varnames(::typeof(cons2all), ::PerturbationEulerEquations2DAuxVars) = ("rho", "v1", "v2", "theta_prime", "rho_mean", "v1_mean", "v2_mean", "e_mean", "theta_mean", "theta_total")


@inline function Trixi.flux(u, aux, orientation::Integer, equations::PerturbationEulerEquations2DAuxVars)
    # perturbated and background values 
    rho_prime, v1, v2, p_prime = cons2prim(u, aux, equations)      
    rho_mean, v1_mean, v2_mean, e_mean = aux
    # total 
    rho = rho_prime + rho_mean
    rho_e = u[4] + (rho * e_mean)
    p = (equations.gamma - 1) * (rho_e - 0.5f0 * rho * (v1^2 + v2^2))
    
    if orientation == 1
        f1 = rho * v1
        f2 = rho * v1 * v1 + p_prime
        f3 = rho * v1 * v2 
        f4 = (rho_e + p) * v1
    else  # orientation == 2
        f1 = rho * v2
        f2 = rho * v2 * v1 
        f3 = rho * v2 * v2 + p_prime 
        f4 = (rho_e + p) * v2
    end 
    return SVector(f1, f2, f3, f4)
end 

@inline function Trixi.flux(u, aux, normal_direction::AbstractVector,
       equations::PerturbationEulerEquations2DAuxVars)
    # perturbated and background values 
    rho_prime, v1, v2, p_prime = cons2prim(u, aux, equations)      
    rho_mean, v1_mean, v2_mean, e_mean = aux
    # total 
    rho = rho_prime + rho_mean
    rho_e = u[4] + (rho * e_mean)

    p = (equations.gamma - 1) * (rho_e - 0.5f0 * rho * (v1^2 + v2^2))
    
    # normal vectors
    v_normal = v1 * normal_direction[1] + v2 * normal_direction[2]
    rho_v_normal = rho  * v_normal

    #flux
    f1 = rho_v_normal
    f2 = rho_v_normal * v1 + p_prime * normal_direction[1]
    f3 = rho_v_normal * v2 + p_prime * normal_direction[2]
    f4 = (rho_e + p) * v_normal
    return SVector(f1, f2, f3, f4)
end

@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, aux_ll, aux_rr, orientation::Integer,
                                     equations::PerturbationEulerEquations2DAuxVars)
    rho_ll, v1_ll, v2_ll, p_ll = cons2primtotal(u_ll, aux_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2primtotal(u_rr, aux_rr, equations)   

    if orientation == 1
        v_ll = v1_ll
        v_rr = v1_rr
    else # orientation == 2
        v_ll = v2_ll
        v_rr = v2_rr
    end

    # Calculate sound speeds with total values for p and rho 
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)
         
    return max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr)
end 

@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, aux_ll, aux_rr, normal_direction::AbstractVector,
                                     equations::PerturbationEulerEquations2DAuxVars)
    rho_ll, v1_ll, v2_ll, p_ll = cons2primtotal(u_ll, aux_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2primtotal(u_rr, aux_rr, equations)   
    
    # Calculate normal velocities and sound speed
    # left
    v_ll = (v1_ll * normal_direction[1]
            +
            v2_ll * normal_direction[2])

    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    # right
    v_rr = (v1_rr * normal_direction[1]
            +
            v2_rr * normal_direction[2])

    c_rr = sqrt(equations.gamma * p_rr / rho_rr)    

    return max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr) * norm(normal_direction)
end 

# Maximum wave speeds in each direction for CFL calculation
@inline function Trixi.max_abs_speeds(u, aux, equations::PerturbationEulerEquations2DAuxVars)
    rho, v1, v2, p = cons2primtotal(u, aux, equations)
    c = sqrt(equations.gamma * p / rho)

    return abs(v1) + c, abs(v2) + c
end


@inline function max_abs_speed(u_ll, u_rr, aux_ll, aux_rr,
                               normal_direction::AbstractVector,
                               equations::PerturbationEulerEquations2DAuxVars)

    rho_ll, v1_ll, v2_ll, p_ll = cons2primtotal(u_ll, aux_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2primtotal(u_rr, aux_rr, equations)

    # Calculate normal velocities and sound speeds
    # left
    v_ll = (v1_ll * normal_direction[1]
            +
            v2_ll * normal_direction[2])

    if equations.gamma * p_ll / rho_ll <= 0
        @show p_ll, rho_ll
    elseif equations.gamma * p_rr / rho_rr <= 0
        @show p_rr, rho_rr
    end

    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    # right
    v_rr = (v1_rr * normal_direction[1]
            +
            v2_rr * normal_direction[2])
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    norm_ = norm(normal_direction)
    return max(abs(v_ll) + c_ll * norm_,
               abs(v_rr) + c_rr * norm_)
end


@inline function boundary_condition_slip_wall_toro_aux(u_inner, aux_inner, normal_direction::AbstractVector, 
                                                x, t, surface_flux_function, 
                                                equations::PerturbationEulerEquations2DAuxVars)
    norm_ = norm(normal_direction)
    # Normalize the vector without using `normalize` since we need to multiply by the `norm_` later
    normal = normal_direction / norm_

    rho, v1, v2, p = cons2primtotal(u_inner, aux_inner, equations)
    rho_mean, v1_mean, v2_mean, e_mean = aux_inner
    p_mean = (equations.gamma - 1) * rho_mean * (e_mean - 0.5f0  * (v1_mean^2 + v2_mean^2))


    # rotate velocities
    v_normal  = v1 * normal[1] + v2 * normal[2] # wall normal
    v_tangent = -v1 * normal[2] + v2 * normal[1]# wall tangent

    # Get the solution of the pressure Riemann problem
    # See Section 6.3.3 of
    # Eleuterio F. Toro (2009)
    # Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction
    # [DOI: 10.1007/b79761](https://doi.org/10.1007/b79761)
    if v_normal <= 0
        sound_speed = sqrt(equations.gamma * p / rho) # local sound speed
        p_star = p *
                 (1 + 0.5f0 * (equations.gamma - 1) * v_normal / sound_speed)^(2 *
                                                                               equations.gamma *
                                                                               equations.inv_gamma_minus_one)
    else # v_normal > 0
        A = 2 / ((equations.gamma + 1) * rho)
        B = p * (equations.gamma - 1) / (equations.gamma + 1)
        p_star = p +
                 0.5f0 * v_normal / A *
                 (v_normal + sqrt(v_normal^2 + 4 * A * (p + B)))
    end
    p_star_prime = p_star - p_mean
    # For the slip wall we directly set the flux as the normal velocity is zero
    return SVector(0,
                   p_star_prime * normal[1],
                   p_star_prime * normal[2],
                   0) * norm_
end

@inline function boundary_condition_slip_wall_aux(u_inner, aux_inner,
                                              normal_direction::AbstractVector,
                                              x, t,
                                              surface_flux_function,
                                              equations::PerturbationEulerEquations2DAuxVars)
    # normalize the outward pointing direction
    normal = normal_direction / norm(normal_direction)
    
    # compute the normal velocity
    u_normal = normal[1] * u_inner[2] + normal[2] * u_inner[3]

    # create the "external" boundary solution state
    u_boundary = SVector(u_inner[1],
                         u_inner[2] - 2 * u_normal * normal[1],
                         u_inner[3] - 2 * u_normal * normal[2],
                         u_inner[4])

    # compute the normal background velocity
    aux_normal = normal[1] * aux_inner[2] + normal[2] * aux_inner[3]

    # create the "external" boundary solution for background state
    aux_boundary = SVector(aux_inner[1],
                         aux_inner[2] - 2 * aux_normal * normal[1],
                         aux_inner[3] - 2 * aux_normal * normal[2],
                         aux_inner[4])

    # calculate the boundary flux
    flux = surface_flux_function(u_inner, u_boundary, aux_inner, aux_boundary, normal_direction, equations)
    return flux
end


    

end