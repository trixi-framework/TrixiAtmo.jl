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


@inline n_aux_node_vars(::PerturbationEulerEquations2DAuxVars{}) = 2

varnames(::typeof(cons2cons), ::PerturbationEulerEquations2DAuxVars) = ("rho_prime", "rhov1", "rhov2", "rhoe_prime")
varnames(::typeof(cons2prim), ::PerturbationEulerEquations2DAuxVars) = ("rho_prime", "v1", "v2", "p_prime")
varnames(::typeof(cons2aux), ::PerturbationEulerEquations2DAuxVars) = ("rho_mean", "rhoE_mean")



@inline function Trixi.cons2prim(u, aux, equations::PerturbationEulerEquations2DAuxVars)  
    rho_prime, rho_v1, rho_v2, rhoe_prime = u
    rho_mean, rhoe_mean = aux 
    
    rho = rho_prime + rho_mean
    rho_e = rhoe_prime + rhoe_mean 

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho

    p = (equations.gamma - 1) * (rho_e - 0.5f0 *  rho * (v1^2 + v2^2))
    p_mean = (equations.gamma - 1) * (rhoe_mean - 0.5f0 *  rho_mean * (v1^2 + v2^2))
    p_prime = p - p_mean #p_prime as defined in giraldo under eq (6)
    return SVector(rho_prime, v1, v2, p_prime)
end

@inline function cons2primtotal(u, aux, equations::PerturbationEulerEquations2DAuxVars)
    rho_prime, rho_v1, rho_v2, rhoe_prime = u
    rho_mean, rhoe_mean = aux 
    
    rho = rho_prime + rho_mean
    rho_e = rhoe_mean + rhoe_prime

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho

    p = (equations.gamma - 1) * (rho_e - 0.5f0 *  rho * (v1^2 + v2^2))   
    return SVector(rho, v1, v2, p)
end 

@inline function cons2constotal(u, aux, equations::PerturbationEulerEquations2DAuxVars)
    rho_prime, rho_v1, rho_v2, rhoe_prime = u
    rho_mean, rhoe_mean = aux 

    rho = rho_prime + rho_mean
    rho_e = rhoe_prime + rhoe_mean
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


@inline function Trixi.flux(u, aux, orientation::Integer, equations::PerturbationEulerEquations2DAuxVars)
    # perturbated and background values 
    rho_prime, v1, v2, p_prime = cons2prim(u, aux, equations)      
    rho_mean, rhoe_mean = aux

    # total 
    rho = rho_prime + rho_mean
    rho_e = u[4] + rhoe_mean
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
    rho_mean, rhoe_mean = aux 

    # total 
    rho = rho_prime + rho_mean
    rho_e = u[4] + rhoe_mean
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

end