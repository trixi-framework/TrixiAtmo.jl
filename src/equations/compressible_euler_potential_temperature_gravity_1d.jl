using Trixi
using Trixi: ln_mean, stolarsky_mean, AbstractCompressibleEulerEquations, inv_ln_mean
import Trixi: varnames, cons2cons, cons2prim, cons2entropy, entropy, energy_total, energy_kinetic, pressure,
	flux_ec, initial_condition_density_wave, max_abs_speeds

@muladd begin

#! format: noindent
struct CompressibleEulerPotentialTemperatureEquationsWithGravity1D{RealT <: Real} <:
	   AbstractCompressibleEulerEquations{1, 4}
	p_0::RealT
	c_p::RealT
	c_v::RealT
	g::RealT
	R::RealT
	gamma::RealT
	a::RealT
	inv_gamma_minus_one::RealT
	K::RealT
	stolarsky_factor::RealT
end

function CompressibleEulerPotentialTemperatureEquationsWithGravity1D(; g = 9.81, RealT = Float64)
	p_0 = 100_000.0
	c_p = 1004.0
	c_v = 717.0
	R = c_p - c_v
	gamma = c_p / c_v
	a = 340.0
	inv_gamma_minus_one = inv(gamma - 1)
	K = p_0 * (R / p_0)^gamma
	stolarsky_factor = (gamma - 1.0)/gamma
	return CompressibleEulerPotentialTemperatureEquationsWithGravity1D{RealT}(p_0, c_p, c_v, g, R, gamma, a, inv_gamma_minus_one, K, stolarsky_factor)
end

function varnames(::typeof(cons2cons), ::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
	("rho", "rho_v1", "rho_theta", "phi")
end

varnames(::typeof(cons2prim), ::CompressibleEulerPotentialTemperatureEquationsWithGravity1D) = ("rho", "v1", "p1", "phi")

Trixi.have_nonconservative_terms(::CompressibleEulerPotentialTemperatureEquationsWithGravity1D) = Trixi.True()

@inline function flux(u, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
	rho, rho_v1, rho_theta = u
	v1 = rho_v1 / rho
	p = equations.p_0 * (equations.R * rho_theta / equations.p_0)^equations.gamma
	p = equations.K * exp(log(rho_theta^equations.gamma))

	f1 = rho_v1
	f2 = rho_v1 * v1 + p
	f3 = rho_theta * v1

	return SVector(f1, f2, f3, 0)
end

@inline function (flux_lmars::FluxLMARS)(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
	a = flux_lmars.speed_of_sound
	rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
	rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)

	rho = 0.5f0 * (rho_ll + rho_rr)

	p_interface = 0.5f0 * (p_ll + p_rr) - 0.5f0 * a * rho * (v1_rr - v1_ll)
	v_interface = 0.5f0 * (v1_ll + v1_rr) - 1 / (2 * a * rho) * (p_rr - p_ll)

	if (v_interface > 0)
		f1, f2, f3 = u_ll * v_interface
	else
		f1, f2, f3 = u_rr * v_interface
	end

	return SVector(f1,
		f2 + p_interface,
		f3, zero(eltype(u)))
end

"""
	flux_tec(u_ll, u_rr, orientation_or_normal_direction,
						equations::CompressibleEulerEquationsPotentialTemperature1D)

	Total energy conservative two-point flux by
-  Artiano et al. (2025), pre-print
   Structure-Preserving High-Order Methods for the Compressible Euler Equations 
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_tec(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
	# Unpack left and right state
	rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
	rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)
	_, _, rho_theta_ll = u_ll
	_, _, rho_theta_rr = u_rr

	# Compute the necessary mean values
	rho_mean = ln_mean(rho_ll, rho_rr)

	gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * v1_avg
	f2 = f1 * v1_avg + p_avg
	f3 = gammamean * v1_avg
	return SVector(f1, f2, f3, zero(eltype(u)))
end

"""
	flux_ec(u_ll, u_rr, orientation_or_normal_direction,
						equations::CompressibleEulerEquationsPotentialTemperature1D)

Entropy conservative two-point flux by
-  Artiano et al. (2025), pre-print
   Structure-Preserving High-Order Methods for the Compressible Euler Equations 
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_ec(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
	# Unpack left and right state
	rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
	rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)
	_, _, rho_theta_ll = u_ll
	_, _, rho_theta_rr = u_rr
	# Compute the necessary mean values
	rho_mean = ln_mean(rho_ll, rho_rr)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * v1_avg
	f2 = f1 * v1_avg + p_avg
	f3 = inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr) * f1
	return SVector(f1, f2, f3, zero(eltype(u)))
end

"""
	flux_etec(u_ll, u_rr, orientation_or_normal_direction,
						equations::CompressibleEulerEquationsPotentialTemperature1D)

Entropy and total energy conservative two-point flux by
-  Artiano et al. (2025), pre-print
   Structure-Preserving High-Order Methods for the Compressible Euler Equations 
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_etec(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
	# Unpack left and right state
	rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
	rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)
	_, _, rho_theta_ll = u_ll
	_, _, rho_theta_rr = u_rr
	# Compute the necessary mean values
	gammamean = stolarsky_mean_opt(rho_theta_ll, rho_theta_rr, equations.gamma, p_ll, p_rr)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	# Calculate fluxes depending on normal_direction
	f3 = gammamean * v1_avg
	f1 = f3 * ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
	f2 = f1 * v1_avg + p_avg
	return SVector(f1, f2, f3, zero(eltype(u)))
end

@inline function prim2cons(prim, equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
	rho, v1, p, phi = prim
	rho_v1 = rho * v1
	rho_theta = equations.p_0 / equations.R * exp(1/equations.gamma * log(p/equations.p_0))
	return SVector(rho, rho_v1, rho_theta, phi)
end

@inline function cons2prim(u, equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
	rho, rho_v1, rho_theta = u
	v1 = rho_v1 / rho
	p = equations.K * exp(equations.gamma * log(rho_theta))
	return SVector(rho, v1, p, phi)
end

@inline function cons2cons(u, equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
	return u
end

@inline function cons2entropy_rhoe(u, equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
	rho, rho_v1, rho_theta = u

	w1 = -0.5f0 * rho_v1^2 / (rho)^2
	w2 = rho_v1 / rho
	w3 = equations.gamma * equations.inv_gamma_minus_one * exp((equations.gamma-1)*log(rho_theta))

	return SVector(w1, w2, w3, zero(eltype(u)))
end

@inline function cons2entropy(u, equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
	rho, rho_v1, rho_theta = u

	w1 = log(equations.K * (rho_theta / rho)^equations.gamma) - equations.gamma
	w2 = 0.0
	w3 = rho / rho_theta * equations.gamma

	return SVector(w1, w2, w3, zero(eltype(u)))
end

@inline function energy_total(cons, equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
	# Mathematical entropy
	p = equations.p_0 * (equations.R * cons[3] / equations.p_0)^equations.gamma

	U = p / (equations.gamma - 1) + 1 / 2 * (cons[2]^2) / (cons[1]) + cons[1] * cons[4]

	return U
end

# Default entropy is the mathematical entropy
@inline function entropy(cons, equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
	entropy_phys(cons, equations)
end

@inline function entropy_phys(cons, equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)

	p = equations.K * (cons[3])^equations.gamma
	s = log(p) - equations.gamma * log(cons[1])
	S = -s * cons[1] / (equations.gamma - 1.0)
	return S
end

@inline function energy_kinetic(cons, equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
	return 0.5f0 * (cons[2]^2) / (cons[1])
end

@inline function max_abs_speeds(u, equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
	rho, v1, p = cons2prim(u, equations)
	c = sqrt(equations.gamma * p / rho)

	return (abs(v1) + c,)
end

@inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
	rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
	rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)

	# Calculate primitive variables and speed of sound
	v_mag_ll = abs(v1_ll)
	c_ll = sqrt(equations.gamma * p_ll / rho_ll)
	v_mag_rr = abs(v1_rr)
	c_rr = sqrt(equations.gamma * p_rr / rho_rr)

	λ_max = max(v_mag_ll, v_mag_rr) + max(c_ll, c_rr)
end

@inline function pressure(cons, equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
	_, _, p = cons2prim(cons, equations)
	return p

end

end # @muladd
