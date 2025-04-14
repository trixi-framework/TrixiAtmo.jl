using Trixi
using OrdinaryDiffEq

struct HydrostaticSetup
	# Physical constants
	g::Float64       # gravity of earth
	c_p::Float64     # heat capacity for constant pressure (dry air)
	c_v::Float64     # heat capacity for constant volume (dry air)
	gamma::Float64   # heat capacity ratio (dry air)
	p_0::Float64     # atmospheric pressure
	T_0::Float64     # 
	u0::Float64      #
	Nf::Float64      # 
	z_B::Float64     #
	z_T::Float64     #
	function HydrostaticSetup(; g = 9.81, c_p = 1004.0, c_v = 717.0, gamma = c_p / c_v, p_0 = 100_000.0, T_0 = 250.0, u0 = 20.0, z_B = 15000.0, z_T = 30000.0)
		Nf = g / sqrt(c_p * T_0)
		new(g, c_p, c_v, gamma, p_0, T_0, u0, Nf, z_B, z_T)
	end
end

function (setup::HydrostaticSetup)(u, x, t, equations::CompressibleEulerEquations2D)
	@unpack g, c_p, c_v, gamma, p_0, T_0, z_B, z_T, Nf, u0 = setup

	rho, rho_v1, rho_v2, rho_e = u

	R = c_p - c_v

	v1 = rho_v1 / rho
	v2 = rho_v2 / rho
	p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2))
	rho_theta = (p / p_0)^(c_v / c_p) * p_0 / R
	theta = rho_theta / rho

	alfa = 0.1

	xr_B = 80000.0
	xr_T = 120000.0

	if x[2] <= z_B
		S_v = 0.0
	else
		S_v = -alfa * sinpi(0.5 * (x[2] - z_B) / (z_T - z_B))^2
	end
	if x[1] < xr_B
		S_h1 = 0.0
	else
		S_h1 = -alfa * sinpi(0.5 * (x[1] - xr_B) / (xr_T - xr_B))^2
	end

	if x[1] > -xr_B
		S_h2 = 0.0
	else
		S_h2 = -alfa * sinpi(0.5 * (x[1] + xr_B) / (-xr_T + xr_B))^2
	end

	exner = exp(-Nf^2 / g * x[2])

	theta_0 = T_0 / exner

	K = p_0 * (R / p_0)^gamma
	du2 = rho * (v1 - u0) * (S_v + S_h1 + S_h2)
	du3 = rho_v2 * (S_v + S_h1 + S_h2)
	du4 = rho * (theta - theta_0) * (S_v + S_h1 + S_h2) * K * gamma / (gamma - 1.0) * (rho_theta)^(gamma - 1.0) + du2 * v1 + du3 * v2

	return SVector(zero(eltype(u)), du2, du3 - g * rho, du4 - g * rho_v2)

end

function (setup::HydrostaticSetup)(x, t, equations::CompressibleEulerEquations2D)
	@unpack g, c_p, c_v, p_0, T_0, u0, Nf = setup

	# Exner pressure, solves hydrostatic equation for x[2]
	exner = exp(-Nf^2 / g * x[2])
	# pressure
	p_0 = 100_000.0  # reference pressure
	R = c_p - c_v    # gas constant (dry air)
	p = p_0 * exner^(c_p / R)

	# density
	rho = p / (R * T_0)
	v1 = u0
	v2 = 0.0

	return prim2cons(SVector(rho, v1, v2, p), equations)
end

linear_hydrostatic_setup = HydrostaticSetup()

equations = CompressibleEulerEquations2D(linear_hydrostatic_setup.gamma)

boundary_conditions = (x_neg = boundary_condition_periodic,
	x_pos = boundary_condition_periodic,
	y_neg = boundary_condition_slip_wall,
	y_pos = boundary_condition_slip_wall)

polydeg = 3
basis = LobattoLegendreBasis(polydeg)

surface_flux = FluxLMARS(340.0)
volume_flux = flux_kennedy_gruber
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

a = 10000.0
L = 240000.0
H = 30000.0
peak = 1.0
y_b = peak / (1 + (L / 2 / a)^2)
alfa = (H - y_b) * 0.5

f1(s) = SVector(-L / 2, y_b + alfa * (s + 1))
f2(s) = SVector(L / 2, y_b + alfa * (s + 1))
f3(s) = SVector((s + 1 - 1) * L / 2, peak / (1 + ((s + 1 - 1) * L / 2)^2 / a^2))
f4(s) = SVector((s + 1 - 1) * L / 2, H)

cells_per_dimension = (100, 50)
mesh = StructuredMesh(cells_per_dimension, (f1, f2, f3, f4), periodicity = (true, false))

semi = SemidiscretizationHyperbolic(mesh, equations, linear_hydrostatic_setup, solver, source_terms = linear_hydrostatic_setup,
	boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.
T = 5.0
tspan = (0.0, T * 3600.0)  # 1000 seconds final time
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000

analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback,
	analysis_callback,
	alive_callback)

###############################################################################
# run the simulation
sol = solve(ode,
	SSPRK43(thread = OrdinaryDiffEq.True()),
	maxiters = 1.0e7,
	dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
	save_everystep = false, callback = callbacks)

summary_callback()
