
#=
This script demonstrates how to solve the Tolman-Oppenheimer-Volkoff (TOV) equation for neutron stars using the `QSpin` package.
See the documentation for the `QSpin.TOV` for more information.

The parameters, and equation of state, are chosen from https://github.com/vanessagraber/teaching_materials/blob/master/summerschool_CRAQ_2019/mass_radius_relations.ipynb.
=#

using QSpin
using QSpin.Parameters: ParameterType
using QSpin.TOV: TOV_Solve_dimensionless, TOV_ref_units, tov_eq_dimless!
using QSpin.TOV.EquationOfState: EoS_GCA2018, EoS_two_component_polytrope
using QSpin.PhysicalConstants: neutron_mass, mass_sun, hbar
using Plots, LaTeXStrings
import OrdinaryDiffEq as DE
using OrdinaryDiffEqLowOrderRK
tovUnits = TOV_ref_units()

Sim_Input = (
    ρ0 = 6.5e14 / tovUnits.rho_ref, # Initial central density in g/cm^3, above 1e15 seems to be unstable
    dr = 0.0001*1e5/tovUnits.length_ref, # Radial step in cm
    Dr = 0.005*1e5/tovUnits.length_ref, # Radial interval for recording values in cm
    r_end = 20e5/tovUnits.length_ref, # Maximum radius to solve up to in cm
    r_beg = 0.e5/tovUnits.length_ref,
);

# Calling EoS_GCA2018
EoS, EoS_inv = EoS_GCA2018();
# Polytropic EoS parameters
EoS_Param_Stiff = (
    K_crust = (3*π^2)^(2/3) * (hbar*1e7)^2 / (5*(neutron_mass*1e3)^(8/3)),
    γ_crust = 5.0/3.0,
    γ_core = 3.0,
    ρ_b = 3e14,
)
#EoS, EoS_inv = EoS_two_component_polytrope(EoS_Param_Stiff);

u0 = [EoS(Sim_Input.ρ0*tovUnits.rho_ref)/tovUnits.pressure_ref; 0.0];
tov! = tov_eq_dimless!(EoS_inv);
condition(u, t, integrator) = u[1] < 0
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)

problem = ODEProblem(tov!, u0, (0.0, Sim_Input.r_end); callback = cb)
sol = DE.solve(problem, OrdinaryDiffEqLowOrderRK.DP5(); reltol = 1e-12)
#sol = DE.solve(problem, DE.Tsit5(); reltol = 1e-12)

ur = Array(sol)
r = sol.t
Pr = ur[1, :]
mr = ur[2, :]
ρr = EoS_inv.(Pr*tovUnits.pressure_ref)
plot(r*tovUnits.length_ref*1e-5, ρr)


#plot(TOV_sol.r*tovUnits.length_ref*1e-5, TOV_sol.Pr*tovUnits.pressure_ref)
