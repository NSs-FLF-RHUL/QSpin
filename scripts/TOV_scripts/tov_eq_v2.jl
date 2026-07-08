
#=
This script demonstrates how to solve the Tolman-Oppenheimer-Volkoff (TOV) equation for neutron stars using the `QSpin` package.
See the documentation for the `QSpin.TOV` for more information.

The parameters, and equation of state, are chosen from https://github.com/vanessagraber/teaching_materials/blob/master/summerschool_CRAQ_2019/mass_radius_relations.ipynb.
=#

using QSpin
using QSpin.Parameters: ParameterType
using QSpin.TOV: TOV_Solve_dimensionless, TOV_ref_units, tov_eq!
using QSpin.TOV.EquationOfState: EoS_GCA2018, EoS_two_component_polytrope
using QSpin.PhysicalConstants: neutron_mass, mass_sun, hbar
using Plots, LaTeXStrings
import CommonSolve as DE
using OrdinaryDiffEqLowOrderRK

Sim_Input = (
    ρ0 = 5e14, # Initial central density in g/cm^3, above 1e15 seems to be unstable
    dr = 0.0001*1e5, # Radial step in cm
    Dr = 0.005*1e5, # Radial interval for recording values in cm
    r_beg = 0.e5,
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
tovUnits = TOV_ref_units()
u0 = [EoS(Sim_Input.ρ0*tovUnits.rho_ref)/tovUnits.pressure_ref; 0.0];
# Building the dimensionless TOV equation function using the provided inverse EoS function
tov! = tov_eq!(EoS_inv);
tovUnits = TOV_ref_units();
dt = Sim_Input.dr / tovUnits.length_ref;
Dr = Sim_Input.Dr / tovUnits.length_ref;
r_beg = Sim_Input.r_beg / tovUnits.length_ref;
r_max = 20*1e3*1e2 / tovUnits.length_ref;
# Callback setup to terminate the integration when the pressure drops below zero
condition(u, t, integrator) = u[1] < 0
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)

sol_tov =
    evolve(tov!, u0, 0.0, r_max; alg, reltol, callback = cb, dt, saveat, solver_options...)
ur = Array(sol_tov)
r = sol_tov.t * tovUnits.length_ref
Pr = ur[1, :] * tovUnits.pressure_ref
mr = ur[2, :] * tovUnits.mass_ref
plot(r*tovUnits.length_ref*1e-5, ρr)


#plot(TOV_sol.r*tovUnits.length_ref*1e-5, TOV_sol.Pr*tovUnits.pressure_ref)
