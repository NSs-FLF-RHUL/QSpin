
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
    ρ0 = 2e15, # Initial central density in g/cm^3, above 1e15 seems to be unstable
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
u0 = [EoS(Sim_Input.ρ0)/tovUnits.pressure_ref; 0.0];
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
problem = ODEProblem(tov!, u0, (r_beg, r_max); callback = cb)
sol = DE.solve(
    problem,
    OrdinaryDiffEqLowOrderRK.DP5();
    saveat = Dr,
    reltol = 1e-13,
    abstol = 1e-13,
)

ur = Array(sol)
r = sol.t * tovUnits.length_ref
Pr = ur[1, :] * tovUnits.pressure_ref
mr = ur[2, :] * tovUnits.mass_ref
ρr = EoS_inv.(Pr)
R_index = findfirst(x->x<0, Pr)
TOV_sol = (;
    r,
    Pr,
    mr,
    ρr,
    R = r[isnothing(R_index) ? end : (R_index - 1)],
    M = mr[isnothing(R_index) ? end : (R_index - 1)],
)

plot(r[1:(end)]*1e-5, ρr[1:(end)])


#plot(TOV_sol.r*tovUnits.length_ref*1e-5, TOV_sol.Pr*tovUnits.pressure_ref)
