
#=
This script demonstrates how to solve the Tolman-Oppenheimer-Volkoff (TOV) equation for neutron stars using the `QSpin` package.
See the documentation for the `QSpin.TOV` for more information.

The parameters, and equation of state, are chosen from https://github.com/vanessagraber/teaching_materials/blob/master/summerschool_CRAQ_2019/mass_radius_relations.ipynb.
=#

using QSpin
using QSpin.Parameters: ParameterType
using QSpin.TOV: TOV_Solve_dimensionless, TOV_ref_units
using QSpin.TOV.EquationOfState: EoS_GCA2018
using QSpin.PhysicalConstants: neutron_mass, mass_sun
using Plots, LaTeXStrings
import OrdinaryDiffEq as DE
using OrdinaryDiffEqLowOrderRK
tovUnits = TOV_ref_units()

Sim_Input = (
    ρ0 = 0.08*(1e39)*neutron_mass*1e3 / tovUnits.rho_ref, # Initial central density in g/cm^3, above 1e15 seems to be unstable
    dr = 0.002*1e5/tovUnits.length_ref, # Radial step in cm
    Dr = 0.01*1e5/tovUnits.length_ref, # Radial interval for recording values in cm
    r_end = 20e5/tovUnits.length_ref, # Maximum radius to solve up to in cm
    r_beg = 0.e5/tovUnits.length_ref,
);

# Function Setup for inverse EoS and TOV equation for the solver
EoS, EoS_inv = EoS_GCA2018();

# Setup initial condition; central pressure and enclosed mass
u0 = [EoS(Sim_Input.ρ0*tovUnits.rho_ref;)/tovUnits.pressure_ref; 0.0];

# Solve the TOV equation using `QSpin` for the stiff and soft EoSs.
@time TOV_sol = TOV_Solve_dimensionless(
    u0,
    Sim_Input.dr,
    Sim_Input.Dr,
    Sim_Input.r_beg,
    Sim_Input.r_end,
    EoS_inv;
    alg = OrdinaryDiffEqLowOrderRK.DP5(),
    reltol = 1e-8,
    abstol = [1.0, 1e15],
)
plot(TOV_sol.r*tovUnits.length_ref*1e-5, TOV_sol.ρr*tovUnits.rho_ref)
