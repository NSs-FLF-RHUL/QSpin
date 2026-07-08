
#=
This script demonstrates how to solve the Tolman-Oppenheimer-Volkoff (TOV) equation for neutron stars using the `QSpin` package.
See the documentation for the `QSpin.TOV` for more information.

The parameters, and equation of state, are chosen from https://github.com/vanessagraber/teaching_materials/blob/master/summerschool_CRAQ_2019/mass_radius_relations.ipynb.
=#

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
    ρ0 = 6.5e14, # Initial central density in g/cm^3, above 1e15 seems to be unstable
    dr = 0.0001*1e5, # Radial step in cm
    Dr = 0.005*1e5, # Radial interval for recording values in cm
    r_beg = 0.e5/tovUnits.length_ref,
    units = "CGS",
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

#plot(TOV_sol.r*tovUnits.length_ref*1e-5, TOV_sol.Pr*tovUnits.pressure_ref)

# Solve the TOV equation using `QSpin` for the stiff and soft EoSs.
@time TOV_sol = TOV_Solve_dimensionless(
    u0,
    Sim_Input.dr,
    Sim_Input.Dr,
    Sim_Input.r_beg,
    Sim_Input.r_end,
    EoS_inv;
    alg = OrdinaryDiffEqLowOrderRK.DP5(),
    reltol = 1e-15,
    abstol = [1.0, 1e15],
)
plot(
    TOV_sol.r*tovUnits.length_ref*1e-5,
    TOV_sol.ρr*tovUnits.rho_ref,
    xlabel = "Radius (km)",
    ylabel = "Density (g/cm^3)",
    title = "TOV Solution for Neutron Star",
    legend = false,
)
