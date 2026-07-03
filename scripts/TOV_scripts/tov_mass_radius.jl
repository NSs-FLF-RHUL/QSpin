#=
This example script provides the mass-radius relations for neutron starts by solving the TOV equation with two-component polytropic EoSs in Soft and Stiff format using QSpin package. The benchmark solutions can be referred from https://github.com/vanessagraber/teaching_materials/blob/master/summerschool_CRAQ_2019/mass_radius_relations.ipynb.

The TOV equation is solved in the dimensionless format with characteristic length scales set by the nuclear density, 2.8e14 g/cm^3, and the nuclear mass, 1.67e-24 g * cm^-3, which can be read from QSpin.TOV.TOV_ref_units. The units can be chosen arbitary by assignining different referecne density.

The equation of states are set in the CGS units.

=#

using QSpin
using QSpin.Parameters: ParameterType
using QSpin.TOV: TOV_Solve_dimensionless, TOV_ref_units, tov_eq_dimless!
using QSpin.TOV.EquationOfState: EoS_two_component_polytrope
using QSpin.PhysicalConstants: neutron_mass, mass_sun, hbar
using Plots, LaTeXStrings
import OrdinaryDiffEq as DE
using OrdinaryDiffEqLowOrderRK

tovUnits = TOV_ref_units()

# Polytropic EoS parameters
EoS_Param_Stiff = (
    K_crust = (3*π^2)^(2/3) * (hbar*1e7)^2 / (5*(neutron_mass*1e3)^(8/3)),
    γ_crust = 5.0/3.0,
    γ_core = 3.0,
    ρ_b = 3e14,
)
EoS_Param_Soft = (
    K_crust = (3*π^2)^(2/3) * (hbar*1e7)^2 / (5*(neutron_mass*1e3)^(8/3)),
    γ_crust = 5.0/3.0,
    γ_core = 5.0/2.0,
    ρ_b = 3e14,
)

Sim_Input = (
    ρ0 = 1e15 / tovUnits.rho_ref, # Initial central density in g/cm^3, above 1e15 seems to be unstable
    dr = 10/tovUnits.length_ref, # Radial step in cm
    Dr = 1e2/tovUnits.length_ref, # Radial interval for recording values in cm
    r_beg = 0.e5/tovUnits.length_ref,
);

# Function Setup for inverse EoS and TOV equation for the solver
EoS_Stiff, EoS_inv_Stiff = EoS_two_component_polytrope(EoS_Param_Stiff);
EoS_Soft, EoS_inv_Soft = EoS_two_component_polytrope(EoS_Param_Soft);

# Getting two primary examples for the same input parameters but different EoSs.
## Setup initial condition; central pressure and enclosed mass
u0_Stiff = [EoS_Stiff(Sim_Input.ρ0*tovUnits.rho_ref)/tovUnits.pressure_ref; 0.0];
u0_Soft = [EoS_Soft(Sim_Input.ρ0*tovUnits.rho_ref)/tovUnits.pressure_ref; 0.0];

# Solving the TOV equation for the stiff EoS with the same input parameters.
@time TOV_sol_Stiff = TOV_Solve_dimensionless(
    u0_Stiff,
    Sim_Input.dr,
    Sim_Input.Dr,
    Sim_Input.r_beg,
    EoS_inv_Stiff;
    alg = OrdinaryDiffEqLowOrderRK.DP5(),
    reltol = 1e-15,
)
# Solving the TOV equation for the soft EoS with the same input parameters.
@time TOV_sol_Soft = TOV_Solve_dimensionless(
    u0_Soft,
    Sim_Input.dr,
    Sim_Input.Dr,
    Sim_Input.r_beg,
    EoS_inv_Soft;
    alg = OrdinaryDiffEqLowOrderRK.DP5(),
    reltol = 1e-15,
)

Pc_Stiff = EoS_Stiff(Sim_Input.ρ0*tovUnits.rho_ref)
Pc_Soft = EoS_Soft(Sim_Input.ρ0*tovUnits.rho_ref)
plt1 = plot(
    TOV_sol_Stiff.r*tovUnits.length_ref/1e5,
    TOV_sol_Stiff.Pr*tovUnits.pressure_ref/Pc_Stiff,
    label = string(L"\gamma_\mathrm{core}=", EoS_Param_Stiff.γ_core),
    xlabel = "Radius (km)",
    ylabel = L"P/P_c",
    framestyle = :box,
    linewidth = 2,
)
plot!(
    plt1,
    TOV_sol_Soft.r*tovUnits.length_ref/1e5,
    TOV_sol_Soft.Pr*tovUnits.pressure_ref/Pc_Soft,
    label = string(L"\gamma_\mathrm{core}=", EoS_Param_Soft.γ_core),
    linewidth = 2,
)
