
#=
This script demonstrates how to solve the Tolman-Oppenheimer-Volkoff (TOV) equation for neutron stars using the `QSpin` package.
See the documentation for the `QSpin.TOV` for more information.

The equation of state is a two-component polytrope, with the parameters chosen from https://github.com/vanessagraber/teaching_materials/blob/master/summerschool_CRAQ_2019/mass_radius_relations.ipynb
=#

using QSpin
using QSpin.Parameters: ParameterType
using QSpin.TOV: TOV_Solve
using QSpin.TOV.EquationOfState: EoS_GCA2018, EoS_two_component_polytrope
using QSpin.PhysicalConstants:
    gravitational_constant, speed_of_light_vacuum, neutron_mass, hbar
using Plots, LaTeXStrings

ħ = hbar;
mn = neutron_mass;

CGS_Input = (
    ρ0 = 1e15, # Initial central density in g/cm^3, above 2e16 seems to be unstable
    dr = 0.0001*1e5, # Radial step in cm
    Dr = 0.005*1e5, # Radial interval for recording values in cm
    r_beg = 0.e5,
    units = "CGS", # CGS or CGS_dim
    K_crust = (3*π^2)^(2/3) * (ħ*1e7)^2 / (5*(neutron_mass*1e3)^(8/3)), # Polytropic constant for the crust
    γ_crust = 5.0/3.0, # Polytropic index for the crust
    γ_core = 3.0, # Polytropic index for the core
    ρ_b = 3e14, # Transition density between crust and core in kg/m^3
);

SI_Input = (
    ρ0 = 1e18, # Initial central density in kg/m^3, above 2e18 seems to be unstable
    dr = 0.0001*1e3, # Radial step in m
    Dr = 0.005*1e3, # Radial interval for recording values in m
    r_beg = 0.e3,
    units = "SI", # SI or SI_dim
    K_crust = (3*π^2)^(2/3) * ħ^2 / (5*neutron_mass^(8/3)), # Polytropic constant for the crust
    γ_crust = 5.0/3.0, # Polytropic index for the crust
    γ_core = 3.0, # Polytropic index for the core
    ρ_b = 3e14*1e3, # Transition density between crust and core in kg/m^3
);

Sim_Input = CGS_Input; # Choose between CGS_Input and SI_Input for sovling the TOV equation.

# Calling EoS_two_component_polytrope
EoS, EoS_inv = EoS_two_component_polytrope(Sim_Input);;
u0 = [EoS(Sim_Input.ρ0); 0.0];

# Callback setup to terminate the integration when the pressure drops below zero
@time TOV_sol = TOV_Solve(
    EoS_inv,
    u0,
    Sim_Input.dr,
    Sim_Input.Dr,
    Sim_Input.r_beg;
    reltol = 1e-13,
    abstol = 1e-13,
    input_units = Sim_Input.units,
    rho_ref = 2.8e14, # nuclear saturation density in g/cm^3
)

if Sim_Input.units == "CGS" || Sim_Input.units == "CGS_dim"
    plot!(TOV_sol.r*1e-5, TOV_sol.ρr)
else
    plot!(TOV_sol.r*1e-3, TOV_sol.ρr*1e-3)
end
