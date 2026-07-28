
#=
This script demonstrates how to solve the Tolman-Oppenheimer-Volkoff (TOV) equation for neutron stars using the `QSpin` package.
See the documentation for the `QSpin.TOV` for more information.

The equation of state is from Negele & Vautherin (1973) and Graber, Cumming & Anderson (2018) in a two-component setup for the consideration of crust.
It is implemented in the `EoS_GCA2018` function in the CGS units only.
=#

using QSpin
using QSpin.Parameters: ParameterType
using QSpin.TOV: TOV_Solve
using QSpin.TOV.EquationOfState: EoS_GCA2018
using Plots, LaTeXStrings
using QSpin.PhysicalConstants: neutron_mass, mass_sun
Sim_Input = (
    ρ0 = 0.08*1e39*neutron_mass*1e3, # Initial central density in g/cm^3, above 2e16 seems to be unstable
    dr = 0.0001*1e5, # Radial step in cm
    Dr = 0.005*1e5, # Radial interval for recording values in cm
    r_beg = 10*1e5,
    M_core = 1.4 * mass_sun * 1e3,
    units = "CGS", # optional
);

# Calling EoS_GCA2018
EoS, EoS_inv = EoS_GCA2018();
u0 = [EoS(Sim_Input.ρ0); Sim_Input.M_core];

# Callback setup to terminate the integration when the pressure drops below zero
@time TOV_sol = TOV_Solve(
    EoS_inv,
    u0,
    Sim_Input.dr,
    Sim_Input.Dr,
    Sim_Input.r_beg;
    reltol = 1e-13,
    input_units = Sim_Input.units, # optional
    rho_ref = 2.8e14, # optional - nuclear saturation density in g/cm^3
)

plot(TOV_sol.r*1e-5, TOV_sol.ρr)
