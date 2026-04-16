using QSpin
using QSpin.Parameters: ParameterType
using QSpin.PhysicalConstants
using Plots, LaTeXStrings

ħ = hbar;
mn = neutron_mass;

# Polytropic EoS parameters
EoS_Param_Stiff = (
    K_crust = (3*π^2)^(2/3) * ħ^2 / (5*mn^(8/3)), # Polytropic constant for the crust
    γ_crust = 5.0/3.0, # Polytropic index for the crust
    γ_core = 3.0, # Polytropic index for the core
    ρ_b = 3e14*1e3, # Transition density between crust and core in kg/m^3
)

# Input Parameters for the TOV solver
TOV_Input = (
    ρ0 = 1e15*1e3, # Initial central density in kg/m^3 (1e15 g/cm^3)
    dr = 0.01*1e3, # Radial step in meters
    Dr = 0.01*1e3, # Radial interval for recording values in meters
    r_end = 20e3, # Maximum radius to solve up to in meters
);

Glitch_Raiser_Input = (
    B = 5e-4, # Mutual Friction Parameter
    Ω_crust = 70.34, # Initial angular velocity of the crust in rad/s
    Ω_sf = 70.34 - 6.3e-3, # Initial angular velocity of the superfluid in rad/s
    Ω_core = 70.34, # Initial angular velocity of the core in rad/s
    dt = 1e-3, # Time step for the ODE solver in seconds
    Dt = 1e-1, # Time interval for recording values in the glitch model in seconds
    t_start = 0.0, # Start time for the glitch model simulation in seconds
    t_end = 120.0, # End time for the glitch model simulation in seconds
);

# Fucntion Setup for inverse EoS and TOV equation for the solver
EoS_Stiff, EoS_inv_Stiff = QSpin.TOV.EoS_two_component_polytrope(EoS_Param_Stiff);
# Setting up initial condition accordingly to the EoS for a given central density ρ0.
u0_Stiff = [EoS_Stiff(TOV_Input.ρ0); 0.0]; # Initial conditions: central pressure and enclosed mass
# Sovling the TOV equation using the RK4 method with QSpin OdeSolve module for Stiff and Soft EoSs.
@time Pr, mr, ρr, M, R, r = QSpin.OdeSolve.TOV_Solve_rk4(
    u0_Stiff,
    TOV_Input.dr,
    TOV_Input.Dr,
    TOV_Input.r_end,
    EoS_inv_Stiff,
);
