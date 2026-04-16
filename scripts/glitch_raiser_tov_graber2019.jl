using QSpin
using QSpin.Parameters: ParameterType
using QSpin.PhysicalConstants
using Plots, LaTeXStrings

ħ = hbar;
mn = neutron_mass;

# Polytropic EoS parameters
EoS_Param_Stiff = (
    γcore = 3.0, # Polytropic index for the core
    ρb = 3e14*1e3, # Transition density between crust and core in kg/m^3
)

# Simulation Input Parameters
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

function eom!(dψ::AbstractArray, ψ::AbstractArray, parameters::ParameterType, time::Float64)
    B = Glitch_Raiser_Input.B
    dΩ[1] = -B * (Ω[1] - Ω[2]) # dΩ_crust/dt
    dΩ[2] = B * (Ω[1] - Ω[2]) - B * (Ω[2] - Ω[3]) # dΩ_sf/dt
    dΩ[3] = B * (Ω[2] - Ω[3]) # dΩ_core/dt
end

Ωt = QSpin.OdeSolve.evolve(eom!, ψ0, t_start, t_end; alg = DE.Tsit5(), dt = dt, saveat = Dt)
