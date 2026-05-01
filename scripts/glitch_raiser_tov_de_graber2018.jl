using QSpin
using QSpin.Parameters: ParameterType
using QSpin.PhysicalConstants
using Plots, LaTeXStrings
import OrdinaryDiffEq as DE
ħ = hbar;
mn = neutron_mass;
Msun = mass_sun;

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
    r_end = 15e3, # Maximum radius to solve up to in meters
);



# Fucntion Setup for inverse EoS and TOV equation for the solver
EoS_Stiff, EoS_inv_Stiff = QSpin.TOV.EoS_two_component_polytrope(EoS_Param_Stiff);
# Setting up initial condition accordingly to the EoS for a given central density ρ0.
u0_Stiff = [EoS_Stiff(TOV_Input.ρ0); 0.0]; # Initial conditions: central pressure and enclosed mass
tov! = QSpin.TOV.tov_eq!(EoS_inv_Stiff);
# Sovling the TOV equation using the RK4 method with QSpin OdeSolve module for Stiff and Soft EoSs.

@time sol_tov = QSpin.OdeSolve.evolve(
    tov!,
    u0_Stiff,
    0.0,
    TOV_Input.r_end,
    ;
    dt = TOV_Input.dr,
    saveat = TOV_Input.Dr,
);
ur = Array(sol_tov);
r = sol_tov.t;
Pr = ur[1, :];
mr = ur[2, :];
ρr = EoS_inv_Stiff.(Pr);


Glitch_Raiser_Input = (
    B_core = 5e-4, # Mutual Friction Parameter
    B_sf = 1e-3, # Mutual Friction Parameter
    Ω_crust = 70.34, # Initial angular velocity of the crust in rad/s
    Ω_sf = 70.34 + 6.3e-3, # Initial angular velocity of the superfluid in rad/s
    Ω_core = 70.34, # Initial angular velocity of the core in rad/s
    I_crust = 4.5e30, # Moment of inertia of the crust in kg m^2
    I_core = 0.8 * 4.5e30, # Moment of inertia of the core in kg m
    N_ext = 0.0,
    dt = 1e-4, # Time step for the ODE solver in seconds
    Dt = 1.0, # Time interval for recording values in the glitch model in seconds
    t_start = 0.0, # Start time for the glitch model simulation in seconds
    t_end = 120.0, # End time for the glitch model simulation in seconds
    ρr = ρr,
    r = r,
);

function eom!(δΩ::AbstractArray, Ω::AbstractArray, Param::ParameterType, time::Float64)
    dΩ = zeros(size(Ω));
    # dΩ_sf/dt
    Ω_sf = Ω[3:end];

    dΩ_sfdr = [diff(Ω_sf) ./ diff(r); 0.0];
    dΩ[3:end] = Param.B_sf .* (2 * Ω_sf + r .* dΩ_sfdr) .* (Ω[1] .- Ω_sf);
    dΩ_sf_net = 4 .* π .* sum(((r .^ 2) .* Param.ρr .* dΩ[3:end]) .* [diff(r); diff(r)[1]]);
    # dΩ_core/dt
    dΩ[2] = 2 * Param.B_core * Ω[2] * (Ω[1] - Ω[2]);
    # dΩ_crust/dt
    dΩ[1] =
        -Param.N_ext/Param.I_crust - Param.I_core/Param.I_crust * dΩ[2] -
        dΩ_sf_net/Param.I_crust;
    return dΩ
end

Ω0 = [
    Glitch_Raiser_Input.Ω_crust;
    Glitch_Raiser_Input.Ω_core;
    Glitch_Raiser_Input.Ω_sf*ones(length(r)) + 5 * (rand(length(r)) .- 0.5) * 1e-6; # Adding small random perturbations to the superfluid angular velocity
];

@time sol = QSpin.OdeSolve.evolve(
    eom!,
    Ω0,
    0.0,
    Glitch_Raiser_Input.t_end,
    Glitch_Raiser_Input;
    dt = Glitch_Raiser_Input.dt,
    saveat = Glitch_Raiser_Input.Dt,
)
Ωt = Array(sol);
t = sol.t;


plot(
    plot(
        t,
        [Ωt[1, :] Ωt[2, :] Ωt[end, :]],
        label = ["Crust" "Core" "Superfluid"],
        xlabel = "Time (s)",
        ylabel = L"\Omega\;(\mathrm{rad/s})",
        framestyle = :box,
        linewidth = 2,
    ),
    heatmap(
        t,
        r/1e3,
        Ωt[3:end, :],
        framestyle = :box,
        xlabel = "Time (s)",
        ylabel = "Radius (km)",
    ),
    layout = (1, 2),
)
