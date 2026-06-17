using QSpin
using QSpin.Parameters: ParameterType
using QSpin.TOV: TOV_Solve
using QSpin.TOV.EquationOfState: EoS_two_component_polytrope

using QSpin.PhysicalConstants
using Plots, LaTeXStrings
import OrdinaryDiffEq as DE
using OrdinaryDiffEqLowOrderRK

file_path = "scripts/mutual_friction_input.json"
output = QSpin.MFriction.VNparaGraber2018(file_path)

ħ = hbar * 1e7;
mn = neutron_mass * 1e3;
Msun = mass_sun;

# Polytropic EoS parameters
EoS_Param_Stiff = (
    K_crust = (3*π^2)^(2/3) * ħ^2 / (5*mn^(8/3)), # Polytropic constant for the crust
    γ_crust = 5.0/3.0, # Polytropic index for the crust
    γ_core = 3.0, # Polytropic index for the core
    ρ_b = 3e14, # Transition density between crust and core in kg/m^3
)

# Input Parameters for the TOV solver
TOV_Input = (
    ρ0 = 1e15, # Initial central density in kg/m^3 (1e15 g/cm^3)
    dr = 0.005*1e5, # Radial step in meters
    Dr = 0.01*1e5, # Radial interval for recording values in meters
    r_end = 25e5, # Maximum radius to solve up to in meters
    r_beg = 0.0,
);

# Fucntion Setup for inverse EoS and TOV equation for the solver
EoS, EoS_inv = EoS_two_component_polytrope(EoS_Param_Stiff);
#EoS, EoS_inv = QSpin.TOV.EoS_GCA2018();
# Setting up initial condition accordingly to the EoS for a given central density ρ0.
u0 = [EoS_Stiff(TOV_Input.ρ0); 0.0]; # Initial conditions: central pressure and enclosed mass

# Sovling the TOV equation using the RK4 method with QSpin OdeSolve module for Stiff and Soft EoSs.

@time TOV_sol = TOV_Solve(
    u0,
    TOV_Input.dr,
    TOV_Input.Dr,
    TOV_Input.r_beg,
    TOV_Input.r_end,
    EoS_inv;
    alg = OrdinaryDiffEqLowOrderRK.DP5(),
    reltol = 1e-8,
    abstol = [1.0, 1e15],
)
Bs = QSpin.MFriction.MutualFrictionCoefficients(
    (ρs = TOV_sol.ρr, r = TOV_sol.r, Beb_core = 1e-4, Bj_core = 1e-4),
    output.Beb_itp,
    output.Bj_itp;
    Rcci = 1.2e4,
)
yBeb = Bs.Beb
yBj = Bs.Bj


Glitch_Raiser_Input = (
    B_core = 5e-4, # Mutual Friction Parameter
    B_sf = yBeb[1:(end-1)], # Mutual Friction Parameter
    Ω_crust = 70.34, # Initial angular velocity of the crust in rad/s
    Ω_sf = 70.34 + 6.3e-3, # Initial angular velocity of the superfluid in rad/s
    Ω_core = 70.34, # Initial angular velocity of the core in rad/s
    I_crust = 4.5e30, # Moment of inertia of the crust in kg m^2
    I_core = 0.8 * 4.5e30, # Moment of inertia of the core in kg m
    N_ext = 0.0,
    dt = 1e-6, # Time step for the ODE solver in seconds
    Dt = 0.1, # Time interval for recording values in the glitch model in seconds
    t_start = 0.0, # Start time for the glitch model simulation in seconds
    t_end = 120.0, # End time for the glitch model simulation in seconds
    ρr = TOV_sol.ρr[1:(end-1)],
    r = TOV_sol.r[1:(end-1)],
);

function eom!(dΩ::AbstractArray, Ω::AbstractArray, Param::ParameterType, time::Float64)
    # dΩ_sf/dt
    Ω_sf = Ω[3:end];
    #Bsf = Param.B_sf * Param.ρr ./ maximum(Param.ρr); # Scaling B_sf with the local density profile
    dΩ_sfdr = [diff(Ω_sf) ./ diff(Param.r); 0.0];
    dΩ[3:end] = Param.B_sf .* (2 * Ω_sf + Param.r .* dΩ_sfdr) .* (Ω[1] .- Ω_sf);
    dΩ_sf_net =
        4 .* π .* sum(
            ((Param.r .^ 2) .* Param.ρr .* dΩ[3:end]) .*
            [diff(Param.r); diff(Param.r)[end]],
        );
    # dΩ_core/dt
    dΩ[2] = 2 * Param.B_core * Ω[2] * (Ω[1] - Ω[2]);
    # dΩ_crust/dt
    dΩ[1] =
        -Param.N_ext/Param.I_crust - Param.I_core/Param.I_crust * dΩ[2] -
        dΩ_sf_net/Param.I_crust;
end

Ω0 = [
    Glitch_Raiser_Input.Ω_crust;
    Glitch_Raiser_Input.Ω_core;
    Glitch_Raiser_Input.Ω_sf*ones(length(TOV_sol.r)-1); # Adding small random perturbations to the superfluid angular velocity
];

@time sol = QSpin.OdeSolve.evolve(
    eom!,
    Ω0,
    0.0,
    Glitch_Raiser_Input.t_end,
    Glitch_Raiser_Input;
    alg = DE.Tsit5(), #KenCarp4() for stiff porblems, but may not be stable either.
    dt = Glitch_Raiser_Input.dt,
    saveat = Glitch_Raiser_Input.Dt,
    reltol = 1e-14,
)
Ωt = Array(sol);
t = sol.t;

plt1 = plot(
    t,
    [Ωt[1, :] Ωt[2, :] Ωt[length(TOV_sol.r)-75, :]],
    label = ["Crust" "Core" "Superfluid"],
    xlabel = "Time (s)",
    ylabel = L"\Omega\;(\mathrm{rad/s})",
    framestyle = :box,
    linewidth = 2,
)

plt2 = heatmap(
    t,
    TOV_sol.r[1:(end-1)]/1e5,
    Ωt[3:end, :],
    framestyle = :box,
    xlabel = "Time (s)",
    ylabel = "Radius (km)",
    ylims = (TOV_Input.r_beg/1e5, TOV_sol.R/1e5),
)

ρc_scan = exp10.(range(13, stop = log10(2e17), length = 150))*1e-3 # in g * cm^-3, which is equivalent to 1e-3 times the input in kg * fm^-3
yBewS = exp10.(output.Beb_itp(log10.(ρc_scan)))
yBjS = exp10.(output.Bj_itp(log10.(ρc_scan)))
plt3 = plot(
    ρc_scan,
    yBewS,
    ls = :dashdot,
    lc = :blue,
    linewidth = 2,
    xaxis = (:log10, [3.5e11, 2e14]),
    yaxis = (:log10, [1e-5, 1e-1]),
    label = L"\textrm{(B)} \; \mathcal{B}_\mathrm{EB} \; \textrm{with }\; E_{p} ",
)
plot!(plt3, scatter!(output.ρs, output.Beb, label = false, mc = :blue, lc = :blue))
plot!(
    plt3,
    plot!(
        ρc_scan,
        yBjS,
        ls = :dash,
        lc = RGB(0.94, 0.65, 0.25),
        linewidth = 2,
        xaxis = (:log10, [3.5e11, 2e14]),
        yaxis = (:log10, [1e-5, 1e-1]),
        label = L"\textrm{(C)} \; \mathcal{B}_\mathrm{J} \; \textrm{with }\; E_{p} ",
    ),
)
plot!(
    plt3,
    scatter!(
        output.ρs,
        output.Bj,
        label = false,
        mc = RGB(0.94, 0.65, 0.25),
        lc = RGB(0.94, 0.65, 0.25),
        framestyle = :box,
        legend = :bottomleft,
    ),
)
xlabel!(L"\rho_s \; (\textrm{g} \; \textrm{cm}^{-3})")
ylabel!(L"\mathcal{B}")

plt4 = plot(
    TOV_sol.r[TOV_sol.ρr .>= 4e11]*1e-5,
    yBeb[TOV_sol.ρr .>= 4e11],
    lc = :blue,
    label = L"\mathcal{B}_{EB}",
    xlabel = "Radius (km)",
    ylabel = L"\mathcal{B}",
    framestyle = :box,
    linewidth = 2,
    legend = :topleft,
    yaxis = (:log10),
)
plot!(
    plt4,
    plot!(
        TOV_sol.r[TOV_sol.ρr .>= 4e11]*1e-5,
        yBj[TOV_sol.ρr .>= 4e11],
        lc = RGB(0.94, 0.65, 0.25),
        label = L"\mathcal{B}_{J}",
        xlabel = "Radius (km)",
        ylabel = L"\mathcal{B}",
        framestyle = :box,
        linewidth = 2,
        legend = :topleft,
        yaxis = (:log10),
    ),
)

plot(plt1, plt2, plt3, plt4, layout = (2, 2))
