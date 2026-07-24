using QSpin
using QSpin.Parameters: ParameterType
using QSpin.TOV: TOV_Solve
using QSpin.TOV.EquationOfState: EoS_GCA2018
using QSpin.PhysicalConstants
using QSpin.GlitchModels: ThreeCompGCA2018!
using QSpin.PhysicalConstants: neutron_mass
using Plots, LaTeXStrings
import QSpin.OdeSolve: evolve
import CommonSolve as DE
using OrdinaryDiffEqLowOrderRK
import OrdinaryDiffEqTsit5: Tsit5

file_path = "scripts/mutual_friction_input.json"
output = QSpin.MFriction.VNparaGraber2018(file_path)

# Input Parameters for the TOV solver
Sim_Input = (
    # gltich model parameters
    B_core = 5e-4, # Mutual Friction Parameter
    Ω_crust = 70.34, # Initial angular velocity of the crust in rad/s
    Ω_sf = 70.34 + 6.3e-3, # Initial angular velocity of the superfluid in rad/s
    Ω_core = 70.34, # Initial angular velocity of the core in rad/s
    I_crust = 4.5e30, # Moment of inertia of the crust in kg m^2
    I_core = 0.8 * 4.5e30, # Moment of inertia of the core in kg m
    N_ext = 0.0,
    # Glitch model solver setup
    dt = 1e-6, # Time step for the ODE solver in seconds
    Dt = 0.1, # Time interval for recording values in the glitch model in seconds
    t_start = 0.0, # Start time for the glitch model simulation in seconds
    t_end = 120.0, # End time for the glitch model simulation in seconds
    # TOV solver parameters
    ρ0 = 0.08*1e39*neutron_mass*1e3, # Initial central density in g/cm^3 (GCS units)
    dr = 0.005*1e5, # Radial step in centimeters
    Dr = 0.01*1e5, # Radial interval for recording values in centimeters
    r_beg = 1e6, # Starting radius in centimeters
    M_core = 1.4e3 * mass_sun, # Mass of the neutron star core in grams
    tov_units = "CGS",
);

# Fucntion Setup for inverse EoS and TOV equation for the solver
EoS, EoS_inv = EoS_GCA2018();
# Setting up initial condition accordingly to the EoS for a given central density ρ0.
u0 = [EoS(Sim_Input.ρ0); 0.0]; # Initial conditions: central pressure and enclosed mass

# Sovling the TOV equation using the RK4 method with QSpin OdeSolve module for Stiff and Soft EoSs.
@time TOV_sol = TOV_Solve(
    EoS_inv,
    u0,
    Sim_Input.dr,
    Sim_Input.Dr,
    Sim_Input.r_beg;
    reltol = 1e-13,
    abstol = 1e-13,
    input_units = Sim_Input.tov_units, # optional
    # rho_ref = 2.8e14, # optional - nuclear saturation density in g/cm^3
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

Ω_ini = [Sim_Input.Ω_crust; Sim_Input.Ω_sf*ones(size(TOV_sol.r)); Sim_Input.Ω_core]

EoMSetup = (
    rho = TOV_sol.ρr,
    r = TOV_sol.r,
    M_NS = TOV_sol.M+1.4*mass_sun,
    R_NS = TOV_sol.R,
    R_cci = 1e6, # 10 km in cm
    B_core = 5e-4, # Mutual Friction Parameter
    B_sf = yBeb, # Mutual Friction Parameter
)

EoM! = ThreeCompGCA2018!(EoMSetup)
sol = evolve(
    EoM!,
    Ω_ini,
    0.0,
    Sim_Input.t_end,
    Sim_Input;
    alg = Tsit5(),
    saveat = Sim_Input.Dt,
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
    ylims = (Sim_Input.r_beg/1e5, TOV_sol.R/1e5),
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
