
using QSpin
using QSpin.Parameters: ParameterType
using QSpin.TOV: TOV_Solve
using QSpin.TOV.EquationOfState: EoS_Type
using QSpin.GlitchModels: ThreeCompGeneral!, gm_input
using Plots, LaTeXStrings
using QSpin.MFriction: mfrictionGraber2016
using QSpin.OdeSolve: evolve
using OrdinaryDiffEqTsit5: Tsit5
file_name = "src/TOV/EoS/WNewton/eos_SkXi450_28.0_40.00_-100.00_glitch.dat"

Sim_Input = (
    ρ0 = 1e15, # Initial central density in g/cm^3, above 2e16 seems to be unstable
    dr = 0.0001*1e6, # Radial step in cm
    Dr = 0.0005*1e6, # Radial interval for recording values in cm
    r_beg = 0.e5,
    units = "CGS", # optional
    B_sf = 1e-3, # Mutual Friction Parameter for the SF component
    N_ext = 0,
    Ω_crust = 2*π*11.195, # Initial angular velocity of the crust in rad/s
    Ω_sf = 2*π*11.195+6.3e-3, # Initial angular velocity of the superfluid in rad/s
    Ω_core = 2*π*11.195, # Initial angular velocity of the core in rad
    ρ_drip = 4e11, # cgs units
    t_end = 1e4, # Total simulation time in seconds
    Dt = 1e-1, # Time step for saving data in seconds
);

EoS, EoS_inv, EoS_nb, EoS_ρN, EoS_Yp, EoS_mp, EoS_kFe, EoS_kFn, ρ_cci =
    EoS_Type("LinterpSkyrme", FileInput = file_name);
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
    input_units = Sim_Input.units, # optional
    rho_ref = 2.8e14, # optional - nuclear saturation density in g/cm^3
);

mp = EoS_mp(TOV_sol.ρr)
mp[mp .< 0] .= 0.0;
kFn = EoS_kFn(TOV_sol.ρr)
kFn[kFn .< 0.0] .= 0.0;
kFe = EoS_kFe(TOV_sol.ρr)
kFe[kFe .< 0.0] .= 0.0;
Yp = EoS_Yp(TOV_sol.ρr)
Yp[Yp .< 0.0] .= 0.0;
nb = EoS_nb(TOV_sol.ρr)
nb[nb .< 0.0] .= 0.0;

Params =
    (kFn = kFn, kFe = kFe, nb = nb, Yp = Yp, ρ = TOV_sol.ρr, B_sf = 4e-4, a = 4.2, b = 2.7)

Bs = mfrictionGraber2016("p", Params)

#i_cci = findfirst(x->x<=ρ_cci, TOV_sol.ρr)
#i_drip = findfirst(x->x<=Sim_Input.ρ_drip, TOV_sol.ρr)
i_drip = argmin(abs.(TOV_sol.ρr .- Sim_Input.ρ_drip))
i_cci = argmin(abs.(TOV_sol.ρr .- ρ_cci))
# Setting initial conditions for the glitch model
Ω_ini = [
    Sim_Input.Ω_crust;
    Sim_Input.Ω_sf*ones(size(TOV_sol.r[(i_cci+1):i_drip]));
    Sim_Input.Ω_core*ones(size(TOV_sol.r[1:i_cci]))
]



# Wrap the parameters in the glith model
EoMSetup = (
    rho = TOV_sol.ρr,
    rhoN = EoS_ρN(TOV_sol.ρr),
    r = TOV_sol.r,
    M_NS = TOV_sol.M,
    R_NS = TOV_sol.R,
    ρ_cci = TOV_sol.ρr[i_cci], # 10 km in cm
    B_sf = Sim_Input.B_sf, # Mutual Friction Parameter
    B_core = Bs.B_core[1:i_cci], # Mutual Friction Parameter
    N_ext = Sim_Input.N_ext,
    ρ_drip = Sim_Input.ρ_drip, # Drip density in g/cm^3
)

EoM! = ThreeCompGeneral!(
    EoMSetup;
    ρ_drip = Sim_Input.ρ_drip,
    ρ_cci = EoMSetup.ρ_cci,
    value_check = true,
)
@time sol = evolve(
    EoM!,
    Ω_ini,
    0.0,
    Sim_Input.t_end,
    Sim_Input;
    alg = Tsit5(),
    saveat = Sim_Input.Dt,
    reltol = 1e-14,
    abstol = 1e-14,
)

# Data Reading and Plotting
Ωt = Array(sol);
t = sol.t;


plt1 = plot(
    t[2:end],
    (Ωt[1, 2:end] .- Ωt[1, 1])*1e6/2/π,
    xaxis = :log10,
    yaxis = :log10,
    xlabel = L"t\;(\textrm{s})",
    ylabel = L"\Delta\nu\;(\textrm{Hz})",
    xticks = [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4],
    framestyle = :box,
    lw = 2,
    ms = 2,
    color = :blue,
    label = "",
)
plt2 = heatmap(
    t[2:i_drip -i_cci+1],
    TOV_sol.r[i_cci:i_drip]*1e-5,
    Ωt[2:(i_drip-i_cci+1), 2:end],
    xaxis = :log10,
    yaxis = :log10,
    xlabel = L"t\;(\textrm{s})",
    ylabel = L"r\;(\textrm{mm})",
    xticks = [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4],
    framestyle = :box,
    color = :viridis,
    label = "",
)
