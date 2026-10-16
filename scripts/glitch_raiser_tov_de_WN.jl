
using QSpin
using QSpin.Parameters: ParameterType
using QSpin.TOV: TOV_Solve
using QSpin.TOV.EquationOfState: EoS_Type
using QSpin.GlitchModels: ThreeCompGeneral!, gm_input
using Plots, LaTeXStrings
using QSpin.MFriction: mfrictionGraber2016
file_name = "src/TOV/EoS/WNewton/eos_SkXi450_28.0_40.00_-100.00_glitch.dat"

Sim_Input = (
    ρ0 = 1e15, # Initial central density in g/cm^3, above 2e16 seems to be unstable
    dr = 0.0001*1e5, # Radial step in cm
    Dr = 0.005*1e5, # Radial interval for recording values in cm
    r_beg = 0.e5,
    units = "CGS", # optional
    B_sf = 1e-3, # Mutual Friction Parameter for the SF component
    N_ext = 0,
    Ω_crust = 70.0, # Initial angular velocity of the crust in rad/s
    Ω_sf = 70.0, # Initial angular velocity of the superfluid in rad/s
    Ω_core = 70.0, # Initial angular velocity of the core in rad
    ρ_drip = 4e11, # cgs units
    t_end = 1e6, # Total simulation time in seconds
    Dt = 1e3, # Time step for saving data in seconds
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

icci = findfirst(x->x<=ρ_cci, TOV_sol.ρr)
idrip = findfirst(x->x<=Sim_Input.ρ_drip, TOV_sol.ρr)

# Setting initial conditions for the glitch model
Ω_ini = [
    Sim_Input.Ω_crust;
    Sim_Input.Ω_sf*ones(size(TOV_sol.r[(icci+1):idrip]));
    Sim_Input.Ω_core*ones(size(TOV_sol.r[1:icci]))
]



# Wrap the parameters in the glith model
EoMSetup = (
    rho = TOV_sol.ρr,
    rhoN = TOV_sol.ρr,
    r = TOV_sol.r,
    M_NS = TOV_sol.M,
    R_NS = TOV_sol.R,
    R_cci = TOV_sol.r[icci], # 10 km in cm
    B_sf = Sim_Input.B_sf, # Mutual Friction Parameter
    B_core = Bs.B_core, # Mutual Friction Parameter
    N_ext = Sim_Input.N_ext,
    ρ_drip = Sim_Input.ρ_drip, # Drip density in g/cm^3
)

EoM! = ThreeCompGeneral!(EoMSetup; ρ_drip = Sim_Input.ρ_drip, value_check = true)
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
