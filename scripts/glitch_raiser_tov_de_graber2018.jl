using QSpin
using QSpin.Parameters: ParameterType
using QSpin.TOV: TOV_Solve
using QSpin.TOV.EquationOfState: EoS_GCA2018
using QSpin.PhysicalConstants
using QSpin.GlitchModels: ThreeCompGCA2018!, gm_input
using QSpin.MFriction: MutualFrictionCoefficients
using QSpin.PhysicalConstants: neutron_mass
using Plots, LaTeXStrings
import QSpin.OdeSolve: evolve
import CommonSolve as DE
using JSON
import OrdinaryDiffEqTsit5: Tsit5

save_path = "local_tests"
gm_input_path = "scripts/glitch_riser_input.json"
mf_input_path = "scripts/mutual_friction_input.json"

# Input Parameters for the TOV solver
Sim_Input = gm_input(gm_input_path)
MF_Output = QSpin.MFriction.VNparaGraber2018(mf_input_path)

# Setup for inverse EoS and TOV equation for the solver
EoS, EoS_inv = EoS_GCA2018();
# Setting up initial condition accordingly to the EoS for a given central density ρ0.
u0 = [EoS(Sim_Input.ρ0); Sim_Input.M_core]; # Initial conditions: central pressure and enclosed mass

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
    (
        ρs = TOV_sol.ρr*1e3,
        r = TOV_sol.r*1e-2,
        Beb_core = Sim_Input.B_core,
        Bj_core = Sim_Input.B_core,
        input_units = Sim_Input.tov_units,
    ),
    MF_Output.B_itp;
    R_cci = Sim_Input.R_cci,
    ρ_b = MF_Output.ρ[1],
)

# Setting initial conditions for the glitch model
Ω_ini = [Sim_Input.Ω_crust; Sim_Input.Ω_sf*ones(size(TOV_sol.r)); Sim_Input.Ω_core]

# Wrap the parameters in the glith model
EoMSetup = (
    rho = TOV_sol.ρr,
    r = TOV_sol.r,
    M_NS = TOV_sol.M,
    R_NS = TOV_sol.R,
    R_cci = Sim_Input.R_cci, # 10 km in cm
    B_core = Sim_Input.B_core, # Mutual Friction Parameter
    B_sf = Bs[Sim_Input.B_sf_type], # Mutual Friction Parameter
    N_ext = Sim_Input.N_ext,
)

EoM! = ThreeCompGCA2018!(EoMSetup; value_check = true)
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

# Data JSON output
jwrap_out = Dict(
    "t" => t,
    "r" => TOV_sol.r,
    "ρ" => TOV_sol.ρr,
    "Ω_crust" => Ωt[1, :],
    "Ω_sf" => Ωt[2:(length(TOV_sol.r)+1), :],
    "Ω_core" => Ωt[end, :],
    "B_core" => Sim_Input.B_core,
)
jsave_name = string(save_path, "/", "glitch_riser_output.json")
open(jsave_name, "w") do io
    JSON.print(io, jwrap_out, 4)
end

# Data Plotting
include("plotting/glitch_riser_summary_plot.jl")
