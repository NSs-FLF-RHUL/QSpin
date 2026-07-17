
using QSpin
using QSpin.Parameters: ParameterType
using QSpin.TOV: TOV_Solve
using QSpin.TOV.EquationOfState: EoS_LInterp
using Plots, LaTeXStrings

file_name = "src/TOV/EoS/WNewton/eos_SkXi450_28.0_40.00_-100.00_glitch.dat"

Sim_Input = (
    ρ0 = 5e15, # Initial central density in g/cm^3, above 2e16 seems to be unstable
    dr = 0.0001*1e5, # Radial step in cm
    Dr = 0.005*1e5, # Radial interval for recording values in cm
    r_beg = 0.e5,
    units = "CGS", # optional
);

EoS, EoS_inv = EoS_LInterp(file_name, (1, 2));
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
)

plot!(TOV_sol.r*1e-5, TOV_sol.ρr)
