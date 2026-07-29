using QSpin
using QSpin.Parameters: ParameterType
using QSpin.TOV: TOV_Solve
using QSpin.TOV.EquationOfState: EoS_GCA2018
using QSpin.PhysicalConstants
using QSpin.GlitchModels: ThreeCompGCA2018!
using QSpin.MFriction: MutualFrictionCoefficients
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
    B_core = 1e-2, # Mutual Friction Parameter
    Ω_crust = 70.34, # Initial angular velocity of the crust in rad/s
    Ω_sf = 70.34 + 6.3e-3, # Initial angular velocity of the superfluid in rad/s
    Ω_core = 70.34, # Initial angular velocity of the core in rad/s
    I_crust = 4.5e30, # Moment of inertia of the crust in kg m^2
    I_core = 0.8 * 4.5e30, # Moment of inertia of the core in kg m
    N_ext = 0.0,
    # Glitch model solver setup
    dt = 1e-6, # Time step for the ODE solver in seconds
    Dt = 0.02, # Time interval for recording values in the glitch model in seconds
    t_start = 0.0, # Start time for the glitch model simulation in seconds
    t_end = 120.0, # End time for the glitch model simulation in seconds
    # TOV solver parameters
    ρ0 = 0.08*1e39*neutron_mass*1e3, # Initial central density in g/cm^3 (GCS units)
    dr = 0.0005*1e5, # Radial step in centimeters
    Dr = 0.001*1e5, # Radial interval for recording values in centimeters
    r_beg = 1e6, # Starting radius in centimeters
    M_core = 1.4e3 * mass_sun, # Mass of the neutron star core in grams
    tov_units = "CGS",
);

# Fucntion Setup for inverse EoS and TOV equation for the solver
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
        Beb_core = 1e-4,
        Bj_core = 1e-4,
        input_units = "CGS",
    ),
    output.Beb_itp,
    output.Bj_itp;
    Rcci = 1e4,
)
yBeb = Bs.Beb
yBj = Bs.Bj


Glitch_Raiser_Input = (
    B_core = Sim_Input.B_core, # Mutual Friction Parameter
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
    N_ext = 0.0,
)

EoM! = ThreeCompGCA2018!(EoMSetup; value_check = true)
sol = evolve(
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

Ωt = Array(sol);
t = sol.t;


Ω_sf = Ωt[2:length(TOV_sol.r), :]
t_plots = [0.12, 0.6, 3.0, 30]
idt = 1
plt1 = plot(
    TOV_sol.r[1:(end-1)]*1e-5,
    Ω_sf[:, idt],
    xflip = true,
    label = string(L"t=0"),
    linewidth = 2,
)
for pp = 1:length(t_plots)
    idt = Int64(t_plots[pp]/Sim_Input.Dt)+1
    plot!(
        plt1,
        TOV_sol.r[1:(end-1)]*1e-5,
        Ω_sf[:, idt],
        linewidth = 2,
        label = string(L"t=", t_plots[pp], L"\;\textrm{s}"),
        xflip = true,
        xlabel = L"r (\mathrm{km})",
        ylabel = L"Ω_\mathrm{sf}\;(\textrm{s}^{-1})",
    )
end
plot!(plt1, xflip = true, xlims = (10.0, 10.43))


plt2 = plot(t, Ωt[1, :], label = L"\Omega_\mathrm{crust}", linewidth = 2)
plot!(
    plt2,
    t,
    Ωt[end, :],
    linewidth = 2,
    label = L"\Omega_\mathrm{core}",
    xlabel = L"t\;(\mathrm{s})",
    ylabel = L"Ω\;(\textrm{s}^{-1})",
)

plt3 = heatmap(
    t,
    TOV_sol.r[1:(end-1)]*1e-5,
    Ω_sf,
    xlabel = L"t\quad(\textrm{s})",
    ylabel = L"r\quad(\textrm{km})",
)


l = @layout [a b; c]
plot(plt2, plt1, plt3, layout = l)
