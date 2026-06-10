#=
This script demonstrates how to solve the Tolman-Oppenheimer-Volkoff (TOV) equation for neutron stars using the `QSpin` package.
See the documentation for the `QSpin.TOV` for more information.

The parameters, and equation of state, are chosen from https://github.com/vanessagraber/teaching_materials/blob/master/summerschool_CRAQ_2019/mass_radius_relations.ipynb.
=#

using QSpin
using QSpin.Parameters: ParameterType
using QSpin.PhysicalConstants: hbar, neutron_mass
using QSpin.TOV: TOV_Solve
using QSpin.TOV.EquationOfState: EoS_two_component_polytrope
using Plots, LaTeXStrings
import OrdinaryDiffEq as DE

# Polytropic EoS parameters
EoS_Param_Stiff = (
    K_crust = (3*π^2)^(2/3) * hbar^2 / (5*neutron_mass^(8/3)),
    γ_crust = 5.0/3.0,
    γ_core = 3.0,
    ρ_b = 3e14*1e3,
)
EoS_Param_Soft = (
    K_crust = (3*π^2)^(2/3) * hbar^2 / (5*neutron_mass^(8/3)),
    γ_crust = 5.0/3.0,
    γ_core = 5.0/2.0,
    ρ_b = 3e14*1e3,
)

# Simulation Input Parameters
Sim_Input = (
    ρ0 = 1e15*1e3, # Initial central density in kg/m^3 (1e15 g/cm^3)
    dr = 0.005*1e3, # Radial step in meters
    Dr = 0.01*1e3, # Radial interval for recording values in meters
    r_end = 20e3, # Maximum radius to solve up to in meters
);

# Function Setup for inverse EoS and TOV equation for the solver
EoS_Stiff, EoS_inv_Stiff = EoS_two_component_polytrope(EoS_Param_Stiff);
EoS_Soft, EoS_inv_Soft = EoS_two_component_polytrope(EoS_Param_Soft);

# Setup initial condition; central pressure and enclosed mass
u0_Stiff = [EoS_Stiff(Sim_Input.ρ0); 0.0];
u0_Soft = [EoS_Soft(Sim_Input.ρ0); 0.0];

# Solve the TOV equation using `QSpin` for the stiff and soft EoSs.
@time TOV_sol_Stiff = TOV_Solve(
    u0_Stiff,
    Sim_Input.dr,
    Sim_Input.Dr,
    15e3,
    EoS_inv_Stiff;
    alg = DE.Tsit5(),
    reltol = 1e-8,
    abstol = [1.0, 1e15],
)
@time TOV_sol_Soft = TOV_Solve(
    u0_Soft,
    Sim_Input.dr,
    Sim_Input.Dr,
    15e3,
    EoS_inv_Soft;
    alg = DE.Tsit5(),
    reltol = 1e-8,
    abstol = [1.0, 1e15],
)

# Compute the M-R relation by varying the central density, and solving the TOV equation for each case.
# The radius is defined by the first point that the pressure becomes negative.
ρc_scan = exp10.(range(17.5, stop = 20, length = 150))
M_StiffScan = zeros(length(ρc_scan));
R_StiffScan = zeros(length(ρc_scan));
M_SoftScan = zeros(length(ρc_scan));
R_SoftScan = zeros(length(ρc_scan));

for cc = 1:length(ρc_scan)
    ρc = ρc_scan[Int.(cc)];
    println(cc, ": Solving TOV for central density ρc = ", ρc/1e18, "1e18 kg/m^3")

    TOV_sol = TOV_Solve(
        [EoS_Stiff(ρc); 0.0],
        Sim_Input.dr,
        Sim_Input.Dr,
        Sim_Input.r_end,
        EoS_inv_Stiff;
        alg = DE.Tsit5(),
        reltol = 1e-8,
        abstol = [1.0, 1e15],
    );
    M_StiffScan[Int.(cc)] = TOV_sol.M
    R_StiffScan[Int.(cc)] = TOV_sol.R

    TOV_sol = TOV_Solve(
        [EoS_Soft(ρc); 0.0],
        Sim_Input.dr,
        Sim_Input.Dr,
        Sim_Input.r_end,
        EoS_inv_Soft;
        alg = DE.Tsit5(),
        reltol = 1e-8,
        abstol = [1.0, 1e15],
    );
    M_SoftScan[Int.(cc)] = TOV_sol.M
    R_SoftScan[Int.(cc)] = TOV_sol.R
end

# Create plots to display the results.

# Fetch reference core pressure Pc
Pc_Stiff = EoS_Stiff(Sim_Input.ρ0)
Pc_Soft = EoS_Soft(Sim_Input.ρ0)

plt1 = plot(
    TOV_sol_Stiff.r/1e3,
    TOV_sol_Stiff.Pr,
    label = string(L"\gamma_\mathrm{core}=", EoS_Param_Stiff.γ_core),
    xlabel = "Radius (km)",
    ylabel = L"P/P_c",
    framestyle = :box,
    linewidth = 2,
)
plot!(
    plt1,
    TOV_sol_Soft.r/1e3,
    TOV_sol_Soft.Pr,
    label = string(L"\gamma_\mathrm{core}=", EoS_Param_Soft.γ_core),
    linewidth = 2,
)
plt2 = plot(
    TOV_sol_Stiff.r/1e3,
    TOV_sol_Stiff.ρr,
    label = string(L"\gamma_\mathrm{core}=", EoS_Param_Stiff.γ_core),
    xlabel = "Radius (km)",
    ylabel = L"\rho/\rho_c",
    framestyle = :box,
    linewidth = 2,
)
plot!(
    plt2,
    TOV_sol_Soft.r/1e3,
    TOV_sol_Soft.ρr,
    label = string(L"\gamma_\mathrm{core}=", EoS_Param_Soft.γ_core),
    linewidth = 2,
)
plt3 = plot(
    TOV_sol_Stiff.r/1e3,
    TOV_sol_Stiff.mr,
    label = string(L"\gamma_\mathrm{core}=", EoS_Param_Stiff.γ_core),
    xlabel = "Radius (km)",
    ylabel = L"\rho/\rho_c",
    framestyle = :box,
    linewidth = 2,
)
plot!(
    plt3,
    TOV_sol_Soft.r/1e3,
    TOV_sol_Soft.mr,
    label = string(L"\gamma_\mathrm{core}=", EoS_Param_Soft.γ_core),
    linewidth = 2,
)
plt4 = scatter(
    [R_StiffScan R_SoftScan]/1e3,
    [M_StiffScan M_SoftScan]/mass_sun;
    zcolor = log.([ρc_scan ρc_scan])/log(10), # color the data by the core density
    markershape = [:circle :diamond],
    bg = :linen,
    markersize = [2 2],
    label = ["Stiff" "Soft"],
    xlabel = "Radius (km)",
    ylabel = L"\textrm{Mass}\;(M_\odot)",
    framestyle = :box,
    linewidth = 2,
)
plot(plt1, plt2, plt3, plt4, layout = (2, 2))
