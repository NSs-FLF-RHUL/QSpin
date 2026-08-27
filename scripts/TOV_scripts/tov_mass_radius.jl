#=
This example script provides the mass-radius relations for neutron starts by solving the TOV equation with two-component polytropic EoSs in Soft and Stiff format using QSpin package. The benchmark solutions can be referred from https://github.com/vanessagraber/teaching_materials/blob/master/summerschool_CRAQ_2019/mass_radius_relations.ipynb.

The TOV equation is solved in the dimensionless format with characteristic length scales set by the nuclear density, 2.8e14 g/cm^3, and the nuclear mass, 1.67e-24 g * cm^-3, which can be read from QSpin.TOV.TOV_ref_units. The units can be chosen arbitary by assignining different referecne density.

The equation of states are set in the CGS units.

=#

using QSpin
using QSpin.Parameters: ParameterType
using QSpin.TOV: TOV_Solve
using QSpin.TOV.EquationOfState: EoS_two_component_polytrope
using QSpin.PhysicalConstants: neutron_mass, mass_sun, hbar
using Plots, LaTeXStrings
import CommonSolve as DE
using OrdinaryDiffEqLowOrderRK

# Polytropic EoS parameters
EoS_Param_Stiff = (
    K_crust = (3*π^2)^(2/3) * (hbar*1e7)^2 / (5*(neutron_mass*1e3)^(8/3)),
    γ_crust = 5.0/3.0,
    γ_core = 3.0,
    ρ_b = 3e14,
)

EoS_Param_Soft = (
    K_crust = (3*π^2)^(2/3) * (hbar*1e7)^2 / (5*(neutron_mass*1e3)^(8/3)),
    γ_crust = 5.0/3.0,
    γ_core = 5.0/2.0,
    ρ_b = 3e14,
)

Sim_Input = (
    ρ0 = 1e15, # Initial central density in g/cm^3, above 1e15 seems to be unstable
    dr = 0.1, # Radial step in cm
    Dr = 1e2, # Radial interval for recording values in cm
    r_beg = 0.0,
);

# Function Setup for inverse EoS and TOV equation for the solver
EoS_Stiff, EoS_inv_Stiff = EoS_two_component_polytrope(EoS_Param_Stiff);
EoS_Soft, EoS_inv_Soft = EoS_two_component_polytrope(EoS_Param_Soft);

# Getting two primary examples for the same input parameters but different EoSs.
## Setup initial condition; central pressure and enclosed mass

u0_Stiff = [EoS_Stiff(Sim_Input.ρ0); 0.0];
u0_Soft = [EoS_Soft(Sim_Input.ρ0); 0.0];

# Solving the TOV equation for the stiff EoS with the same input parameters.
@time TOV_sol_Stiff = TOV_Solve(
    EoS_inv_Stiff,
    u0_Stiff,
    Sim_Input.dr,
    Sim_Input.Dr,
    Sim_Input.r_beg;
    alg = OrdinaryDiffEqLowOrderRK.DP5(),
    reltol = 1e-15,
)
# Solving the TOV equation for the soft EoS with the same input parameters.
@time TOV_sol_Soft = TOV_Solve(
    EoS_inv_Soft,
    u0_Soft,
    Sim_Input.dr,
    Sim_Input.Dr,
    Sim_Input.r_beg;
    alg = OrdinaryDiffEqLowOrderRK.DP5(),
    reltol = 1e-15,
)

# Compute the M-R relation by varying the central density, and solving the TOV equation for each case.
# The radius is defined by the first point that the pressure becomes negative.
ρc_scan = exp10.(range(14.5, stop = 17, length = 150))
M_StiffScan = zeros(length(ρc_scan));
R_StiffScan = zeros(length(ρc_scan));
M_SoftScan = zeros(length(ρc_scan));
R_SoftScan = zeros(length(ρc_scan));

for cc = 1:length(ρc_scan)
    ρc = ρc_scan[Int.(cc)]
    println(cc, ": Solving TOV eq. for central density ρc = ", ρc/1e14, "1e14 g/cm^3")
    TOV_sol = TOV_Solve(
        EoS_inv_Stiff,
        [EoS_Stiff(ρc); 0.0],
        Sim_Input.dr,
        Sim_Input.Dr,
        Sim_Input.r_beg;
        alg = OrdinaryDiffEqLowOrderRK.DP5(),
        reltol = 1e-15,
    )
    M_StiffScan[Int.(cc)] = TOV_sol.M
    R_StiffScan[Int.(cc)] = TOV_sol.R

    TOV_sol = TOV_Solve(
        EoS_inv_Soft,
        [EoS_Soft(ρc); 0.0],
        Sim_Input.dr,
        Sim_Input.Dr,
        Sim_Input.r_beg;
        alg = OrdinaryDiffEqLowOrderRK.DP5(),
        reltol = 1e-15,
    )
    M_SoftScan[Int.(cc)] = TOV_sol.M
    R_SoftScan[Int.(cc)] = TOV_sol.R
end

# Plotting two primary stiff and soft examples
Pc_Stiff = EoS_Stiff(Sim_Input.ρ0)
Pc_Soft = EoS_Soft(Sim_Input.ρ0)
plt1 = plot(
    TOV_sol_Stiff.r/1e5,
    TOV_sol_Stiff.Pr/Pc_Stiff,
    label = string(L"\gamma_\mathrm{core}=", EoS_Param_Stiff.γ_core),
    xlabel = "Radius (km)",
    ylabel = L"P/P_c",
    framestyle = :box,
    linewidth = 2,
)
plot!(
    plt1,
    TOV_sol_Soft.r/1e5,
    TOV_sol_Soft.Pr/Pc_Soft,
    label = string(L"\gamma_\mathrm{core}=", EoS_Param_Soft.γ_core),
    linewidth = 2,
)
plt2 = plot(
    TOV_sol_Stiff.r/1e5,
    TOV_sol_Stiff.ρr/Sim_Input.ρ0,
    label = string(L"\gamma_\mathrm{core}=", EoS_Param_Stiff.γ_core),
    xlabel = "Radius (km)",
    ylabel = L"\rho/\rho_c",
    framestyle = :box,
    linewidth = 2,
)
plot!(
    plt2,
    TOV_sol_Soft.r/1e5,
    TOV_sol_Soft.ρr/Sim_Input.ρ0,
    label = string(L"\gamma_\mathrm{core}=", EoS_Param_Soft.γ_core),
    linewidth = 2,
)
plt3 = plot(
    TOV_sol_Stiff.r/1e5,
    TOV_sol_Stiff.mr*1e-3/mass_sun,
    label = string(L"\gamma_\mathrm{core}=", EoS_Param_Stiff.γ_core),
    xlabel = "Radius (km)",
    ylabel = L"M/M_\odot",
    framestyle = :box,
    linewidth = 2,
)
plot!(
    plt3,
    TOV_sol_Soft.r/1e5,
    TOV_sol_Soft.mr*1e-3/mass_sun,
    label = string(L"\gamma_\mathrm{core}=", EoS_Param_Soft.γ_core),
    linewidth = 2,
)
plt4 = scatter(
    [R_StiffScan R_SoftScan]/1e5,
    [M_StiffScan M_SoftScan]*1e-3/mass_sun,
    ;
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
