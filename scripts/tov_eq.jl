using QSpin
using QSpin.Parameters: ParameterType
using QSpin.PhysicalConstants
using Plots, LaTeXStrings
import OrdinaryDiffEq as DE
# This script is demonstrating how to solve the Tolman–Oppenheimer–Volkoff (TOV) equation for neutron stars using the QSpin package with the built-in TOV solver using Runge-Kutta 4th order method.
# The parameters are chosen from https://github.com/vanessagraber/teaching_materials/blob/master/summerschool_CRAQ_2019/mass_radius_relations.ipynb.
ħ = hbar;
mn = neutron_mass;
# Polytropic EoS parameters
EoS_Param_Stiff = (
    K_crust = (3*π^2)^(2/3) * ħ^2 / (5*mn^(8/3)), # Polytropic constant for the crust
    γ_crust = 5.0/3.0, # Polytropic index for the crust
    γ_core = 3.0, # Polytropic index for the core
    ρ_b = 3e14*1e3, # Transition density between crust and core in kg/m^3
)
EoS_Param_Soft = (
    K_crust = (3*π^2)^(2/3) * ħ^2 / (5*mn^(8/3)), # Polytropic constant for the crust
    γ_crust = 5.0/3.0, # Polytropic index for the crust
    γ_core = 5.0/2.0, # Polytropic index for the core
    ρ_b = 3e14*1e3, # Transition density between crust and core in kg/m^3
)
# Simulation Input Parameters
Sim_Input = (
    ρ0 = 1e15*1e3, # Initial central density in kg/m^3 (1e15 g/cm^3)
    dr = 0.005*1e3, # Radial step in meters
    Dr = 0.01*1e3, # Radial interval for recording values in meters
    r_end = 20e3, # Maximum radius to solve up to in meters
);

# Fucntion Setup for inverse EoS and TOV equation for the solver
EoS_Stiff, EoS_inv_Stiff = QSpin.TOV.EoS_two_component_polytrope(EoS_Param_Stiff);
EoS_Soft, EoS_inv_Soft = QSpin.TOV.EoS_two_component_polytrope(EoS_Param_Soft);

# Setting up initial condition accordingly to the EoS for a given central density ρ0.
u0_Stiff = [EoS_Stiff(Sim_Input.ρ0); 0.0]; # Initial conditions: central pressure and enclosed mass
u0_Soft = [EoS_Soft(Sim_Input.ρ0); 0.0]; # Initial conditions: central pressure and enclosed mass

# Sovling the TOV equation using the RK4 method with QSpin OdeSolve module for Stiff and Soft EoSs.


@time TOV_sol_Stiff = QSpin.TOV.TOV_Solve(
    u0_Stiff,
    Sim_Input.dr,
    Sim_Input.Dr,
    15e3,
    EoS_inv_Stiff;
    alg = DE.Tsit5(),
    reltol = 1e-8,
    abstol = [1.0, 1e15],
)

@time TOV_sol_Soft = QSpin.TOV.TOV_Solve(
    u0_Soft,
    Sim_Input.dr,
    Sim_Input.Dr,
    15e3,
    EoS_inv_Soft;
    alg = DE.Tsit5(),
    reltol = 1e-8,
    abstol = [1.0, 1e15],
)


# Computing the M-R relation by varying the central density and solving the TOV equation for each case. The radius is defined by the first point that the pressure becomes negative.
ρc_scan = exp10.(range(17.5, stop = 20, length = 150))
M_StiffScan = zeros(length(ρc_scan));
R_StiffScan = zeros(length(ρc_scan));
M_SoftScan = zeros(length(ρc_scan));
R_SoftScan = zeros(length(ρc_scan));

for cc = 1:length(ρc_scan)

    ρc = ρc_scan[Int.(cc)];
    println(cc, ": Solving TOV for central density ρc = ", ρc/1e18, "1e18 kg/m^3")
    TOV_sol = QSpin.TOV.TOV_Solve(
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
    # println("Mass: ", M/PhysConst.Msun, " M_sun, Radius: ", R/1e3, " km")
    TOV_sol = QSpin.TOV.TOV_Solve(
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

# Plotting
Pc_Stiff = EoS_Stiff(Sim_Input.ρ0) # Getting the reference core pressure for the stiff case.
Pc_Soft = EoS_Soft(Sim_Input.ρ0) # Getting the reference core pressure for the soft case.


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
