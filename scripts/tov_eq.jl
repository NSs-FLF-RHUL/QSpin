import OrdinaryDiffEq as DE

using QSpin
using QSpin.Parameters: ParameterType
using QSpin.PhysicalConstants: hbar, sun_mass, neutron_mass
using Plots, LaTeXStrings
# This script is demonstrating how to solve the Tolman–Oppenheimer–Volkoff (TOV) equation for neutron stars using the QSpin package with the built-in TOV solver using Runge-Kutta 4th order method.
# The parameters are chosen from https://github.com/vanessagraber/teaching_materials/blob/master/summerschool_CRAQ_2019/mass_radius_relations.ipynb.

"""
Helper function to return TOV-specific fields from `QSpin.OdeSolve.evolve` solution object.
"""
function TOV_quantities(u, EoS_rho_from_P)
    Pr = u[1, :]
    mr = u[2, :]
    rho_r = EoS_rho_from_P.(Pr)
    r_index = findfirst(x->x<0, Pr)
    M = mr[r_index]
    R = u.t[r_index]
    return Pr, mr, rho_r, M, R
end

# Polytropic EoS parameters
EoS_Param_Stiff = (
    K_crust = (3*π^2)^(2/3) * hbar^2 / (5*neutron_mass^(8/3)),
    gamma_crust = 5.0/3.0,
    gamma_core = 3.0, # Polytropic index for the core
    rho_b = 3e14*1e3, # Transition density between crust and core in kg/m^3
)
EoS_Param_Soft = (
    K_crust = (3*π^2)^(2/3) * hbar^2 / (5*neutron_mass^(8/3)),
    gamma_crust = 5.0/3.0,
    gamma_core = 5.0/2.0, # Polytropic index for the core
    rho_b = 3e14*1e3, # Transition density between crust and core in kg/m^3
)
# Simulation Input Parameters
Sim_Input = (
    ρ0 = 1e15*1e3, # Initial central density in kg/m^3 (1e15 g/cm^3)
    dr = 0.01*1e3, # Radial step in meters
    Dr = 0.01*1e3, # Radial interval for recording values in meters
    r_end = 20e3, # Maximum radius to solve up to in meters
);

# Function Setup for inverse EoS and TOV equation for the solver
EoS_P_from_rho_stiff, EoS_rho_from_P_stiff =
    QSpin.TOV.EoS_two_component_polytrope(EoS_Param_Stiff, true);
EoS_P_from_rho_soft, EoS_rho_from_P_soft =
    QSpin.TOV.EoS_two_component_polytrope(EoS_Param_Soft, true);

# Setting up initial condition accordingly to the EoS for a given central density ρ0.
u0_Stiff = [EoS_P_from_rho_stiff(Sim_Input.ρ0); 0.0]; # Initial conditions: central pressure and enclosed mass
u0_Soft = [EoS_P_from_rho_soft(Sim_Input.ρ0); 0.0]; # Initial conditions: central pressure and enclosed mass

# Solve with stiff parameters
tov_stiff = QSpin.TOV.tov_eq!(EoS_rho_from_P_stiff)
@time u_stiff = QSpin.OdeSolve.evolve(
    tov_stiff,
    u0_Stiff,
    0.0,
    Sim_Input.r_end;
    alg = DE.Tsit5(),
    dt = Sim_Input.dr,
    saveat = Sim_Input.Dr,
)
Pr_Stiff, mr_Stiff, ρr_Stiff, M_Stiff, R_Stiff =
    TOV_quantities(u_stiff, EoS_rho_from_P_stiff)
r = u_stiff.t

# Solve with soft parameters
tov_soft = QSpin.TOV.tov_eq!(EoS_rho_from_P_soft)
@time u_soft = QSpin.OdeSolve.evolve(
    tov_soft,
    u0_Soft,
    0.0,
    Sim_Input.r_end;
    alg = DE.Tsit5(),
    dt = Sim_Input.dr,
    saveat = Sim_Input.Dr,
)
Pr_Soft, mr_Soft, ρr_Soft, M_Soft, R_Soft = TOV_quantities(u_soft, EoS_rho_from_P_soft)

# Computing the M-R relation by varying the central density and solving the TOV equation for each case. The radius is defined by the first point that the pressure becomes negative.
ρc_scan = exp10.(range(17.5, stop = 20, length = 150))
M_StiffScan = zeros(length(ρc_scan));
R_StiffScan = zeros(length(ρc_scan));
M_SoftScan = zeros(length(ρc_scan));
R_SoftScan = zeros(length(ρc_scan));
for cc = 1:length(ρc_scan)
    ρc = ρc_scan[Int.(cc)];

    Pr, mr, ρr, M, R = QSpin.OdeSolve.TOV_Solve_rk4(
        [EoS_P_from_rho_stiff(ρc); 0.0],
        Sim_Input.dr,
        Sim_Input.Dr,
        Sim_Input.r_end,
        EoS_rho_from_P_stiff,
    );
    M_StiffScan[Int.(cc)] = M
    R_StiffScan[Int.(cc)] = R
    # println("Mass: ", M/sun_mass, " M_sun, Radius: ", R/1e3, " km")
    Pr, mr, ρr, M, R = QSpin.OdeSolve.TOV_Solve_rk4(
        [EoS_P_from_rho_soft(ρc); 0.0],
        Sim_Input.dr,
        Sim_Input.Dr,
        Sim_Input.r_end,
        EoS_rho_from_P_soft,
    );
    M_SoftScan[Int.(cc)] = M
    R_SoftScan[Int.(cc)] = R
end

# Plotting
Pc_Stiff = EoS_P_from_rho_stiff(Sim_Input.ρ0) # Getting the reference core pressure for the stiff case.
Pc_Soft = EoS_P_from_rho_soft(Sim_Input.ρ0) # Getting the reference core pressure for the soft case.
output_fig = plot(
    plot(
        r/1e3,
        [Pr_Stiff/Pc_Stiff Pr_Soft/Pc_Soft],
        label = [string(L"\gamma_\mathrm{core}=", EoS_Param_Stiff.gamma_core) string(
            L"\gamma_\mathrm{core}=",
            EoS_Param_Soft.gamma_core,
        )],
        xlabel = "Radius (km)",
        ylabel = L"P/P_c",
        framestyle = :box,
        linewidth = 2,
    ),
    plot(
        r/1e3,
        [mr_Stiff mr_Soft]/sun_mass,
        label = [string(L"\gamma_\mathrm{core}=", EoS_Param_Stiff.gamma_core) string(
            L"\gamma_\mathrm{core}=",
            EoS_Param_Soft.gamma_core,
        )],
        xlabel = "Radius (m)",
        ylabel = L"m/M_\odot",
        framestyle = :box,
        linewidth = 2,
    ),
    plot(
        r/1e3,
        [ρr_Stiff ρr_Soft]/1e18,
        label = [string(L"\gamma_\mathrm{core}=", EoS_Param_Stiff.gamma_core) string(
            L"\gamma_\mathrm{core}=",
            EoS_Param_Soft.gamma_core,
        )],
        xlabel = "Radius (m)",
        ylabel = L"\rho\;(\times10^{18}\;\textrm{kg/m}^3)",
        framestyle = :box,
        linewidth = 2,
    ),
    scatter(
        [R_StiffScan R_SoftScan]/1e3,
        [M_StiffScan M_SoftScan]/sun_mass;
        zcolor = log.([ρc_scan ρc_scan])/log(10), # color the data by the core density
        markershape = [:circle :diamond],
        markersize = [2 2],
        label = ["Stiff" "Soft"],
        xlabel = "Radius (km)",
        ylabel = L"\textrm{Mass}\;(M_\odot)",
        framestyle = :box,
        linewidth = 2,
    ),
    layout = (2, 2),
)
savefig(output_fig, "outputs/tov-output.png")
