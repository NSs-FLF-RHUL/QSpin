using QSpin
using QSpin.Parameters: ParameterType
using Plots, LaTeXStrings
# This script is demonstrating how to solve the Tolman–Oppenheimer–Volkoff (TOV) equation for neutron stars using the QSpin package with the built-in TOV solver using Runge-Kutta 4th order method.
# The parameters are chosen from https://github.com/vanessagraber/teaching_materials/blob/master/summerschool_CRAQ_2019/mass_radius_relations.ipynb.


include("PhysConsts.jl")

# Polytropic EoS parameters
EoS_Param_Stiff = (
    γcore = 3.0, # Polytropic index for the core
    ρb = 3e14*1e3, # Transition density between crust and core in kg/m^3
)
EoS_Param_Soft = (
    γcore = 5.0/2.0, # Polytropic index for the core
    ρb = 3e14*1e3, # Transition density between crust and core in kg/m^3
)
# Simulation Input Parameters
Sim_Input = (
    ρ0 = 1e15*1e3, # Initial central density in kg/m^3 (1e15 g/cm^3)
    dr = 0.01*1e3, # Radial step in meters
    Dr = 0.01*1e3, # Radial interval for recording values in meters
    r_end = 20e3, # Maximum radius to solve up to in meters
);

"""
Set up the EoS functions for a given set of parameters and physical constants.

# Arguments
- `EoS_Param::ParameterType`: The parameters for the polytropic EoS.
- `PhysConst::ParameterType`: The physical constants in SI units.

# Returns
- `EoS_P::Function`: The EoS function that gives pressure as a function of density.
- `EoS_Inv::Function`: The inverse EoS function that gives density as a function of pressure.

# A polytropic EoS is used in this script.

"""
function EoS_Func(EoS_Param::ParameterType, PhysConst::ParameterType)
    Kcrust = (3*π^2)^(2/3) * PhysConst.ħ^2 / (5*PhysConst.mn^(8/3)); # Crust EoS constant
    γcrust = 5/3; # Crust EoS polytropic
    Kcore = Kcrust * EoS_Param.ρb^(γcrust-EoS_Param.γcore); # Core EoS constant to ensure continuity at ρb

    """
    Define the EoS function giving the P-ρ relation.

    # Arguments
    - `ρ::Union{Float64,AbstractArray,Array{Float64},Vector{Float64}}`: The density or array of densities for which to compute the pressure.

    # Returns
    - `P::Union{Float64,AbstractArray,Array{Float64},Vector{Float64}}`: The corresponding pressure or array of pressures computed from the EoS.
    """
    function EoS_P(ρ::Union{Float64,AbstractArray,Array{Float64},Vector{Float64}})
        if ρ < 0
            P = 0.0; # Ensure non-negative pressure
        elseif ρ <= EoS_Param.ρb
            P = Kcrust * ρ .^ γcrust;
        else
            P = Kcore * ρ .^ EoS_Param.γcore;
        end
        return P;
    end
    Pb = EoS_P(EoS_Param.ρb); # Pressure at the crust-core transition for continuity check
    println("Pressure at crust-core transition (Pb): ", Pb)

    """
    Define the inverse EoS function giving the ρ-P relation.
    # Arguments
    - `P::Union{Float64,AbstractArray,Array{Float64},Vector{Float64}}`: The pressure or array of pressures for which to compute the density.

    # Returns
    - `ρ::Union{Float64,AbstractArray,Array{Float64},Vector{Float64}}`: The corresponding density or array of densities computed from the EoS.
    """
    function EoS_Inv(P::Union{Float64,AbstractArray,Array{Float64},Vector{Float64}})
        if P < 0
            ρ = 0.0; # Ensure non-negative density
        elseif P <= Pb
            ρ = (P / Kcrust) .^ (1/γcrust);
        else
            ρ = (P / Kcore) .^ (1/EoS_Param.γcore);
        end
        return ρ
    end

    return EoS_P, EoS_Inv
end

# Fucntion Setup for inverse EoS and TOV equation for the solver
EoS_Stiff, EoS_inv_Stiff = EoS_Func(EoS_Param_Stiff, PhysConst);
EoS_Soft, EoS_inv_Soft = EoS_Func(EoS_Param_Soft, PhysConst);

# Setting up initial condition accordingly to the EoS for a given central density ρ0.
u0_Stiff = [EoS_Stiff(Sim_Input.ρ0); 0.0]; # Initial conditions: central pressure and enclosed mass
u0_Soft = [EoS_Soft(Sim_Input.ρ0); 0.0]; # Initial conditions: central pressure and enclosed mass

# Sovling the TOV equation using the RK4 method with QSpin OdeSolve module
@time Pr_Stiff, mr_Stiff, ρr_Stiff, M_Stiff, R_Stiff, r = QSpin.OdeSolve.TOV_Solve_rk4(
    u0_Stiff,
    Sim_Input.dr,
    Sim_Input.Dr,
    Sim_Input.r_end,
    EoS_inv_Stiff,
    PhysConst,
);

@time Pr_Soft, mr_Soft, ρr_Soft, M_Soft, R_Soft = QSpin.OdeSolve.TOV_Solve_rk4(
    u0_Soft,
    Sim_Input.dr,
    Sim_Input.Dr,
    Sim_Input.r_end,
    EoS_inv_Soft,
    PhysConst,
);

# Plotting
Pc_Stiff = EoS_Stiff(Sim_Input.ρ0) # Check continuity at crust-core transition
Pc_Soft = EoS_Soft(Sim_Input.ρ0) # Check continuity at crust-core transition

ρc_scan = exp10.(range(17.5, stop = 20, length = 150))


M_StiffScan = zeros(length(ρc_scan));
R_StiffScan = zeros(length(ρc_scan));
M_SoftScan = zeros(length(ρc_scan));
R_SoftScan = zeros(length(ρc_scan));
for cc = 1:length(ρc_scan)
    ρc = ρc_scan[Int.(cc)];

    Pr, mr, ρr, M, R = QSpin.OdeSolve.TOV_Solve_rk4(
        [EoS_Stiff(ρc); 0.0],
        Sim_Input.dr,
        Sim_Input.Dr,
        Sim_Input.r_end,
        EoS_inv_Stiff,
        PhysConst,
    );
    M_StiffScan[Int.(cc)] = M
    R_StiffScan[Int.(cc)] = R
    # println("Mass: ", M/PhysConst.Msun, " M_sun, Radius: ", R/1e3, " km")
    Pr, mr, ρr, M, R = QSpin.OdeSolve.TOV_Solve_rk4(
        [EoS_Soft(ρc); 0.0],
        Sim_Input.dr,
        Sim_Input.Dr,
        Sim_Input.r_end,
        EoS_inv_Soft,
        PhysConst,
    );
    M_SoftScan[Int.(cc)] = M
    R_SoftScan[Int.(cc)] = R
end

plot(
    plot(
        r/1e3,
        [Pr_Stiff/Pc_Stiff Pr_Soft/Pc_Soft],
        label = [string(L"\gamma_\mathrm{core}=", EoS_Param_Stiff.γcore) string(
            L"\gamma_\mathrm{core}=",
            EoS_Param_Soft.γcore,
        )],
        xlabel = "Radius (km)",
        ylabel = L"P/P_c",
        framestyle = :box,
        linewidth = 2,
    ),
    plot(
        r/1e3,
        [mr_Stiff mr_Soft]/PhysConst.Msun,
        label = [string(L"\gamma_\mathrm{core}=", EoS_Param_Stiff.γcore) string(
            L"\gamma_\mathrm{core}=",
            EoS_Param_Soft.γcore,
        )],
        xlabel = "Radius (m)",
        ylabel = L"m/M_\odot",
        framestyle = :box,
        linewidth = 2,
    ),
    plot(
        r/1e3,
        [ρr_Stiff ρr_Soft]/1e18,
        label = [string(L"\gamma_\mathrm{core}=", EoS_Param_Stiff.γcore) string(
            L"\gamma_\mathrm{core}=",
            EoS_Param_Soft.γcore,
        )],
        xlabel = "Radius (m)",
        ylabel = L"\rho\;(\times10^{18}\;\textrm{kg/m}^3)",
        framestyle = :box,
        linewidth = 2,
    ),
    scatter(
        [R_StiffScan R_SoftScan]/1e3,
        [M_StiffScan M_SoftScan]/PhysConst.Msun;
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
#plot!(size=(800,400))
