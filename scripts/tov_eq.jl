using QSpin
using QSpin.Parameters: ParameterType
using Plots, LaTeXStrings

# Physical constants in SI units
PhysConst = (
    ħ = 1.0545718 * 1e-34, # m^2*kg / s
    Msun = 1.9891 * 1e30,      # kg
    c = 299792458,         # m / s
    G = 6.67408 * 1e-11,   # m^3 / (kg * s^2)
    kpc = 3.08567758 * 1e19, # m
    eV = 1.782662 * 1e-36,  # kg
    Gyear = 31556926 * 1e9,   # s
    mn = 1.674927471 * 1e-27, # kg
)

# Polytropic EoS parameters
EoS_Param = (
    γcore = 5/2.0, # Polytropic index for the core
    ρb = 3e14*1e3, # Transition density between crust and core in kg/m^3
)

# Simulation Input Parameters
Sim_Input = (
    ρ0 = 1e15*1e3, # Initial central density in kg/m^3 !! CHECKING UNITS
    dr = 0.01*1e3, # Radial step in meters
    Dr = 0.1*1e3, # Radial interval for recording values in meters
    r_end = 15e3, # Maximum radius to solve up to in meters
);
# Polytropic EoS for P-ρ relation
function EoS_Func(EoS_Param::ParameterType, PhysConst::ParameterType)
    Kcrust = (3*π^2)^(2/3) * PhysConst.ħ^2 / (5*PhysConst.mn^(8/3)); # Crust EoS constant
    γcrust = 5/3; # Crust EoS polytropic
    Kcore = Kcrust * EoS_Param.ρb^(γcrust-EoS_Param.γcore); # Core EoS constant to ensure continuity at ρb
    function EoS_P(ρ::Float64)
        if ρ < 0
            P = 0.0; # Ensure non-negative pressure
        elseif ρ <= EoS_Param.ρb
            P = Kcrust * ρ^γcrust;
        else
            P = Kcore * ρ^EoS_Param.γcore;
        end
        return P;
    end
    Pb = EoS_P(EoS_Param.ρb); # Pressure at the crust-core transition for continuity check
    println("Pressure at crust-core transition (Pb): ", Pb)
    function EoS_Inv(P::Float64)
        if P < 0
            ρ = 0.0; # Ensure non-negative density
        elseif P <= Pb
            ρ = (P / Kcrust)^(1/γcrust);
        else
            ρ = (P / Kcore)^(1/EoS_Param.γcore);
        end
        return ρ
    end

    return EoS_P, EoS_Inv
end

# ToV equation setup
function TOV_Eq(Eos_inv::Function, PhysConst::ParameterType)
    function tov_eq(u::AbstractArray, r::Float64)
        P = u[1];
        m = u[2];
        ρ = EoS_inv(P);
        if r == 0.0
            dPdr = 0;
        else
            dPdr =
                -PhysConst.G / r^2 * (ρ + P/PhysConst.c^2) * (m + 4*π*r^3*P/PhysConst.c^2) /
                (1 - 2*PhysConst.G*m/(r*PhysConst.c^2));
        end
        dmdr = 4*π*r^2*ρ;
        return [dPdr; dmdr]
    end
    return tov_eq
end

# Fucntion Setup for inverse EoS and TOV equation for the solver
EoS, EoS_inv = EoS_Func(EoS_Param, PhysConst);
TOV = TOV_Eq(EoS_inv, PhysConst);

# Setting up initial condition accordingly to the EoS for a given central density ρ0.
P0 = EoS(Sim_Input.ρ0); # Initial central pressure from EoS
m0 = 0.0; # Initial enclosed mass at the center
u0 = [P0; m0]; # Initial conditions: central pressure and enclosed mass

# Sovling the TOV equation using the RK4 method with QSpin OdeSolve module
@time ur, r =
    QSpin.OdeSolve.evolve_rk4(u0, Sim_Input.dr, Sim_Input.Dr, Sim_Input.r_end, TOV);

# Plotting
Pc = EoS(Sim_Input.ρ0) # Check continuity at crust-core transition
plot!(
    r/1e3,
    ur[1, :]/Pc,
    label = string(L"\gamma_{core}=", EoS_Param.γcore),
    xlabel = "Radius (km)",
    ylabel = L"P/P_c",
    title = string("TOV Solution for ", L"ρ_{c}=", Sim_Input.ρ0, L"\textrm{kg/m}^3"),
    xticks = range(0.0, stop = 15.0, length = 7),
    yticks = range(0.0, stop = 1.0, length = 6),
)
