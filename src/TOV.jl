"""
    Sumbodule containing the helper fucntion for solvign the Tolman–Oppenheimer–Volkoff (TOV) equation, taht describes the equation of state for the neutron star interior under the non-relativistic approximation limit.

    The TOV equation is a ordinary differential equation of pressure as a function of radius r and is written by

    ``
    \\frac{\\mathrm{d}P}{\\mathrm{d}r} = - \\frac{G}{r^2}
    \\left[ \\rho(r) + \\frac{P(r)}{c^2} \\right]
    \\left[ m(r) + 4\\pi r^3 \\frac{P(r)}{c^2} \\right]
    \\left[ 1 - \\frac{2Gm(r)}{c^2 r} \\right]^{-1},
    ``

    where ``G`` and ``c`` represent the gravitational constant and speed of light in vacuum, respectively.

    It is coupled to the mass m(r) within radius r, which is given by

    ``
    \\frac{\\mathrm{d}m}{\\mathrm{d}r} = 4\\pi\\rho(r)r^2.
    ``

    where ρ is the density as a function of r and is assumeed to be spherically symmetric.

    To solve the TOV equation, we need to specify an equation of state (EoS) that relates the pressure P to the density ρ, that is we must provide ``P = P(\\rho)`` (and the inverse). Here, to solve the TOV equation numerically, we use the vector presentation for the variables, `u = (P, m)`` and solving a first-order system ODE in ``u``.
"""
module TOV

#using Dierckx
using ..OdeSolve: ode_rk4
using ..OdeSolve: DESolve
using ..Parameters: ParameterType
using ..PhysicalConstants: hbar, gravitational_constant, neutron_mass, speed_of_light_vacuum

function tov_eq!(EoS_rho_from_P::Function)
    function tov_inner!(du, u, paras, r)
        P = u[1];
        m = u[2];
        ρ = EoS_rho_from_P(P);
        if r == 0.0
            du[1] = 0.0;
        else
            du[1] =
                -gravitational_constant / r^2 *
                (ρ + P/speed_of_light_vacuum^2) *
                (m + 4*pi*r^3*P/speed_of_light_vacuum^2) /
                (1 - 2*gravitational_constant*m/(r*speed_of_light_vacuum^2));
        end
        du[2] = 4*pi*r^2*ρ;
    end
    return tov_inner!
end


@doc raw"""

    Create EoS functions for two-component polytrope neutron stars.

    Such stars are assumed to be spheres composed of non-interacting, degenerate neutrons that are characterised by two EoS, in the core and crust regions.
    The boundary between the regions is located at some density ``\rho_b``;

    ```math
    P = \begin{cases}
    K_\mathrm{crust} \rho^{\gamma_{\mathrm{crust}}} & \rho < \rho_b, \\
    K_\mathrm{core} \rho^{\gamma_{\mathrm{core}}} & \rho > \rho_b,
    \end{cases}
    \quad
    K_\mathrm{core} = K_\mathrm{crust} \rho_b^{ \gamma_{\mathrm{crust}} - \gamma_{\mathrm{core}} },
    ```

    with the equation for ``K_\mathrm{core}`` being a result of imposing pressure continuity at the crust-core interface.

    Assuming the crust is pure degenerate neutron matter, the non-relativistic EoS is given by

    ``
    P_\mathrm{crust} = \frac{(3\pi^2)^{2/3}}{5} \frac{\bar{h}^2}{m_n^{8/3}} \rho^{5/3}.
    ``

    namely,

    ``
    K_\mathrm{crust} = \frac{(3\pi^2)^{2/3}}{5} \frac{\bar{h}^2}{m_n^{8/3}}.
    ``

    To keep the generality of the function, we keep ``K_\mathrm{crust}`` as an input parameter here.

    `parameters` is assumed to define the following quantities:

    - ``K_\mathrm{crust}``, under the name `K_crust`.
    - ``\gamma_\mathrm{crust}``, under the name `gamma_crust`.
    - ``\gamma_\mathrm{core}``, under the name `gamma_core`.
    - ``\rho_b``, under the name `rho_b`.

    # Arguments
    - `parameters::ParameterType`: Parameter values for the two-component polytrope model.
    - `report_transition_pressure::Bool`: If `true`, the transition pressure at the interface will be printed.

    # Returns
    - `EoS_P_from_rho::Function`: ``P(\rho)`` as a function called with `(rho,)`.
    - `EoS_rho_from_P::Function`: ``\rho(P)`` as a function called with `(P,)`.
"""
function EoS_two_component_polytrope(
    ParamIn::ParameterType,
    report_transition_pressure::Bool = false,
)
    K_core = ParamIn.K_crust * ParamIn.ρ_b ^ (ParamIn.γ_crust - ParamIn.γ_core);

    function EoS_P_from_rho(ρ)
        if ρ < 0
            P = 0.0;
        elseif ρ <= ParamIn.ρ_b
            P = ParamIn.K_crust * ρ .^ ParamIn.γ_crust;
        else
            P = K_core * ρ .^ ParamIn.γ_core;
        end
        return P;
    end

    # Pressure at the crust-core transition for continuity check
    P_b = EoS_P_from_rho(ParamIn.ρ_b);

    if report_transition_pressure
        println("Pressure at crust-core transition (Pb): ", P_b)
    end

    function EoS_rho_from_P(P)
        if P < 0
            ρ = 0.0;
        elseif P <= P_b
            ρ = (P / ParamIn.K_crust) .^ (1/ParamIn.γ_crust);
        else
            ρ = (P / K_core) .^ (1/ParamIn.γ_core);
        end
        return ρ
    end

    return EoS_P_from_rho, EoS_rho_from_P
end



function TOV_Solve_rk4(
    u0::Union{AbstractArray,Array{Float64}},
    dr::Float64,
    Dr::Float64,
    r_end::Float64,
    EoS_inv::Function,
)
    # Physical constants in SI units

    """
    Set up the TOV equation for a given inverse EoS function and physical constants.

    # Arguments
    - `EoS_inv::Function`: The inverse EoS function that gives density as
    - 'u::AbstractArray': The current state of the system, where u[1] is pressure P and u[2] is enclosed mass m.
    -'r::Float64': The current radius at which the TOV equation is being evaluated.

    # Returns
    - `tov_eq::Function`: A function of [dP/dr; dm/dr] for the TOV equation.

    """
    function TOV_Eq(Eos_inv::Function)
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
        function tov_eq(u::AbstractArray, r::Float64)
            P = u[1];
            m = u[2];
            ρ = EoS_inv(P);
            if r == 0.0
                dPdr = 0;
            else
                dPdr =
                    -PhysConst.G / r^2 *
                    (ρ + P/PhysConst.c^2) *
                    (m + 4*π*r^3*P/PhysConst.c^2) / (1 - 2*PhysConst.G*m/(r*PhysConst.c^2));
            end
            dmdr = 4*π*r^2*ρ;
            return [dPdr; dmdr]
        end
        return tov_eq
    end
    TOV = TOV_Eq(EoS_inv);
    ΔNr = floor(Int, Dr / dr)
    Nr = floor(Int, r_end / Dr)
    uall = zeros(eltype(u0), size(u0)..., Nr + 1)
    dims = ndims(u0)
    radial_dimension_index = dims + 1
    selectdim(uall, radial_dimension_index, 1) .= u0
    rspan = zeros(Nr + 1)

    r = 0.0
    ucurrent = u0
    save_number = 1
    step_number = 0

    @inbounds while r < r_end
        ucurrent = ode_rk4(ucurrent, dr, r, TOV)
        r += dr
        step_number += 1
        if sum(isnan.(ucurrent[:]))>0
            println(
                "NaN detected in the field at radius ",
                r,
                ". Radius Step could be too big.",
            )
            break
        end
        if mod(step_number, ΔNr) == 0
            selectdim(uall, radial_dimension_index, save_number + 1) .= ucurrent
            rspan[save_number+1] = r
            save_number += 1
        end
    end

    Pr = uall[1, :];
    mr = uall[2, :];
    ρr = zeros(length(rspan))
    for rr = 1:length(rspan)
        ρr[rr] = EoS_inv(Pr[rr]) # Density profile from the inverse EoS
    end

    R_index = findfirst(x->x<0, Pr);
    M = mr[R_index];
    R = rspan[R_index];
    return Pr, mr, ρr, M, R, rspan

end


function MutualFriction(ParamIn::ParameterType)

    if typeof(ParamIn.B_core) == Float64
        B_core = ParamIn.B_core;
    elseif typeof(ParamIn.B_core) == Matrix{Float64}
        println("B_core Interpolation")
        spl = Spline1D(ParamIn.ρr_core, ParamIn.B_core, k = 3, bc = "extrapolate")
        B_core = spl(ParamIn.ρr);
    else
        tyopeof(ParamIn.B_core)
        println("B_core is a function of density")
        B_core = ParamIn.B_core(ParamIn.ρr);
    end

    if typeof(ParamIn.B_sf) == Float64
        B_sf = ParamIn.B_sf;
    elseif typeof(ParamIn.B_sf) == Matrix{Float64}
        println("B_sf Interpolation")
        spl = Spline1D(ParamIn.ρr_sf, ParamIn.B_sf, k = 3, bc = "extrapolate")
        B_sf = spl(ParamIn.ρr);
    else
        tyopeof(ParamIn.B_sf)
        println("B_sf is a function of density")
        B_sf = ParamIn.B_sf(ParamIn.ρr);
    end

    B_sf[B_sf .< ParamIn.ρ_dripping] .= 0;
    B_core[B_core .< ParamIn.ρ_dripping] .= 0;
    return B_crust, B_core, B_sf
end

end
