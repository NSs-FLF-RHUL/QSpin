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

using ..OdeSolve: evolve
using ..Parameters: ParameterType
using ..PhysicalConstants: gravitational_constant, speed_of_light_vacuum
using OrdinaryDiffEqLowOrderRK: OrdinaryDiffEqLowOrderRK
using SciMLBase: DiscreteCallback, terminate!

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


function TOV_Solve(
    u0::Union{AbstractArray,Array{Float64}},
    dr::Float64,
    Dr::Float64,
    r_max::Float64,
    EoS_inv::Function;
    alg = OrdinaryDiffEqLowOrderRK.DP5(),
    dt = dr,
    saveat = Dr,
    reltol = 1e-8,
    solver_options...,
)
    tov! = tov_eq!(EoS_inv);
    condition(u, t, integrator) = u[1] < 0
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition, affect!)
    sol_tov = evolve(
        tov!,
        u0,
        0.0,
        r_max;
        alg,
        callback = cb,
    );
    ur = Array(sol_tov);
    r = sol_tov.t;
    Pr = ur[1, :];
    mr = ur[2, :];

    TOV_sol = (r = r, Pr = Pr, mr = mr, ρr = EoS_inv.(Pr), R = r[end], M = mr[end])
    return TOV_sol
end


end
