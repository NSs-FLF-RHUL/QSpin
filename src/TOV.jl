"""
Submodule containing helper functions for solving the TOV equation.

The TOV equation reads

``
\\frac{\\mathrm{d}P}{\\mathrm{d}r} &=
- \\frac{G}{r^2}
\\left[ \\rho(r) + \\frac{P(r)}{c^2} \\right]
\\left[ m(r) + 4\\pi r^3 \frac{P(r)}{c^2} \\right]
\\left[ 1 - \frac{2Gm(r)}{c^2 r} \\right]^{-1},
``

where ``G`` and ``c`` represent the gravitational constant and speed of light in vacuum, respectively.
``r`` is the radial coordinate, with ``\\rho(r)`` and ``P(r)`` the mass density and pressure as position ``r``, respectively.
Finally, ``m(r)`` is the gravitational mass satisfying ``m(0) = 0`` and is found by solving the continuity equation

``
\\frac{\\mathrm{d}m}{\\mathrm{d}r} &= 4\\pi\\rho(r)r^2.
``

To solve the TOV equation we must also provide an equation of state (EoS) that relates the pressure and density of the neutron star, that is we must provide ``P = P(\\rho)`` (and the inverse).

Note that we can solve the TOV equation(s) numerically using the vector ``u = (P, m)`` and solving a first-order system in ``u``.
"""

using ..Parameters: ParameterType
using ..PhysicalConstants: hbar, gravitational_constant, neutron_mass, speed_of_light_vacuum

"""
Return a function that evaluates the TOV equation.

Returns a function `tov!(du, u, params, r)` that evaluates the TOV equation (in vector form) given the EoS relationship ``\\rho(P)``.

`EoS_rho_to_P` is assumed to take `(P, parameters, r)` as inputs, where `P = u[1]` is the pressure.
"""
function tov_eq!(EoS_rho_to_P::Function)
    function tov_inner!(du, u, params, r)
        P = u[1]
        m = u[2]
        rho = EoS_rho_to_P(P);
        if r == 0.0
            du[1] .= 0.0;
        else
            du[1] .=
                - gravitational_constant / r^2 *
                (rho + P / speed_of_light_vacuum^2) *
                (m + 4*π*r^3*P / speed_of_light_vacuum^2) /
                (1 - 2*gravitational_constant*m / (r*speed_of_light_vacuum^2));
        end
        du[2] .= 4*π*r^2*rho;
    end
    return tov_inner!
end

"""
Create EoS functions for two-component polytrope neutron stars.

Such stars are assumed to be spheres composed of non-interacting, degenerate neutrons that are characterised by two EoS, in the core and crust regions.
The boundary between the regions is located at some radial coordinate ``r = \\rho_b``;

``
P =
\\begin{cases}
K_{crust} \\rho^{\\gamma_{crust}} & \\rho < \\rho_b, \\
K_{core} \\rho^{\\gamma_{core}} & \\rho > \\rho_b,
\\end{cases}
&\\qquad
K_{core} = K_{crust} \\rho_b^{\\gamma_{crust}-\\gamma_{core}},
``

with the equation for ``K_{core}`` being a result of imposing pressure continuity at the crust-core interface.

We further assume that the crust is pure degenerate neutron matter;

``
P_{crust} = \\frac{(3\\pi^2)^{2/3}}{5} \\frac{\\bar{h}^2}{m_n^{8/3}} \\rho^{5/3}.
``

`parameters` is assumed to define the following quantities:

- ``K_{crust}``, under the name ``K_crust``.
- ``\\gamma_{crust}``, under the name `gamma_crust`.
- ``\\gamma_{core}``, under the name `gamma_core`.
- ``\rho_b``, under the name `rho_b`.

# Arguments
- `parameters::ParameterType`: Parameter values for the two-component polytrope model.
- `report_transition_pressure::Bool`: If `true`, the transition pressure at the interface will be printed.

# Returns
- `EoS_P_to_rho::Function`: ``P(\\rho)`` as a function called with `(rho, params, r)`.
- `EoS_rho_to_P::Function`: ``\\rho(P)`` as a function called with `(P, params, r)`.
"""
function EoS_two_component_polytrope(
    parameters::ParameterType,
    report_transition_pressure::Bool = false,
)
    K_core =
        parameters.K_crust *
        parameters.rho_b ^ (parameters.gamma_crust - parameters.gamma_core);

    function EoS_P_to_rho(rho, params, r)
        if rho < 0
            P = 0.0;
        elseif rho <= parameters.rho_b
            P = parameters.K_crust * rho .^ parameters.gamma_crust;
        else
            P = K_core * rho .^ parameters.gamma_core;
        end
        return P;
    end

    # Pressure at the crust-core transition for continuity check
    P_b = EoS_P_to_rho(parameters.rho_b, (), 0.0);

    if report_transition_pressure
        println("Pressure at crust-core transition (Pb): ", P_b)
    end

    function EoS_rho_to_P(P, params, r)
        if P < 0
            rho = 0.0;
        elseif P <= P_b
            rho = (P / parameters.K_crust) .^ (1/parameters.gamma_crust);
        else
            rho = (P / K_core) .^ (1/parameters.gamma_core);
        end
        return rho
    end

    return EoS_P_to_rho, EoS_rho_to_P
end
