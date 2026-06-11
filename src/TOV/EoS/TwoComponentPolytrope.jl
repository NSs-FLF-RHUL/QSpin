using DocStringExtensions: TYPEDSIGNATURES
using Roots: find_zero
using ....Parameters: ParameterType

"""
$(TYPEDSIGNATURES)

Create EoS functions for two-component polytrope neutron stars.

Such stars are assumed to be spheres composed of non-interacting, degenerate neutrons that are characterised by two EoS, in the core and crust regions.
The boundary between the regions is located at some density ``\\rho_b``;

```math
P = \\begin{cases}
K_\\mathrm{crust} \\rho^{\\gamma_{\\mathrm{crust}}} & \\rho < \\rho_b, \\\\
K_\\mathrm{core} \\rho^{\\gamma_{\\mathrm{core}}} & \\rho > \\rho_b,
\\end{cases}
\\quad
K_\\mathrm{core} = K_\\mathrm{crust} \\rho_b^{ \\gamma_{\\mathrm{crust}} - \\gamma_{\\mathrm{core}} },
```

with the equation for ``K_\\mathrm{core}`` being a result of imposing pressure continuity at the crust-core interface.

Assuming the crust is pure degenerate neutron matter, the non-relativistic EoS is given by

``
P_\\mathrm{crust} = \\frac{(3\\pi^2)^{2/3}}{5} \\frac{\\bar{h}^2}{m_n^{8/3}} \\rho^{5/3}.
``

namely,

``
K_\\mathrm{crust} = \\frac{(3\\pi^2)^{2/3}}{5} \\frac{\\bar{h}^2}{m_n^{8/3}}.
``

To keep the generality of the function, we keep ``K_\\mathrm{crust}`` as an input parameter here.

`parameters` is assumed to define the following quantities:

- ``K_\\mathrm{crust}``, under the name `K_crust`.
- ``\\gamma_\\mathrm{crust}``, under the name `gamma_crust`.
- ``\\gamma_\\mathrm{core}``, under the name `gamma_core`.
- ``\\rho_b``, under the name `rho_b`.

# Arguments
- `parameters::ParameterType`: Parameter values for the two-component polytrope model.
- `report_transition_pressure::Bool`: If `true`, the transition pressure at the interface will be printed.

# Returns
- `EoS_P_from_rho::Function`: ``P(\\rho)`` as a function called with `(rho,)`.
- `EoS_rho_from_P::Function`: ``\\rho(P)`` as a function called with `(P,)`.
"""
function EoS_two_component_polytrope(
    ParamIn::ParameterType,
    report_transition_pressure::Bool = false,
)
    K_core = ParamIn.K_crust * ParamIn.ρ_b ^ (ParamIn.γ_crust - ParamIn.γ_core)

    function EoS_P_from_rho(ρ)
        P = if ρ < 0
            0.0
        elseif ρ <= ParamIn.ρ_b
            ParamIn.K_crust * ρ .^ ParamIn.γ_crust
        else
            K_core * ρ .^ ParamIn.γ_core
        end
        return P
    end

    # Pressure at the crust-core transition for continuity check
    P_b = EoS_P_from_rho(ParamIn.ρ_b)

    if report_transition_pressure
        println("Pressure at crust-core transition (Pb): ", P_b)
    end

    function EoS_rho_from_P(P)
        ρ = if P < 0
            0.0
        elseif P <= P_b
            (P / ParamIn.K_crust) .^ (1/ParamIn.γ_crust)
        else
            (P / K_core) .^ (1/ParamIn.γ_core)
        end
        return ρ
    end

    return EoS_P_from_rho, EoS_rho_from_P
end
