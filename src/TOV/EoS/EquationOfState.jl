"""
Submodule containing specific Equations of State that may be used in the TOV equations.
"""
module EquationOfState

using DocStringExtensions: TYPEDSIGNATURES

include("GraberCummingAnderson2018.jl")
include("NegeleVautherin1973.jl")
include("TwoComponentPolytrope.jl")

"""

$(TYPEDSIGNATURES)

Non-dimensionalise functions returned from an equation of state, according to the characteristic lengths provided.

Characteristic lengths should be provided in a `NamedTuple` of the same format as the output of `QSpin.TOV.characteristic_lengths_TOV`.

`builder` functions (like those provided in `QSpin.TOV.EquationOfState`) return two functions; one that evaluates pressure ``P`` as a function of density ``\\rho`` and another that does the inverse (in that order). This function mutates these outputs to instead by functions of the corresponding dimensionless pressure and dimensionless density. For example, if you were previously doing

```julia
EoS_P_from_rho, EoS_rho_from_P = builder(args...; kwargs...)
```

and want to turn `EoS_P_from_rho` and `EoS_rho_from_P` into their non-dimensional equivalents, using characteristic lengths with ``R = 1`` say, you would run

```julia
using QSpin.TOV: characteristic_lengths_TOV

char_lengths = characteristic_lengths_TOV(; length = 1.0)
nd_EoS_P_from_rho, nd_EoS_rho_from_P = nondimensional_EoS(char_lengths, builder, args...; kwargs...)
```

# Arguments
- `characteristic_lengths::NamedTuple`: Mapping of characteristic length names to values.
- `builder::Function`: The equation of state builder function, that operates on dimension-full fields. Likely a member of `QSpin.TOV.EquationOfState`.
- `args`: Positional arguments passed to `builder`.
- `kwargs`: Keyword arguments passed to `builder`.

# Returns
- `EoS_P_from_rho::Function`: Dimensionless pressure, evaluated as a function of dimensionless density. ``\\hat{P}(\\hat{\\rho})`` as a function called with `(rho,)`.
- `EoS_rho_from_P::Function`: Dimensionless density, evaluated as a function of dimensionless pressure. ``\\hat{\\rho}(\\hat{P})`` as a function called with `(P,)`
"""
function nondimensional_EoS(
    characteristic_lengths::NamedTuple,
    builder::Function,
    args...;
    kwargs...,
)
    EoS_P_from_rho, EoS_rho_from_P = builder(args...; kwargs...)

    Q = characteristic_lengths.Q
    Rho = characteristic_lengths.Rho

    nondim_Eos_P_from_rho = rho -> EoS_P_from_rho(rho * Rho) / Q
    nondim_EoS_rho_from_P = P -> EoS_rho_from_P(P * Q) / Rho

    return nondim_Eos_P_from_rho, nondim_EoS_rho_from_P
end

end
