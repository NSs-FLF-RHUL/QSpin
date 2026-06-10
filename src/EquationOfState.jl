"""
Submodule containing specific Equations of State that are used in the TOV equation.
"""
module EquationOfState

include("_EoS/GraberCummingAnderson2018.jl")
include("_EoS/NegeleVautherin1973.jl")
include("_EoS/TwoComponentPolytrope.jl")

end
