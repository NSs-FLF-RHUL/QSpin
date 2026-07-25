"""
Submodule containing specific Equations of State that may be used in the TOV equations.
"""
module EquationOfState

include("GraberCummingAnderson2018.jl")
include("NegeleVautherin1973.jl")
include("TwoComponentPolytrope.jl")
include("EoS_LInterp.jl")
end
