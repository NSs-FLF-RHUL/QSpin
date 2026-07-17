module QSpin

# Load ParameterType for solvers
include("Parameters.jl")

# Load constants and other static variables
include("PhysicalConstants.jl")

include("OdeSolve.jl")
include("Grids.jl")
include("Hamiltonian.jl")

# TOV module
include("TOV/TOV.jl")

# Equation of State module
include("TOV/EoS/EquationOfState.jl")
end
