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

# Glitch Modol module
include("GlitchModels.jl")
end
