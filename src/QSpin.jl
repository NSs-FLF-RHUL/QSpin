module QSpin

# Load ParameterType for solvers
include("Parameters.jl")

include("OdeSolve.jl")
include("Grids.jl")
include("Hamiltonian.jl")

end
