module QSpin

# Load ParameterType for solvers
include("Parameters.jl")

# Load physical constants that have no dependencies
include("PhysicalConstants.jl")

include("OdeSolve.jl")
include("Grids.jl")

# Load submodules that revolve around particular physical systems
include("Hamiltonian.jl")
include("TOV.jl")

end
