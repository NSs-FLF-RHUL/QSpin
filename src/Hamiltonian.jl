module Hamiltonian

"""
Construct a Hamiltonian from a Kinetic and Potential energy.

Returns a function `H!(dψ, ψ, parameters, time)` that evaluates the Hamiltonian ``\\mathcal{H}``
formed from the given kinetic (``K``) and potential (``V``) energy functions at the given time:

``
\\mathcal{H}(\\psi, t) = \\mathrm{irt} \\times \\left( K(\\psi, t) + V(\\psi, t) \\right),
``

where ``\\mathrm{irt}`` is a given propagation factor.

The `KineticEnergy` and `PotentialEnergy` functions are assumed to take `(ψ, parameters, time)`
as positional arguments, and return the value of the respective energy quantities.

# Arguments
- `KineticEnergy::Function`: Function that computes the kinetic energy term in the equation of motion.
- `PotentialEnergy::Function`: Function that computes the potential energy term in the equation of motion.
- `irt::Complex64`: Propagation factor, which can be -1 for imaginary time evolution, -im for real time evolution,
    or a general complex value for a dissipative evolution, namely, dGPE.

# Returns
- `H!::Function`: The Hamiltonian formed from the kinetic and potential energies
"""
function hamiltonian!(KineticEnergy, PotentialEnergy, irt = im)
    function H!(dψ::AbstractArray, ψ::AbstractArray, parameters::NamedTuple, time::Float64)
        dψ .=
            irt .*
            (KineticEnergy(ψ, parameters, time) + PotentialEnergy(ψ, parameters, time))
    end
    return H!
end

end
