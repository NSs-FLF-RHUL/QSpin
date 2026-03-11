module Hamiltonian
"""
hamiltonian(KineticEnergy,PotentialEnergy,irt)(ψ,time)

This function returns a function H(ψ,time) that computes the Hamiltonian.

:param KineticEnergy: a function that computes the kinetic energy term in the equation of motion, which takes ψ and time as input.
:param PotentialEnergy: a function that computes the potential energy term in the equation of motion, which takes ψ and time as input.
:param irt: the propagation factor, which can be -1 for imaginary time evolution, -im for real time evolution, or a general complex value for a dissipative evolution, namely, dGPE.
:param ψ: the field variable for computing the Hamiltonian, which can be a vector or an array depending on the problem.
:param time: the time variable

"""
    function hamiltonian(KineticEnergy,PotentialEnergy,irt)
        """
        internal_hamiltonian(ψ::Union{AbstractArray,Array{Float64},Array{ComplexF64}},time::Float64)

        Here evalute ψ and time based on the Hamiltonian defined by the input kinetic energy and potential energy functions, as well as the propagation factor irt.

        """
        function internal_hamiltonian(ψ::Union{AbstractArray,Array{Float64},Array{ComplexF64}},time::Float64)
            return irt .* (KineticEnergy(ψ,time) + PotentialEnergy(ψ,time))
        end
        return internal_hamiltonian
    end
end