module Hamiltonian

"""
    hamiltonian!(KineticEnergy, PotentialEnergy, irt = im)

Return a function H!(dψ, ψ, parameters, time) that evaluates the Hamiltonian formed from the
given Kinetic and Potential energy functions, at the given time. The ``KineticEnergy`` and
``PotentialEnergy`` functions are assumed to take ``(ψ, parameters, time)`` as positional arguments.

:param KineticEnergy: a function that computes the kinetic energy term in the equation of motion, which takes ψ and time as input.
:param PotentialEnergy: a function that computes the potential energy term in the equation of motion, which takes ψ and time as input.
:param irt: the propagation factor, which can be -1 for imaginary time evolution, -im for real time evolution, or a general complex value for a dissipative evolution, namely, dGPE.
"""
function hamiltonian!(KineticEnergy, PotentialEnergy, irt = im)
    """
        internal_hamiltonian!(dψ::AbstractArray, ψ::AbstractArray, parameters, time::Float64)

    Evaluate the Hamiltonian at the given time and field ψ in-place, overwriting the input array dψ.

    :param dψ: Array to update in-place.
    :param ψ: Current field values.
    :param parameters: Container holding additional parameters needed by the equations of motion.
    :param time: Current time.
    """
    function internal_hamiltonian!(
        dψ::AbstractArray,
        ψ::AbstractArray,
        parameters::NamedTuple,
        time::Float64,
    )
        dψ .=
            irt .*
            (KineticEnergy(ψ, parameters, time) + PotentialEnergy(ψ, parameters, time))
    end
    return internal_hamiltonian!
end

end
