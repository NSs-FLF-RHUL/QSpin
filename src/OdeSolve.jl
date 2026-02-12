module OdeSolve
using LinearAlgebra
using FFTW
using MAT
using ParallelStencil


# Loading the rk4 propagator, ode_rk4.jl
include("ode_rk4.jl")

function evolve_rk4(ψ0::Array{Float64}, dt::Float64, Dt::Float64, tend::Float64, eom::Function)
    # Method convention is to allocate time slices along a new dimension.
    dims = ndims(ψ0)
    time_dimension_index = dims + 1
    println("Field is ", dims, "D dimensional. Time slices will be along dimension ", time_dimension_index, "\n")

    ΔNt = Int(Dt / dt)
    Nt = Int(tend / Dt)

    ψall = zeros(size(ψ0)..., Nt + 1)
    selectdim(ψall, time_dimension_index, 1) .= ψ0
    tspan = zeros(Nt + 1) # i think these can be preallocated

    t = 0.
    ψcurrent = ψ0
    save_number = 1
    step_number = 0

    while t < tend
        ψcurrent = ode_rk4(ψcurrent, dt, t, eom)
        t += dt
        step_number += 1

        if mod(step_number, ΔNt) == 0
            println("t=", t, " \n")
            selectdim(ψall, time_dimension_index, save_number + 1) .= ψcurrent
            tspan[save_number+1] = t
            save_number += 1
        end
    end
    return ψall, tspan
end

end
