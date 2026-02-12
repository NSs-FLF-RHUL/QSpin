module OdeSolve
using LinearAlgebra
using FFTW
using MAT
using ParallelStencil


# Loading the rk4 propagator, ode_rk4.jl
include("ode_rk4.jl")

function evolve_rk4(ψ0::Array{Float64}, dt::Float64, Dt::Float64, tend::Float64, eom::Function)
    ΔNt = Int(Dt / dt)
    Nt = Int(tend / Dt)
    dims = ndims(ψ0)

    println("Field is ", dims, "D dimensional", " \n")

    ψall = zeros(Nt + 1, size(ψ0)...)
    selectdim(ψall, 1, 1) .= ψ0
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
            selectdim(ψall, 1, save_number + 1) .= ψcurrent
            tspan[save_number+1] = t
            save_number += 1
        end
    end
    return ψall, tspan
end

end
