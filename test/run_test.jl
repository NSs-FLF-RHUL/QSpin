using LinearAlgebra
using FFTW
using Plots
using MAT
using ParallelStencil

# Parameter for the eom and initial condition setup
    a  = 0.05
    ψ0 = 0.1


include("../src/eom.jl")
include("../src/ode_rk4.jl")

function evolve(ψ0::Float64,dt::Float64,Dt::Float64,tend::Float64)
    save_number = 1
    step_number = 0
    ΔNt         = Int.(Dt / dt)
    Nt          = Int.(tend / Dt)
    ψall        = zeros(Nt+1)
    tspan       = zeros(Nt+1)
    ψall[1]   = ψ0
    t           = 0.
    while t < tend
        ψn = ode_rk4(ψ0,dt,t)
        ψ0 = ψn
        t  += dt
        step_number += 1
        if mod(step_number,ΔNt) == 0
            println("t=",t, " \n")
            ψall[save_number+1] = ψ0
            tspan[save_number+1] = t
            save_number += 1
        end
    end
    return ψall, tspan
end



