module OdeSolve
    using LinearAlgebra
    using FFTW
    using MAT
    using ParallelStencil


# Loading the eom.jl and rk4 propagator, ode_rk4.jl
    include("ode_rk4.jl")

    function evolve1d(ψ0::Vector{Float64},dt::Float64,Dt::Float64,tend::Float64,eom::Function)
        save_number = 1
        step_number = 0
        ΔNt         = Int.(Dt / dt)
        Nt          = Int.(tend / Dt)
        ψall        = zeros(length(ψ0),Nt+1)
        tspan       = zeros(Nt+1)
        ψall[:,1]   = ψ0
        t           = 0.
        while t < tend
            ψn = ode_rk4(ψ0,dt,t,eom)
            ψ0 = ψn
            t  += dt
            step_number += 1
            if mod(step_number,ΔNt) == 0
                println("t=",t, " \n")
                ψall[:,save_number+1] = ψ0
                tspan[save_number+1] = t
                save_number += 1
            end
        end
        return ψall, tspan
    end

    function evolve2d(ψ0::Vector{Float64},dt::Float64,Dt::Float64,tend::Float64)
        save_number = 1
        step_number = 0
        ΔNt         = Int.(Dt / dt)
        Nt          = Int.(tend / Dt)
        ψall        = zeros(length(ψ0),Nt+1)
        tspan       = zeros(Nt+1)
        ψall[:,1]   = ψ0
        t           = 0.
        while t < tend
            ψn = ode_rk4(ψ0,dt,t)
            ψ0 = ψn
            t  += dt
            step_number += 1
            if mod(step_number,ΔNt) == 0
                println("t=",t, " \n")

                tspan[save_number+1] = t
                save_number += 1
            end
        end
        return ψall, tspan
    end

end
