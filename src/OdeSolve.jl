module OdeSolve
    using LinearAlgebra
    using FFTW
    using MAT
    using ParallelStencil


# Loading the eom.jl and rk4 propagator, ode_rk4.jl
    include("ode_rk4.jl")

    function evolve_rk4(ψ0::Array{Float64},dt::Float64,Dt::Float64,tend::Float64,eom::Function)
        save_number = 1
        step_number = 0
        ΔNt         = Int(Dt / dt)
        Nt          = Int(tend / Dt)
        dims        = ndims(ψ0)
        
        if dims <= 1
            dims == 1
            println("Effective ", dims,"D problem", " \n")
            ψall        = zeros(Nt+1,length(ψ0))
        end
        tspan       = zeros(Nt+1)
        ψall[1,:]   = ψ0
        t           = 0.
        
        while t < tend
            ψn = ode_rk4(ψ0,dt,t,eom)
            ψ0 = ψn
            t  += dt
            step_number += 1
            if mod(step_number,ΔNt) == 0
                println("t=",t, " \n")
                if dims == 1
                    ψall[save_number+1,:] = ψ0
                end
                tspan[save_number+1] = t
                save_number += 1
            end
        end
    end
end
