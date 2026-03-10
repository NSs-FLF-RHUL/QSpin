using QSpin
using FFTW
using MAT, Random
using Plots, LaTeXStrings

# Parameter setup for the GPE simulation

g = 1
μ = 50
Ω = 0.5
Xmax = 20.
Ymax = 20.
Nx   = 256
Ny   = 256


    x, y, X, Y, kx, ky,Kx, Ky, facx, facy = QSpin.Grids.CartGrid([Xmax, Ymax], [Nx, Ny])
    Kx = Kx * facx;
    Ky = Ky * facy;
    trap = 0.5 * (X.^2 + Y.^2)
    KE_mtx = 0.5 * (Kx.^2 + Ky.^2)

"""
    potential(trap::Array{Float64}, ψ::Array{ComplexF64})
    
    the potential term in the GPE, which includes the external trap and the nonlinear interaction. The chemical potential μ is also included in the potential term for convenience.

    :param trap: the external trap potential, which is a two-dimensional array in this example.
    :param ψ: the field variable, which is a two-dimensional complex array in this example
"""
function potential(trap::Array{Float64}, ψ::Array{ComplexF64})
    pot = (trap.-μ) .* ψ + g .* (abs.(ψ).^2).* ψ
    return pot
end

"""
    hamil(ψ::Array{ComplexF64}, time::Float64,irt::ComplexF64)

Setting the equation of motion for the target problem

    :param ψ: variable/vector/array associated with the problem. In this example, ψ is a two-dimensional complex field.
    :param time: the time of the problem
    :param irt: propagation factor. -1 for imaginary time evolution, -im for real time evolution, or can be a general complex value for a dissipative evolution, namely, dGPE.

"""
irt = -1.

function Hamil(irt::Union{ComplexF64,Float64})
    function hamil(ψ::Array{ComplexF64}, time::Float64)
   
        keψ  = QSpin.Grids.fft_ke(KE_mtx)(ψ)
        potψ = potential(trap, ψ)
        Lzψ  = QSpin.Grids.fft_Lzψ(X,Y,Kx,Ky)(ψ) 

#        if sum(isnan.(keψ))>0
#            println("NaN detected in kinetic energy term at time ", time, ". Check the stability of the simulation and consider reducing the time step.")            
#        end
#            if sum(isnan.(Lzψ))>0
#            println("NaN detected in angular momentum term at time ", time, ". Check the stability of the simulation and consider reducing the time step.")            
#        end
        return irt .* (keψ + potψ - Ω * Lzψ)
    end
    return hamil
end

"""
    gpe2D(ψ0::Union{AbstractArray,Array{Float64},Array{ComplexF64}},dt::Float64,Dt::Float64,tend::Float64)

Setting the equation of motion for the target problem
    :param Dt: the time span on pulling out results    
    :param time: the total running time for the problem

"""
function gpe2D(ψ0::Union{AbstractArray,Array{Float64},Array{ComplexF64}},dt::Float64,Dt::Float64,tend::Float64)    
        ψ, t = QSpin.OdeSolve.evolve_rk4(ψ0,dt,Dt,tend,Hamil(irt))
    return ψ, t 
end

ψ0::Array{ComplexF64} = exp.(-((X).^2+Y.^2)/2) .* X.*Y + randn(ComplexF64, (length(x), length(y)))

ψt, t = gpe2D(ψ0, 5e-4, 0.5, 25.)
println("Fin.")

# Create animation
anim = @animate for i in 1:length(t)
    heatmap!(x,y,abs.(ψt[:,:,i]).^2, 
         title=string(L"\mathrm{Rotating BEC}, \Omega=",Ω, L" \omega, i\omega\tau=",round(t[i]*10)/10),
         xlims=(-Xmax, Xmax), 
         ylims=(-Ymax, Ymax),
         aspect_ratio = 1.,
         clim=(0,50))
end

# Save as GIF
gif(anim, "outputs/gpe_imt.gif", fps=50)
# To save as mp4 instead:
# gif(anim, "sine_wave.mp4", fps=30)