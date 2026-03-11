using QSpin
using FFTW
using MAT, Random
using Plots, LaTeXStrings

# Parameter setup for the GPE simulation
irt = -1.
g = 1
μ = 50
Ω = 0.75
Xmax = 22.
Ymax = 22.
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
function potential(ψ::Array{ComplexF64},time)
    pot = (trap.-μ) .* ψ + g .* (abs.(ψ).^2).* ψ - Ω * QSpin.Grids.fft_Lzψ(X,Y,Kx,Ky)(ψ,time) 
    return pot
end


function KE(ψ::Union{AbstractArray,Array{Float64},Array{ComplexF64}}, time::Float64)
    return QSpin.Grids.fft_ke(KE_mtx)(ψ::Union{AbstractArray,Array{Float64},Array{ComplexF64}}, time::Float64)
end
Hamil = QSpin.Hamiltonian.hamiltonian(KE,potential,irt)

dt = 5e-4
Dt = .5
tend = 25.

ψ0::Array{ComplexF64} = sqrt(50.0).*exp.(-((X).^2+Y.^2)/2) .* X.*Y + 0.01*randn(ComplexF64, (length(x), length(y)))


#function gpe2D(ψ0::Union{AbstractArray,Array{Float64},Array{ComplexF64}},dt::Float64,Dt::Float64,tend::Float64)    
        ψt, t = QSpin.OdeSolve.evolve_rk4(ψ0,dt,Dt,tend,Hamil)
#    return ψ, t 
#end



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
gif(anim, string("outputs/gpe_imt-Om=", Ω, ".gif"), fps=50)
# To save as mp4 instead:
# gif(anim, "sine_wave.mp4", fps=30)