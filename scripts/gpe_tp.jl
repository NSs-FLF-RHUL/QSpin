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

    trap = 0.5 * (X.^2 + Y.^2)
    KE = 0.5 * (Kx.^2 * facx^2 + Ky.^2 * facy^2)
    ke_grid = facx^2



"""
    hamil(ψ::Array{ComplexF64}, time::Float64)

Setting the equation of motion for the target problem

    :param ψ: variable/vector/array associated with the problem. In this example, ψ is a two-dimensional complex field.
    :param time: the time of the problem

"""
function hamil(ψ::Array{ComplexF64}, time::Float64)
    ke = ifft(KE.*fft(ψ))
    pot = (trap.-μ) .* ψ + g .* (abs.(ψ).^2).* ψ - Ω * im .* (Y .* ifft(im .* Kx .* facx .*fft(ψ))-X.*ifft( im .* Ky .* facy .*fft(ψ)))
    dψdt = -1 * (ke+pot)
end

"""
    gpe2D(Dt::Float64,tend::Float64)

Setting the equation of motion for the target problem
    :param Dt: the time span on pulling out results    
    :param time: the total running time for the problem

"""
function gpe2D(dt::Float64,Dt::Float64,tend::Float64)    
    ψ0::Array{ComplexF64} = exp.(-((X).^2+Y.^2)/2) .* X.*Y + randn(ComplexF64, (length(x), length(y)))
    ψ, t = QSpin.OdeSolve.evolve_rk4(ψ0,dt,Dt,tend,hamil)
    return ψ, t 
end

ψt, t = gpe2D(5e-3,0.5,20.)
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
#gif(anim, "outputs/gpe_imt.gif", fps=50)
# To save as mp4 instead:
# gif(anim, "sine_wave.mp4", fps=30)