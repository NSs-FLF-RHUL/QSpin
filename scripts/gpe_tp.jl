using QSpin
using FFTW
using MAT, Random
using Plots, LaTeXStrings
FFTW.set_num_threads(4) 
# Parameter setup for the GPE simulation
g = 1
μ = 50
Ω = 0.75
Xmax = 22.
Ymax = 22.
Nx   = 256
Ny   = 256
irt = -1.
dt = 2.5e-4
Dt = .5
tend = 50.

# Grid setup
    x, y, X, Y, kx, ky,Kx, Ky, facx, facy = QSpin.Grids.CartGrid([Xmax, Ymax], [Nx, Ny])
    trap = 0.5 * (X.^2 + Y.^2)
    KE_mtx = 0.5 * (Kx.^2 + Ky.^2)

# Setting the equation of motion for the target problem
"""
    KE(ψ::Array{ComplexF64}, time::Float64)
    
    Creating the function of kinetic energy term in the Hamiltonian.

"""
KE = QSpin.Grids.fft_ke(KE_mtx)
   
"""
    Pot(trap::Array{Float64}, ψ::Array{ComplexF64})
    
    Creating the fucntion of the potential term in the Hamiltonian.

    Here we take the GPE in the rotating frame as an example, so the potential term includes the external trap potential, the nonlinear interaction term, and the rotation term.

    :param trap: the external trap potential, which is a two-dimensional array in this example.
    :param ψ: the field variable, which is a two-dimensional complex array in this example
"""
#function Pot(ψ::Array{ComplexF64},time)
#    pot = (trap.-μ) .* ψ + g .* (abs.(ψ).^2).* ψ - Ω * QSpin.Grids.fft_Lzψ(X,Y,Kx,Ky)(ψ,time) 
#    return pot
#end

Lz = QSpin.Grids.fft_Lzψ(X,Y,Kx,Ky) 
function Pot(ψ::Array{ComplexF64},time)
    pot = (trap.-μ) .* ψ + g .* (abs.(ψ).^2).* ψ
    ang_mom = Ω * Lz(ψ,time) 
    return pot - ang_mom
end

Hamil = QSpin.Hamiltonian.hamiltonian(KE,Pot,irt)

# Initial condition
ψ0::Array{ComplexF64} = sqrt(50.0).*exp.(-((X).^2+Y.^2)/2) .* X.*Y + 0.01*randn(ComplexF64, (length(x), length(y)))

# Solving the GPE
println("Simulation begins")
@time ψt, t = QSpin.OdeSolve.evolve_rk4(ψ0,dt,Dt,tend,Hamil);
println("Simulation finishes.")

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