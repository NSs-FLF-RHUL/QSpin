using QSpin
using FFTW
using MAT
using Plots, LaTeXStrings

# Parameter setup for the GPE simulation
dt_factor = 0.01
Xmax = 5
Ymax = 4
Nx   = 128 
Ny   = 96 

dx   = 2 * Xmax / Nx
dy   = 2 * Ymax / Ny

x    = range(-Xmax,stop=Xmax-dx,length=Int(Nx))
y    = range(-Ymax,stop=Ymax-dy,length=Int(Ny))

facx = pi / Xmax
facy = pi / Ymax

kx   = [range(0,stop=Nx/2-1,length=Int(Nx/2)) ;range(-Nx/2,stop=-1,length=Int(Nx/2))]
ky   = [range(0,stop=Ny/2-1,length=Int(Ny/2)) ;range(-Ny/2,stop=-1,length=Int(Ny/2))]

X   = repeat(x',Ny,1)
Y   = repeat(y,1,Nx)

Kx = repeat(kx',Ny,1)
Ky = repeat(ky,1,Nx);

trap = 0.5 * (X.^2 + Y.^2)
KE = 0.5 * (Kx.^2 * facx^2 + Ky.^2 * facy^2)
ke_grid = facx^2

"""
    hamil(ψ::Array{ComplexF64}, time::Float64)

Setting the equation of motion for the target problem

    :param ψ: variable/vector/array associated with the problem. In this example, ψ is a two-dimensional complex field.
    "param time: the time of the problem

"""
function hamil(ψ::Array{ComplexF64}, time::Float64)
    ke = ifft(KE.*fft(ψ))
    pot = trap .* ψ
    dψdt = -1im * (ke+pot)
end

"""
    gpe2D(tend::Float64)

Setting the equation of motion for the target problem

    "param time: the total running time for the problem

"""
function gpe2D(tend::Float64)

    save_number = 1
    irt = -1im

    t = 0.
    dt = 1e-3
    Dt = 0.1
    ψ0::Array{ComplexF64} = exp.(-((X.-2).^2+Y.^2)/2)
    #tend = 1
    print("tend= ",tend ,"\n")
    save_number = 1

    ψ, t = QSpin.OdeSolve.evolve_rk4(ψ0,dt,Dt,tend,hamil)
    return ψ, t 
end

ψt, t = gpe2D(2*π)
println("Fin.")

# Create animation
anim = @animate for i in 1:length(t)
    heatmap!(x,y,abs.(ψt[:,:,i]).^2, 
         title="coherent state",
         xlims=(-Xmax, Xmax), 
         ylims=(-Ymax, Ymax))
end

# Save as GIF
gif(anim, "outputs/coherent.gif", fps=30)
# To save as mp4 instead:
# gif(anim, "sine_wave.mp4", fps=30)