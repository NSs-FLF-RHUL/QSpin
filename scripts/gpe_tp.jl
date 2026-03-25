using QSpin
using FFTW
using MAT, Random
using Plots, LaTeXStrings

import OrdinaryDiffEq as DE

Random.seed!(42)
FFTW.set_num_threads(6)

ParaIn = (
    g = 1.0,
    μ = 50.0,
    Ω = 0.75,
    Xmax = 22.0,
    Ymax = 22.0,
    Nx = 256,
    Ny = 256,
    irt = -1.0,
    dt = 5e-4,
    Dt = 0.5,
    tstart = 0.0,
    tend = 1.0,
)
t_span = (ParaIn.tstart, ParaIn.tend)

# Grid setup
x, y, X, Y, kx, ky, Kx, Ky, facx, facy =
    QSpin.Grids.CartGrid([ParaIn.Xmax, ParaIn.Ymax], [ParaIn.Nx, ParaIn.Ny])
trap = 0.5 * (X .^ 2 + Y .^ 2)
KE_mtx = 0.5 * (Kx .^ 2 + Ky .^ 2)

KE = QSpin.Grids.fft_ke(KE_mtx)
Lz = QSpin.Grids.fft_Lzψ(X, Y, Kx, Ky)

function Pot(ψ::Array{ComplexF64}, parameters::NamedTuple, time)
    pot = (trap .- parameters.μ) .* ψ + parameters.g .* (abs.(ψ) .^ 2) .* ψ
    ang_mom = parameters.Ω * Lz(ψ, parameters, time)
    return pot - ang_mom
end

Hamil! = QSpin.Hamiltonian.hamiltonian!(KE, Pot, ParaIn.irt)

ψ0::Array{ComplexF64} =
    sqrt(50.0) .* exp.(-((X) .^ 2 + Y .^ 2) / 2) .* X .* Y +
    0.01 * randn(ComplexF64, (length(x), length(y)))

problem = DE.ODEProblem(Hamil!, ψ0, t_span, ParaIn)

# Solving the GPE
println("Simulation begins")
@time ψt = DE.solve(problem, DE.Tsit5(), dt = ParaIn.dt, saveat = ParaIn.Dt);
println("Simulation finishes.")
t = ψt.t

# Create animation
anim = @animate for i = 1:length(t)
    heatmap!(
        x,
        y,
        abs.(ψt[:, :, i]) .^ 2,
        title = string(
            L"\mathrm{Rotating BEC}, \Omega=",
            ParaIn.Ω,
            L" \omega, i\omega\tau=",
            round(t[i] * 10) / 10,
        ),
        xlims = (-ParaIn.Xmax, ParaIn.Xmax),
        ylims = (-ParaIn.Ymax, ParaIn.Ymax),
        aspect_ratio = 1.0,
        clim = (0, 50),
    )
end

gif(anim, string("outputs/gpe_imt-Om=", ParaIn.Ω, "-DE.gif"), fps = 50)
