using QSpin
using FFTW
using MAT, Random
using Plots, LaTeXStrings

import OrdinaryDiffEq as DE

Random.seed!(42)
FFTW.set_num_threads(6)

GridParams = (Xmax = 22.0, Ymax = 22.0, Nx = 256, Ny = 128)
ParaIn = (g = 1.0, μ = 50.0, Ω = 0.75, irt = -1.0)
SolverParams = (dt = 5e-4, Dt = 0.5, t_start = 0.0, t_end = 1.0)

# Grid setup
x, y, X, Y, kx, ky, Kx, Ky, facx, facy =
    QSpin.Grids.CartGrid([GridParams.Xmax, GridParams.Ymax], [GridParams.Nx, GridParams.Ny])
trap = 0.5 * (X .^ 2 + Y .^ 2)
KE_mtx = 0.5 * (Kx .^ 2 + Ky .^ 2)

# Setting the equation of motion for the target problem
PFFT = plan_fft(zeros(ComplexF64, GridParams.Ny, GridParams.Nx), [1, 2])
PiFFT = plan_ifft(zeros(ComplexF64, GridParams.Ny, GridParams.Nx), [1, 2])

KE = QSpin.Grids.Pfft_ke(KE_mtx, PFFT, PiFFT);
Lz = QSpin.Grids.fft_Lzψ(X, Y, Kx, Ky)

function Pot(ψ::Array{ComplexF64}, parameters::NamedTuple, time::Float64)
    pot = (trap .- parameters.μ) .* ψ + parameters.g .* (abs.(ψ) .^ 2) .* ψ
    ang_mom = parameters.Ω * Lz(ψ, parameters, time)
    return pot - ang_mom
end

Hamil! = QSpin.Hamiltonian.hamiltonian!(KE, Pot, ParaIn.irt)

ψ0::Array{ComplexF64} =
    sqrt(50.0) .* exp.(-((X) .^ 2+Y .^ 2)/2) .* X .* Y +
    0.01*randn(ComplexF64, (length(y), length(x)))

@time ψt = QSpin.OdeSolve.evolve(
    Hamil!,
    ψ0,
    SolverParams.t_start,
    SolverParams.t_end,
    ParaIn;
    alg = DE.Tsit5(),
    dt = SolverParams.dt,
    saveat = SolverParams.Dt,
)

# Create animation
anim = @animate for i = 1:length(ψt.t)
    heatmap!(
        x,
        y,
        abs.(ψt[:, :, i]) .^ 2,
        title = string(
            L"\mathrm{Rotating BEC}, \Omega=",
            ParaIn.Ω,
            L" \omega, i\omega\tau=",
            round(ψt.t[i] * 10) / 10,
        ),
        xlims = (-GridParams.Xmax, GridParams.Xmax),
        ylims = (-GridParams.Ymax, GridParams.Ymax),
        aspect_ratio = 1.0,
        clim = (0, 50),
    )
end

gif(anim, string("outputs/gpe_imt-Om=", ParaIn.Ω, "-DE.gif"), fps = 50)
