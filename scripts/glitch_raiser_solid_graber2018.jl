# This is an example script solve the simple glitch model under solid-body rotation approach
using Plots
using QSpin
using QSpin.Parameters: ParameterType

import OrdinaryDiffEq as DE

M = [-2 1; 1 -2]
ψ0 = [0.1; 0.2]
dt = 1e-3
Dt = 1e-1
t_start = 0.0
t_end = 1.0
t_span = (t_start, t_end)

function eom!(dψ::AbstractArray, ψ::AbstractArray, parameters::ParameterType, time::Float64)
    dψ .= M * ψ
end

u = QSpin.OdeSolve.evolve(eom!, ψ0, t_start, t_end; alg = DE.Tsit5(), dt = dt, saveat = Dt)

output_plot = plot(u.t, u[1, :])
plot!(
    output_plot,
    u.t,
    u[2, :],
    xlabel = "time (A.U.)",
    ylabel = "Rotating Frequency (A.U.)",
    title = "Solving a set of coupled ODEs",
)

savefig(output_plot, "./outputs/glitch-raiser-output.png")
