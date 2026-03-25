# This is an example script solve the simple glitch model under solid-body rotation approach
using Plots
using QSpin

import OrdinaryDiffEq as DE

M = [-2 1; 1 -2]
ψ0 = [0.1; 0.2]
dt = 1e-3
Dt = 1e-1
t_start = 0.0
t_end = 1.0
t_span = (t_start, t_end)

"""
Setting the equation of motion for the target problem

    :param ψ: variable/vector/array associated with the problem.
    "param time: the time of the problem
"""
function eom(dψ::AbstractArray, ψ::AbstractArray, parameters, time::Float64)
    dψ .= M * ψ
end

problem = DE.ODEProblem(eom, ψ0, t_span)
u = DE.solve(problem, DE.Tsit5(), dt = dt, saveat = Dt)

output_plot = plot(u_desolve.t, u_desolve[1, :])
plot!(
    output_plot,
    u.t,
    u[2, :],
    xlabel = "time (A.U.)",
    ylabel = "Rotating Frequency (A.U.)",
    title = "Solving a set of coupled ODEs [DESolve]",
)

savefig(output_plot, "./outputs/glitch-raiser-output.png")
