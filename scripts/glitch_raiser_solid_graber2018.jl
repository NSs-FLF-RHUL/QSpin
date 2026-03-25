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
function eom_desolve(dψ::AbstractArray, ψ::AbstractArray, parameters, time::Float64)
    dψ .= M * ψ
end

de_problem = DE.ODEProblem(eom_desolve, ψ0, t_span)

@time u_desolve = DE.solve(de_problem, DE.Tsit5(), dt = dt, saveat = Dt)
t_desolve = u_desolve.t

desolve_plot = plot(t_desolve, u_desolve[1, :])
plot!(
    desolve_plot,
    t_desolve,
    u_desolve[2, :],
    xlabel = "time (A.U.)",
    ylabel = "Rotating Frequency (A.U.)",
    title = "Solving a set of coupled ODEs [DESolve]",
)

function eom(ψ::AbstractArray, time::Float64)
    return M * ψ
end

@time u_qspin, t_qspin = QSpin.OdeSolve.evolve_rk4([0.1; 0.2], 1e-3, 1e-1, 1.0, eom)

qspin_plot = plot(t_qspin, u_qspin[1, :])
plot!(
    qspin_plot,
    t_qspin,
    u_qspin[2, :],
    xlabel = "time (A.U.)",
    ylabel = "Rotating Frequency (A.U.)",
    title = "Solving a set of coupled ODEs [QSpin]",
)

savefig(desolve_plot, "./outputs/desolve.png")
savefig(qspin_plot, "./outputs/qspin.png")

println("Timestep difference:", maximum(abs2, (t_qspin - t_desolve)))
println("Solution difference:", maximum(abs2, (u_qspin - u_desolve)))
