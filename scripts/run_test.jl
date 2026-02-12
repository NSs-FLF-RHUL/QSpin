# This is an example script solve the simple glitch model under solid-body rotation approach
using Plots
using QSpin
"""
Setting the equation of motion for the targed problem

    :param ψ: vairable/vector/array assocaited with the problem.
    "param time: the time of the problem
"""
function eom(ψ::Array{Float64},time::Float64)
    A = [-2 1;1 -2]
    dψdt = A * ψ
end

# Running and Plotting
@time ut, t = QSpin.OdeSolve.evolve_rk4([.1;.2],1e-3,1e-1,1.,eom)
plot(t,ut[1,:])
plot!(t,ut[2,:],xlabel = "time (A.U.)",ylabel="Rotating Frequency (A.U.)",title="Solving a set of coupled ODE")
