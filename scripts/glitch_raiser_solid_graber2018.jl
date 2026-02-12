# This is an example script solve the simple glitch model under solid-body rotation approach
using Plots
using QSpin

# Parameter Setup according to V. Graber et al., ApJ 865, 23 (2018)
Msun = 1.989e30 # solar mass in kg
R = 10.e3# NS radius in m
BSF = 1e-2 #[1e-3,1e-2]
Bcore = 1e-2 # [3e-3,1e-2]
Itot = 0.35 * Msun * R^2
ISF = 0.95 * Itot
Icore = 0.05 * Itot
Icrust = 0.01599339442028418 * Itot
Next = 1.2e2

Ω0 = [70.34 - 6.3e-3; 70.34; 70.34]

"""
Setting the equation of motion for the target problem

    :param ψ: variable/vector/array associated with the problem.
    "param time: the time of the problem
"""
function eom(ψ::Array{Float64}, time::Float64)
    Mtx = [-2 1; 1 -2]
    return Mtx * ψ
end

# Running and Plotting
@time ut, t = QSpin.OdeSolve.evolve_rk4([0.1; 0.2], 1e-3, 1e-1, 1., eom)
output_plot = plot(t, ut[:, 1])
plot!(output_plot, t, ut[:, 2], xlabel="time (A.U.)", ylabel="Rotating Frequency (A.U.)", title="Solving a set of coupled ODEs")

savefig(output_plot, "./outputs/output-plot.png")

