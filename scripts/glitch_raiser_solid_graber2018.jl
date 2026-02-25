# This is an example script solve the simple glitch model under solid-body rotation approach
using Plots, LaTeXStrings
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

BSF = 2.5e-4
Bcore = 5e-4
EMbrake = 1e-10
ISF = 0.95
Itot = 1;
Icore = 0.05

Ω0 = [70.34; 70.34 + 6.3e-3; 70.34]

"""
    eom(ψ::Array{Float64}, time::Float64)

Setting the equation of motion for the target problem

    :param ψ: variable/vector/array associated with the problem. In this example, ψ=[Ω_crust; Ω_SF Ω_core]
    "param time: the time of the problem

    In this example script, we solve the simpliest three-component glitch model in V. Graber et al., ApJ 865, 23 (2018):
        Ω̇_crust = -N_ext/I_crust - I_core * Ω̇_core / I_crust - I_SF * Ω̇_SF / I_crust
        Ω̇_SF = 2 * B * Ω_SF * (Ω_crust - Ω_SF)
        Ω̇_core = 2 * B_core * Ω_core * (Ω_crust - Ω_core)
    with the initial condition set in Ω0 = [Ω_core(t=0) ; Ω_SF(t=0) ; Ω_core(t=0) ].
"""
function eom(ψ::Array{Float64}, time::Float64)
    Mtx = [0 0 0;1 -1 0; 1 0 -1]
    vec = Mtx * ψ
    Mtx_2 = [0 0 0; 0 2*BSF*ψ[2] 0; 0 0 2*Bcore*ψ[3]]
    vec   = Mtx_2 * vec
    vec[1] = -EMbrake - ISF/Itot * vec[2] - Icore/Itot * vec[3]
    return vec
end

# Running and Plotting
@time ut, t = QSpin.OdeSolve.evolve_rk4(Ω0, 5e-3, 1e-1, 120.0, eom)
output_plot = plot(t, ut[1, :],label=L"$\Omega_\mathrm{crust}$")
#plot!(output_plot,t,Ω0[1].-EMbrake*t)
plot!(output_plot,t,ut[2,:],label=L"$\Omega_\mathrm{SF}$")
plot!(
    output_plot,
    t,
    ut[3, :],
    label=L"$\Omega_\mathrm{core}$",
    xlabel = "time (s)",
    ylabel = "Rotating Frequency (Hz)",
    title = "Simple Glitch Raiser Sim.",
)

savefig(output_plot, "./outputs/output-plot.png")
