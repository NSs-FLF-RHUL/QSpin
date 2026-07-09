# This is an example script solve the simple glitch model under solid-body rotation approach
using QSpin
using QSpin.Parameters: ParameterType
import QSpin.OdeSolve: evolve
import QSpin.GlitchModels: ThreeCompSolid!
import CommonSolve as DE
using Plots

Sim_Input = (
    # Glitch Model Input
    Ω_crust = 70.34, # Initial angular velocity of the crust in rad/s
    Ω_sf = 70.34 + 6.3e-3, # Initial angular velocity of the superfluid in rad/s
    Ω_core = 70.34, # Initial angular velocity of the core in rad/s
    I_crust = 4.5e30, # Moment of inertia of the crust in kg m^2
    I_core = 0.8 * 4.5e30, # Moment of inertia of the core in kg m
    N_ext = 0.0,
    B_core = 5e-4, # Mutual Friction Parameter
    B_sf = 5e-4, # Mutual Friction Parameter
    # Glitch Model Solver Setup
    Dt = 0.1, # Time interval for recording values in the glitch model in seconds
    t_start = 0.0, # Start time for the glitch model simulation in seconds
    t_end = 120.0, # End time for the glitch model simulation in seconds
)

Ω_ini = [Sim_Input.Ω_crust; Sim_Input.Ω_sf; Sim_Input.Ω_core]
sol = evolve(
    ThreeCompSolid!,
    Ω_ini,
    0.0,
    Sim_Input.t_end,
    Sim_Input;
    alg = DE.Tsit5(),
    saveat = Sim_Input.Dt,
)

output_plot = plot(sol.t, sol[1, :])
plot!(
    output_plot,
    u.t,
    u[2, :],
    xlabel = "time (A.U.)",
    ylabel = "Rotating Frequency (A.U.)",
    title = "Solving a set of coupled ODEs",
)

savefig(output_plot, "./outputs/glitch-raiser-output.png")
