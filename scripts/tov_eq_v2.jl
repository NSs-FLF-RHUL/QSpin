import OrdinaryDiffEq as DE
using OrdinaryDiffEqLowOrderRK: OrdinaryDiffEqLowOrderRK
using QSpin.TOV: EoS_two_component_polytrope, tov_eq!, EoS_NegeleVautherin1973, EoS_GCA2018
using QSpin.OdeSolve: evolve
using QSpin.Parameters: ParameterType
using SciMLBase: ODEProblem, DiscreteCallback, terminate!
using QSpin.PhysicalConstants:
    gravitational_constant, speed_of_light_vacuum, neutron_mass, hbar
using Plots

ħ = hbar * 1e3 * 1e4; # convert to g * cm^2 / s
mn = neutron_mass * 1e3; # convert to g

EoS_Param = (
    K_crust = (3*π^2)^(2/3) * ħ^2 / (5*mn^(8/3)), # Polytropic constant for the crust
    γ_crust = 5.0/3.0, # Polytropic index for the crust
    γ_core = 3.0, # Polytropic index for the core
    ρ_b = 3e14, # Transition density between crust and core in kg/m^3
)

# Simulation Input Parameters
TOV_Sol_Input = (
    ρ0 = 1e15, # Initial central density in (1g/cm^3)
    dr = 0.01*1e5, # Radial step in cm
    Dr = 0.01*1e5, # Radial interval for recording values in cm
    r_end = 20*1e5, # Maximum radius to solve up to in cm
);

#EoS, EoS_inv = EoS_two_component_polytrope(EoS_Param);
#EoS, EoS_inv = EoS_NegeleVautherin1973();
EoS, EoS_inv = EoS_GCA2018();
u0 = [EoS(TOV_Sol_Input.ρ0); 0.0]; # Initial conditions: central pressure and enclosed mass
tov! = tov_eq!(EoS_inv);

#condition(u, t, integrator) = u[1] < 0
#affect!(integrator) = terminate!(integrator)
#cb = DiscreteCallback(condition, affect!)

#problem = ODEProblem(tov!, u0, (0.0, TOV_Sol_Input.r_end); callback = cb)
#sol =
#    DE.solve(problem, OrdinaryDiffEqLowOrderRK.DP5(); reltol = 1e-12, abstol = [1.0, 1e15])
#sol = DE.solve(problem, DE.Tsit5(); reltol = 1e-12)

#ur = Array(sol)
#r = sol.t
#Pr = ur[1, :]
#mr = ur[2, :]
#ρr = EoS_inv.(Pr)

TOV_sol = QSpin.TOV.TOV_Solve(
    [EoS(TOV_Sol_Input.ρ0); 0.0],
    TOV_Sol_Input.dr,
    TOV_Sol_Input.Dr,
    TOV_Sol_Input.r_end,
    EoS_inv;
    alg = DE.Tsit5(),
    reltol = 1e-8,
    abstol = [1.0, 1e15],
);



plot(TOV_sol.r*1e-5, TOV_sol.ρr)
