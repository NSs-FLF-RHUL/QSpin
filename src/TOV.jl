"""
Sumbodule containing the helper fucntion for solvign the Tolman–Oppenheimer–Volkoff (TOV) equation, taht describes the equation of state for the neutron star interior under the non-relativistic approximation limit.

The TOV equation is a ordinary differential equation of pressure as a function of radius r and is written by

``
\\frac{\\mathrm{d}P}{\\mathrm{d}r} = - \\frac{G}{r^2}
\\left[ \\rho(r) + \\frac{P(r)}{c^2} \\right]
\\left[ m(r) + 4\\pi r^3 \\frac{P(r)}{c^2} \\right]
\\left[ 1 - \\frac{2Gm(r)}{c^2 r} \\right]^{-1},
``

where ``G`` and ``c`` represent the gravitational constant and speed of light in vacuum, respectively.

It is coupled to the mass m(r) within radius r, which is given by

``
\\frac{\\mathrm{d}m}{\\mathrm{d}r} = 4\\pi\\rho(r)r^2.
``

where ρ is the density as a function of r and is assumeed to be spherically symmetric.

To solve the TOV equation, we need to specify an equation of state (EoS) that relates the pressure P to the density ρ, that is we must provide ``P = P(\\rho)`` (and the inverse). Here, to solve the TOV equation numerically, we use the vector presentation for the variables, `u = (P, m)`` and solving a first-order system ODE in ``u``.
"""
module TOV

using ..OdeSolve: evolve
using ..Parameters: ParameterType
using ..PhysicalConstants:
    gravitational_constant, speed_of_light_vacuum, neutron_mass, electron_volt, hbar
using Roots: find_zero
using DocStringExtensions: TYPEDSIGNATURES
using OrdinaryDiffEqLowOrderRK: OrdinaryDiffEqLowOrderRK
using SciMLBase: DiscreteCallback, terminate!

function tov_eq!(EoS_rho_from_P::Function)
    function tov_inner!(du, u, paras, r)
        P = u[1]
        m = u[2]
        ρ = EoS_rho_from_P(P)
        du[1] = if r == 0.0
            0.0
        else
            -gravitational_constant / r^2 *
            (ρ + P/speed_of_light_vacuum^2) *
            (m + 4*pi*r^3*P/speed_of_light_vacuum^2) /
            (1 - 2*gravitational_constant*m/(r*speed_of_light_vacuum^2))
        end
        du[2] = 4*pi*r^2*ρ
    end
    return tov_inner!
end

"""

$(TYPEDSIGNATURES)

Solving the TOV equation for a given EoS and initial conditions using the `CommonSolve`  package, adapted in OdeSolve.jl of this package.

Now it uses DP5 method, which is a 5th order explicit Runge-Kutta method with an embedded 4th order method for error estimation, suitable for non-stiff problems. For stiff problems, consider using `KenCarp4` or `Rodas5` methods from the `OrdinaryDiffEqSDIRK` package, and this method is also used in [TOVsolver_Julia](https://github.com/jskMNMGCH/TOVsolver_Julia).

# Arguments
- `u0::AbstractArray`: Initial conditions for the TOV equation, given as a vector of the form `[P(0), m(0)]`, where `P(0)`
- 'dr::Float64`: Radial step size for the numerical solver.
- `Dr::Float64`: Radial interval for recording values of the solution.
- `r_max::Float64`: Maximum radius to solve up to.
- `EoS_inv::Function`: Inverse EoS function that takes pressure as input and returns density, i.e., `EoS_inv(P) = ρ`.
- `solver_options...`: Additional keyword arguments to pass to the ODE solver.

# Returns
- `TOV_sol::NamedTuple`: A named tuple containing the solution of the TOV equation, with the following fields:
    - `r`: Radial coordinates at which the solution is evaluated.
    - `Pr`: Pressure as a function of radius.
    - `mr`: Enclosed mass as a function of radius.
    - `ρr`: Density as a function of radius, obtained by applying the inverse EoS to the pressure solution.
    - `R`: The radius of the star, defined as the radius at which the pressure drops to zero.
    - `M`: The total mass of the star, defined as the enclosed mass at the radius `R`.
"""
function TOV_Solve(
    u0::Union{AbstractArray,Array{Float64}},
    dr::Float64,
    Dr::Float64,
    r_max::Float64,
    EoS_inv::Function;
    alg = OrdinaryDiffEqLowOrderRK.DP5(),
    dt = dr,
    saveat = Dr,
    reltol = 1e-12,
    solver_options...,
)
    tov! = tov_eq!(EoS_inv)
    condition(u, t, integrator) = u[1] < 0
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition, affect!)
    sol_tov = evolve(
        tov!,
        u0,
        0.0,
        r_max;
        alg,
        reltol,
        callback = cb,
        dt,
        saveat,
        solver_options...,
    )
    ur = Array(sol_tov)
    r = sol_tov.t
    Pr = ur[1, :]
    mr = ur[2, :]

    R_index = findfirst(x->x<0, Pr)
    TOV_sol = (;
        r,
        Pr,
        mr,
        ρr = EoS_inv.(Pr),
        R = r[R_index-1],
        M = mr[isnothing(R_index) ? end : (R_index - 1)],
    )
    return TOV_sol
end

"""
$(TYPEDSIGNATURES)

Create EoS functions for two-component polytrope neutron stars.

Such stars are assumed to be spheres composed of non-interacting, degenerate neutrons that are characterised by two EoS, in the core and crust regions.
The boundary between the regions is located at some density ``\\rho_b``;

```math
P = \\begin{cases}
K_\\mathrm{crust} \\rho^{\\gamma_{\\mathrm{crust}}} & \\rho < \\rho_b, \\\\
K_\\mathrm{core} \\rho^{\\gamma_{\\mathrm{core}}} & \\rho > \\rho_b,
\\end{cases}
\\quad
K_\\mathrm{core} = K_\\mathrm{crust} \\rho_b^{ \\gamma_{\\mathrm{crust}} - \\gamma_{\\mathrm{core}} },
```

with the equation for ``K_\\mathrm{core}`` being a result of imposing pressure continuity at the crust-core interface.

Assuming the crust is pure degenerate neutron matter, the non-relativistic EoS is given by

``
P_\\mathrm{crust} = \\frac{(3\\pi^2)^{2/3}}{5} \\frac{\\bar{h}^2}{m_n^{8/3}} \\rho^{5/3}.
``

namely,

``
K_\\mathrm{crust} = \\frac{(3\\pi^2)^{2/3}}{5} \\frac{\\bar{h}^2}{m_n^{8/3}}.
``

To keep the generality of the function, we keep ``K_\\mathrm{crust}`` as an input parameter here.

`parameters` is assumed to define the following quantities:

- ``K_\\mathrm{crust}``, under the name `K_crust`.
- ``\\gamma_\\mathrm{crust}``, under the name `gamma_crust`.
- ``\\gamma_\\mathrm{core}``, under the name `gamma_core`.
- ``\\rho_b``, under the name `rho_b`.

# Arguments
- `parameters::ParameterType`: Parameter values for the two-component polytrope model.
- `report_transition_pressure::Bool`: If `true`, the transition pressure at the interface will be printed.

# Returns
- `EoS_P_from_rho::Function`: ``P(\\rho)`` as a function called with `(rho,)`.
- `EoS_rho_from_P::Function`: ``\\rho(P)`` as a function called with `(P,)`.
"""
function EoS_two_component_polytrope(
    ParamIn::ParameterType,
    report_transition_pressure::Bool = false,
)
    K_core = ParamIn.K_crust * ParamIn.ρ_b ^ (ParamIn.γ_crust - ParamIn.γ_core)

    function EoS_P_from_rho(ρ)
        P = if ρ < 0
            0.0
        elseif ρ <= ParamIn.ρ_b
            ParamIn.K_crust * ρ .^ ParamIn.γ_crust
        else
            K_core * ρ .^ ParamIn.γ_core
        end
        return P
    end

    # Pressure at the crust-core transition for continuity check
    P_b = EoS_P_from_rho(ParamIn.ρ_b)

    if report_transition_pressure
        println("Pressure at crust-core transition (Pb): ", P_b)
    end

    function EoS_rho_from_P(P)
        ρ = if P < 0
            0.0
        elseif P <= P_b
            (P / ParamIn.K_crust) .^ (1/ParamIn.γ_crust)
        else
            (P / K_core) .^ (1/ParamIn.γ_core)
        end
        return ρ
    end

    return EoS_P_from_rho, EoS_rho_from_P
end


function EoS_NegeleVautherin1973(
    report_transition_pressure::Bool = false,
    ;
    ci = [
        -4.0,
        2.8822899e-1,
        5.9150523e-1,
        9.0185940e-2,
        -1.1025614e-1,
        2.9377479e-2,
        -3.2618465e-3,
        1.3543555e-4,
    ],
)
    function EoS_P_from_rho(ρ)
        P = if ρ < 0
            0.0
        else
            nb = ρ / (neutron_mass * 1e3); # in the unit of g/cm^3
            nb_scaled = nb * 1e-35
            x = log.(nb_scaled)
            energy_sum = zeros(size(nb))
            pressure_sum = zeros(size(nb))
            for j = 2:length(ci)
                i = j - 1
                energy_sum .+= ci[j] * x .^ (i-1)
                pressure_sum .+= (i-1) * ci[j] * x .^ (i-2)
            end
            P = (nb) .* pressure_sum .* exp.(energy_sum) * 1e6 * electron_volt * 1e7 # Convert to MeV fm^-3
        end
        return P
    end

    function EoS_rho_from_P(P)
        ρ = if P < 0
            0.0
        else
            find_zero(y -> EoS_P_from_rho(y) - P, 4e11)
        end
        return ρ
    end
    return EoS_P_from_rho, EoS_rho_from_P
end


"""
$(TYPEDSIGNATURES)

Create EoS functions according to V. Graber, A. Cumming & N. Anderson, ApJ 865, 23 (2018).
The EoS given by J. Negele & D. Vautherin, Nuclear Physics A 207, 298 (1973) is used for the inner crust regime and the outer crust is consdiered by

"""
function EoS_GCA2018(
    report_transition_pressure::Bool = false,
    ;
    ci = [
        -4.0,
        2.8822899e-1,
        5.9150523e-1,
        9.0185940e-2,
        -1.1025614e-1,
        2.9377479e-2,
        -3.2618465e-3,
        1.3543555e-4,
    ],
    ρ_drip = 4.3e11, # in g/m^3,
    Ye = 0.4,
)
    ħ = hbar * 1e3 * 1e4; # convert to g * cm^2 / s
    mn = neutron_mass * 1e3; # convert to g
    c = speed_of_light_vacuum
    function EoS_P_from_rho(ρ)
        P = if ρ < 0
            0.0
        elseif ρ < ρ_drip
            ħ * c * (3 * π^2 * Ye * ρ / (mn * 1e3))^(4/3) / 12 / π^2
        else
            nb = ρ / (neutron_mass * 1e3); # in the unit of g/cm^3
            nb_scaled = nb * 1e-35
            x = log.(nb_scaled)
            energy_sum = zeros(size(nb))
            pressure_sum = zeros(size(nb))
            for j = 2:length(ci)
                i = j - 1
                energy_sum .+= ci[j] * x .^ (i-1)
                pressure_sum .+= (i-1) * ci[j] * x .^ (i-2)
            end
            P = (nb) .* pressure_sum .* exp.(energy_sum) * 1e6 * electron_volt * 1e7 #
        end

        return P
    end

    function EoS_rho_from_P(P)
        ρ = if P < 0
            0.0
        elseif P < EoS_P_from_rho(ρ_drip)
            (12 * π^2 * P / ħ / c)^(3/4) * (mn * 1e3) / (3 * π^2 * Ye)
        else
            find_zero(y -> EoS_P_from_rho(y) - P, 5e11)
        end
        return ρ
    end
    return EoS_P_from_rho, EoS_rho_from_P
end

end
