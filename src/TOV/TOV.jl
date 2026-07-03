"""
Submodule containing the helper function for solving the Tolman-Oppenheimer-Volkoff (TOV) equation, that describes the equation of state for the neutron star interior under the non-relativistic approximation limit.

The TOV equation is a ordinary differential equation of pressure as a function of radius r and is written by

``
\\frac{ \\mathrm{d}P}{\\mathrm{d}r} = - \\frac{G}{r^2}
\\left[ \\rho(r) + \\frac{P(r)}{c^2} \\right]
\\left[ m(r) + 4\\pi r^3 \\frac{P(r)}{c^2} \\right]
\\left[ 1 - \\frac{2Gm(r)}{c^2 r} \\right]^{-1},
``

where ``G`` and ``c`` represent the gravitational constant and speed of light in vacuum, respectively.

It is coupled to the mass m(r) within radius r, which is given by

``
\\frac{\\mathrm{d}m}{\\mathrm{d}r} = 4\\pi\\rho(r)r^2.
``

where ``\\rho`` is the density as a function of r and is assumed to be spherically symmetric.

To solve the TOV equation, we need to specify an equation of state (EoS) that relates the pressure P to the density ``\\rho``, that is we must provide ``P = P(\\rho)`` (and the inverse). Here, to solve the TOV equation numerically, we use the vector presentation for the variables, ``u = (P, m)`` and solving a first-order system ODE in ``u``.
"""
module TOV

using ..OdeSolve: evolve
using ..Parameters: ParameterType
using ..PhysicalConstants: gravitational_constant, speed_of_light_vacuum
using DocStringExtensions: TYPEDSIGNATURES
using OrdinaryDiffEqLowOrderRK: OrdinaryDiffEqLowOrderRK
using SciMLBase: DiscreteCallback, terminate!, ODEProblem
import OrdinaryDiffEq as DE

# Include the equations of state module as a submodule of the TOV module
include("EoS/EquationOfState.jl")

"""

$(TYPEDSIGNATURES)

Returns a function that evaluates the RHS of the TOV equations, given an equation of state function.

# Arguments
- `EoS_rho_from_P::Function`: Single-argument (inverse) equation-of-state function that maps density (``\\rho``) to pressure (``P``).

# Returns
- `tov_inner!::Function`: Callable as `tov_inner!(du, u, params, r)` that evaluates the RHS of the TOV equations, writing the result to `du`.
"""
function tov_eq!(EoS_rho_from_P::Function)
    function tov_inner!(du, u, params, r)
        P = u[1]
        m = u[2]
        rho = EoS_rho_from_P(P)
        du[1] = if r == 0.0
            0.0
        else
            -gravitational_constant / r^2 *
            (rho + P/speed_of_light_vacuum^2) *
            (m + 4*pi*r^3*P/speed_of_light_vacuum^2) /
            (1 - 2*gravitational_constant*m/(r*speed_of_light_vacuum^2))
        end
        du[2] = 4*pi*r^2*rho
    end
    return tov_inner!
end



function TOV_ref_units(; units::Union{String} = "CGS", rho_ref::Float64 = 2.8e14)
    if units == "CGS"
        println(" CGS unit")
        G0 = gravitational_constant * 1e3
        c0 = speed_of_light_vacuum * 1e2
    elseif units == "SI"
        println(" SI unit")
        G0 = gravitational_constant
        c0 = speed_of_light_vacuum
        rho_ref = rho_ref * 1e3
    else
        error(" !!!! Non-Supported Units !!!!")
    end
    length_ref = c0 / sqrt(G0*rho_ref)
    pressure_ref = rho_ref * c0^2
    mass_ref = rho_ref * length_ref^3
    tovUnits = (; length_ref, pressure_ref, mass_ref, rho_ref)
    return tovUnits
end
"""

$(TYPEDSIGNATURES)

Returns a function that evaluates the RHS of the dimensionless TOV equations, given an equation of state function.

# Arguments
- `EoS_rho_from_P::Function`: Single-argument (inverse) equation-of-state function that maps dimensionless density (``\\hat{\\rho}``) to dimensionless pressure (``\\hat{P}``).

# Returns
- `tov_dimless_inner!::Function`: Callable as `tov_dimless!(du, u, params, r)` that evaluates the RHS of the dimensionless TOV equations, writing the result to `du`.
"""
function tov_eq_dimless!(EoS_rho_from_P::Function; Units::Union{String} = "CGS")
    tovUnits = TOV_ref_units(units = Units)
    P_ref = tovUnits.pressure_ref;
    rho_ref = tovUnits.rho_ref;
    function tov_dimless_inner!(du, u, params, r)
        P = u[1]
        m = u[2]
        rho = EoS_rho_from_P(P*P_ref) / rho_ref
        du[1] = if r == 0.0
            0.0
        else
            -1 / r^2 * (rho + P) * (m + 4*pi*r^3*P) / (1 - 2*m/r)
        end
        du[2] = 4*pi*r^2*rho
    end
    return tov_dimless_inner!
end

"""

$(TYPEDSIGNATURES)

Solve the TOV equation for given equation(s) of state and initial condition(s), wrapping `QSpin.OdeSolve.evolve`.

By default, the solver employs the DP5 method, which is a 5th order explicit Runge-Kutta method with an embedded 4th order method for error estimation, suitable for non-stiff problems. For stiff problems, consider using `KenCarp4` or `Rodas5` methods from the `OrdinaryDiffEqSDIRK` package.

# Arguments
- `u0::AbstractArray`: Initial conditions for the TOV equation, given as a vector of the form `[P(0), m(0)]`.
- `dr::Float64`: Radial step size for the numerical solver.
- `Dr::Float64`: Radial interval for recording values of the solution.
- `r_max::Float64`: Maximum radius to solve up to.
- `EoS_inv::Function`: Single-argument (inverse) equation-of-state function that maps density (``\\rho``) to pressure (``P``).
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
    u0::AbstractArray,
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
    # Building the TOV equation function using the provided inverse EoS function
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


function TOV_Solve_dimensionless(
    u0::AbstractArray,
    dr::Float64,
    Dr::Float64,
    r_beg::Float64,
    r_max::Float64,
    EoS_inv::Function;
    alg = OrdinaryDiffEqLowOrderRK.DP5(),
    dt = dr,
    saveat = Dr,
    reltol = 1e-12,
    solver_options...,
)
    # Building the dimensionless TOV equation function using the provided inverse EoS function
    tov! = tov_eq_dimless!(EoS_inv);
    tovUnits = TOV_ref_units()
    condition(u, t, integrator) = u[1] < 0
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition, affect!)

    problem = ODEProblem(tov!, u0, (r_beg, r_max); callback = cb)
    sol = DE.solve(problem, OrdinaryDiffEqLowOrderRK.DP5(); reltol = 1e-12)

    ur = Array(sol)
    r = sol.t
    Pr = ur[1, :]
    mr = ur[2, :]
    ρr = EoS_inv.(Pr*tovUnits.pressure_ref)
    R_index = findfirst(x->x<0, Pr)
    TOV_sol = (;
        r,
        Pr,
        mr,
        ρr = EoS_inv.(Pr),
        R = r[isnothing(R_index) ? end : (R_index - 1)],
        M = mr[isnothing(R_index) ? end : (R_index - 1)],
    )
    return TOV_sol
end




end
