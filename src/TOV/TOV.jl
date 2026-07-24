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

using ..PhysicalConstants: gravitational_constant, speed_of_light_vacuum
using DocStringExtensions: TYPEDSIGNATURES
using OrdinaryDiffEqLowOrderRK: OrdinaryDiffEqLowOrderRK
using SciMLBase: ContinuousCallback, terminate!, ODEProblem
using CommonSolve: CommonSolve
using DataInterpolations: QuadraticSpline

# Include the equations of state module as a submodule of the TOV module
include("EoS/EquationOfState.jl")

"""

$(TYPEDSIGNATURES)

Returns a named tuple containing the reference units for the TOV equation, given a reference density,rho_ref, in the CGS or SI units.

"""
function TOV_ref_units(; input_units::String = "CGS", rho_ref::Float64 = 2.8e14)
    if input_units in ("CGS", "CGS_dim")
        #println(" CGS unit")
        G0 = gravitational_constant * 1e3
        c0 = speed_of_light_vacuum * 1e2
    elseif input_units in ("SI", "SI_dim")
        #println(" SI unit")
        rho_ref = rho_ref * 1e3
        G0 = gravitational_constant
        c0 = speed_of_light_vacuum
    else
        error(
            "Unsupported units '$(input_units)'.  Only supported units are 'CGS' ('CGS_dim') and 'SI' ('SI_dim')",
        )
    end
    if input_units in ("CGS_dim", "SI_dim")
        length_ref = 1.0
        pressure_ref = 1.0
        mass_ref = 1.0
        rho_ref = 1.0
    else
        length_ref = c0 / sqrt(G0*rho_ref)
        pressure_ref = rho_ref * c0^2
        mass_ref = rho_ref * length_ref^3
        G0 = 1.0
        c0 = 1.0
    end
    tovUnits = (; length_ref, pressure_ref, mass_ref, rho_ref, G0, c0, input_units)
    return tovUnits
end

"""

$(TYPEDSIGNATURES)

Returns a function that evaluates the RHS of the TOV equations, given an equation of state function.

# Arguments
- `EoS_rho_from_P::Function`: Single-argument inverse equation-of-state function that maps pressure (``P``) to density (``\\rho``).

# Returns
- `tov_inner!::Function`: Callable as `tov_inner!(du, u, params, r)` that evaluates the RHS of the TOV equations, writing the result to `du`.
"""
function tov_eq!(
    EoS_rho_from_P::Union{Function,QuadraticSpline};
    units::Union{String} = "CGS",
    rho_ref::Float64 = 2.8e14, # nuclear saturation density in g/cm^3
)
    tovUnits = TOV_ref_units(input_units = units, rho_ref = rho_ref)
    P_ref = tovUnits.pressure_ref;
    rho_ref = tovUnits.rho_ref;
    G0 = tovUnits.G0;
    c0 = tovUnits.c0;

    function tov_inner!(du, u, params, r)
        P = u[1]
        m = u[2]
        rho = EoS_rho_from_P(P*P_ref) / rho_ref
        du[1] = if r == 0.0
            0.0
        else
            -G0 / r^2 * (rho + P/c0^2) * (m + 4*pi*r^3*P/c0^2) / (1 - 2*G0*m/(r*c0^2))
        end
        du[2] = 4*pi*r^2*rho
    end
    return tov_inner!
end

"""

$(TYPEDSIGNATURES)

Solve the TOV equation for given equation(s) of state and initial condition(s), wrapping `QSpin.OdeSolve.evolve`.

By default, the solver employs the DP5 method, which is a 5th order explicit Runge-Kutta method with an embedded 4th order method for error estimation, suitable for non-stiff problems. For stiff problems, consider using `KenCarp4` or `Rodas5` methods from the `OrdinaryDiffEqSDIRK` package.

# Arguments
- `u0::AbstractArray`: Initial conditions for the TOV equation, given as a vector of the form `[P(0), m(0)]`.
- `dr::Float64`: Radial step size for the numerical solver.
- `Dr::Float64`: Radial interval for recording values of the solution.
- `r_max`: Maximum radius in the selected input units. Defaults to 20 km.
- `EoS_inv::Function`: Single-argument inverse equation-of-state function that maps pressure (``P``) to density (``\\rho``).
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
    EoS_inv::Union{Function,QuadraticSpline},
    u0::AbstractArray,
    dr::Float64,
    Dr::Float64,
    r_beg::Float64;
    alg = OrdinaryDiffEqLowOrderRK.DP5(),
    dt = dr,
    saveat = Dr,
    input_units::Union{String} = "CGS",
    rho_ref::Float64 = 2.8e14, # nuclear saturation density in g/cm^3
    r_max::Union{Nothing,Float64} = nothing,
    reltol = 1e-13,
    abstol = 1e-13,
    solver_options...,
)

    is_si = input_units in ("SI", "SI_dim")
    r_max = isnothing(r_max) ? (is_si ? 25e3 : 25e5) : r_max

    if u0[1] <= 0
        r = [r_beg]
        Pr = [u0[1]]
        mr = [u0[2]]
        ρr = EoS_inv.(Pr)
        return (; r, Pr, mr, ρr, R = r_beg, M = u0[2])
    end

    # Building the dimensionless TOV equation function using the provided inverse EoS function
    tov! = tov_eq!(EoS_inv; units = input_units, rho_ref = rho_ref);
    # Scale the initial conditions and radial parameters to dimensionless units
    tovUnits = TOV_ref_units(; input_units = input_units, rho_ref = rho_ref);

    u0 = u0 ./ [tovUnits.pressure_ref, tovUnits.mass_ref]
    dt = dt / tovUnits.length_ref
    saveat = saveat / tovUnits.length_ref
    r_beg = r_beg / tovUnits.length_ref
    r_max = r_max / tovUnits.length_ref

    # Locate the stellar surface by root-finding on P(r) = 0.
    condition(u, t, integrator) = u[1]
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!; save_positions = (false, true))
    # Define the ODE problem and solve it with the requested algorithm and options.
    problem = ODEProblem(tov!, u0, (r_beg, r_max); callback = cb)
    sol = CommonSolve.solve(problem, alg; dt, saveat, reltol, abstol, solver_options...)

    ur = Array(sol)
    r = sol.t * tovUnits.length_ref
    Pr = ur[1, :] * tovUnits.pressure_ref
    mr = ur[2, :] * tovUnits.mass_ref
    ρr = EoS_inv.(Pr)
    # Output arguments in the original units
    TOV_sol = (; r, Pr, mr, ρr, R = r[end], M = mr[end])

    return TOV_sol
end

end
