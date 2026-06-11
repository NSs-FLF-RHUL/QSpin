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
using ..PhysicalConstants: gravitational_constant, speed_of_light_vacuum
using DocStringExtensions: TYPEDSIGNATURES
using OrdinaryDiffEqLowOrderRK: OrdinaryDiffEqLowOrderRK
using SciMLBase: DiscreteCallback, terminate!

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

"""

$(TYPEDSIGNATURES)

Returns a function that evaluates the RHS of the dimensionless TOV equations, given an equation of state function.

# Arguments
- `EoS_rho_from_P::Function`: Single-argument (inverse) equation-of-state function that maps dimensionless density (``\\rho``) to dimensionless pressure (``P``).

# Returns
- `dimensionless_TOV_inner!::Function`: Callable as `dimensionless_TOV_inner!(du, u, params, r)` that evaluates the RHS of the dimensionless TOV equations, writing the result to `du`.
"""
function dimensionless_TOV(EoS_rho_from_P::Function)
    function dimensionless_TOV_inner!(du, u, params, r)
        P = u[1]
        m = u[2]
        rho = EoS_rho_from_P(P)
        du[1] = if r == 0.0
            0.0
        else
            - (rho + P) * (m + 4*pi*r^3*P) / (r * (r - 2*m))
        end
        du[2] = 3*pi*r^2*rho
    end
    return dimensionless_TOV_inner!
end

"""

$(TYPEDSIGNATURES)

Return a `NamedTuple` containing the characteristic length scales for each of the fields in the TOV equation.

The TOV equation can be made dimensionless by using the following characteristic length scales;

```math
r = R \\hat{r}, \\quad
m = \\frac{GR}{c^2} \\hat{m}, \\quad
P = \\frac{M c^2}{R^3} \\hat{P}, \\quad
\\rho = \\frac{M}{R^3} \\hat{\\rho},
```

where ``\\hat{\\cdot}`` denotes the corresponding dimensionless field, ``G`` is the gravitational constant, and ``c`` the speed of light in a vacuum. The characteristic lengths (coefficients of the dimensionless fields, above) have only one degree of freedom, so when one of them is specified the others can be computed.

This function, given the value for any **one** of the characteristic lengths, returns a `NamedTuple` that specifies the values of all the characteristic lengths. Passing in multiple characteristic lengths, or none, will result in an error (even in the case when the characteristic lengths that are passed are consistent with each other).

The values for ``c`` and ``G`` may also be specified to the function. This should be done when wanting to use a non-SI-unit system, for example CGS or ``c = 1`` units.

# Arguments
- `length::Number`: Characteristic length for the space dimension, equal to ``R``.
- `mass::Number`: Characteristic length for the mass dimension, equal to ``\\frac{GR}{c^2}``.
- `pressure::Number`: Characteristic length for the pressure field, equal to ``\\frac{G}{R^2}``.
- `density::Number`: Characteristic length for the density field, equal to ``\\frac{G}{c^2 R^2}``.
- `c::Number`: Speed of light in vacuum. Defaults to value in SI units, ``\\approx 3\times 10^8``.
- `G::Number`: Gravitational constant. Defaults to value in SI units, ``\\approx 6.7\times 10^{-11}``.

# Returns
- `::NamedTuple`: Map of characteristic lengths `R`, `M`, `Q`, `Rho` to their values.
"""
function characteristic_lengths_TOV(;
    length::Union{Number,Nothing} = nothing,
    mass::Union{Number,Nothing} = nothing,
    pressure::Union{Number,Nothing} = nothing,
    density::Union{Number,Nothing} = nothing,
    c::Number = speed_of_light_vacuum,
    G::Number = gravitational_constant,
)
    RMQRho = [length, mass, pressure, density];
    RMQRho_given = .!isnothing.(RMQRho);

    if sum(RMQRho_given) != 1
        throw("Multiple, or no, length scales provided.")
    end

    # Compute the value of R, if it was not the scale given.
    if RMQRho_given[2]
        RMQRho[1] = RMQRho[2] * c^2 / G
    elseif RMQRho_given[3]
        RMQRho[1] = sqrt(G / RMQRho[3])
    elseif RMQRho_given[4]
        RMQRho[1] = sqrt(G / RMQRho[4]) / c
    end

    # Now that R is necessarily defined, compute the missing length scales.
    if !RMQRho_given[2]
        RMQRho[2] = G * RMQRho[1] / c^2
    end
    if !RMQRho_given[3]
        RMQRho[3] = G / RMQRho[1]^2
    end
    if !RMQRho_given[4]
        RMQRho[4] = G / (c * RMQRho[1])^2
    end

    return (R = RMQRho[1], M = RMQRho[2], Q = RMQRho[3], Rho = RMQRho[4])
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

end
