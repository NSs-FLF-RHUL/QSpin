"""
Submodule containing helper functions for solving the TOV equation.

The TOV equation reads

``
\\frac{\\mathrm{d}P}{\\mathrm{d}r} &=
- \\frac{G}{r^2}
\\left[ \\rho(r) + \\frac{P(r)}{c^2} \\right]
\\left[ m(r) + 4\\pi r^3 \frac{P(r)}{c^2} \\right]
\\left[ 1 - \frac{2Gm(r)}{c^2 r} \\right]^{-1},
``

where ``G`` and ``c`` represent the gravitational constant and speed of light in vacuum, respectively.
``r`` is the radial coordinate, with ``\\rho(r)`` and ``P(r)`` the mass density and pressure as position ``r``, respectively.
Finally, ``m(r)`` is the gravitational mass satisfying ``m(0) = 0`` and is found by solving the continuity equation

``
\\frac{\\mathrm{d}m}{\\mathrm{d}r} &= 4\\pi\\rho(r)r^2.
``

To solve the TOV equation we must also provide an equation of state (EoS) that relates the pressure and density of the neutron star, that is we must provide ``P = P(\\rho)`` (and the inverse).

Note that we can solve the TOV equation(s) numerically using the vector ``u = (P, m)`` and solving a first-order system in ``u``.
"""

using ..PhysicalConstants: gravitational_constant, speed_of_light_vacuum

"""
Return a function that evaluates the TOV equation.

Returns a function `tov!(du, u, params, r)` that evaluates the TOV equation (in vector form) given the EoS relationship ``\\rho(P)``.

`Eos_inv` is assumed to take `(P, parameters, r)` as inputs, where `P = u[1]` is the pressure.
"""
function tov_eq!(EoS_inv::Function)
    function tov_inner!(du, u, params, r)
        P = u[1]
        m = u[2]
        rho = EoS_inv(P);
        if r == 0.0
            du[1] .= 0.0;
        else
            du[1] .=
                - gravitational_constant / r^2 *
                (rho + P / speed_of_light_vacuum^2) *
                (m + 4*π*r^3*P / speed_of_light_vacuum^2) /
                (1 - 2*gravitational_constant*m / (r*speed_of_light_vacuum^2));
        end
        du[2] .= 4*π*r^2*rho;
    end
    return tov_inner!
end
