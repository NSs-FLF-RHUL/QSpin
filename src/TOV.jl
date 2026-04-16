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

using ..Parameters: ParameterType
using ..PhysicalConstants: hbar, gravitational_constant, neutron_mass, speed_of_light_vacuum

function tov_eq!(EoS_rho_from_P::Function)
    function tov_inner!(du,u,paras,r)
        P = u[1];
        m = u[2];
        ρ = EoS_rho_from_P(P);
        if r == 0.0
            du[1] = 0.;
        else
            du[1] =
                -gravitational_constant / r^2 *
                (ρ + P/speed_of_light_vacuum^2) *
                (m + 4*pi*r^3*P/speed_of_light_vacuum^2) / (1 - 2*gravitational_constant*m/(r*speed_of_light_vacuum^2));
        end
        du[2] = 4*pi*r^2*ρ;
    end
    return tov_inner!
end
