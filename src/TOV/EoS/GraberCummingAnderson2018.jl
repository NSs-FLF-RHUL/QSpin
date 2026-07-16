using DocStringExtensions: TYPEDSIGNATURES
using ....Parameters: ParameterType
using ....PhysicalConstants: electron_volt, hbar, neutron_mass, speed_of_light_vacuum
using Roots: find_zero
"""
$(TYPEDSIGNATURES)

Create EoS functions according to V. Graber, A. Cumming & N. Anderson, ApJ 865, 23 (2018).

The EoS given by J. Negele & D. Vautherin, Nuclear Physics A 207, 298 (1973) is used for the inner crust regime and the outer crust is considered by

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
    ρ_drip = 4.e11, # in g/m^3,
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
            rho_guess = if log(P) > 76.0
                # ((log(P) - 75.06) / 4.87e-9)^(1/0.6028) # Emperical Guess from a Fitting for density up to about 4e16
                ((log(P) - 68.3) / 2.836e-7)^(1/0.502) # Emperical Guess from a Fitting for density up to about 9e16
            else
                5e12
            end
            find_zero(y -> EoS_P_from_rho(y) - P, rho_guess)
        end
        return ρ
    end
    return EoS_P_from_rho, EoS_rho_from_P
end
