using DocStringExtensions: TYPEDSIGNATURES
using Roots: find_zero
using ....Parameters: ParameterType
using ....PhysicalConstants: electron_volt, neutron_mass

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
            rho_guess = if log(P) > 76.0
                # ((log(P) - 75.06) / 4.87e-9)^(1/0.6028) # Emperical Guess from a Fitting for density up to about 2e16
                ((log(P) - 68.3) / 2.836e-7)^(1/0.502) # Emperical Guess from a Fitting for density up to about 4e16
            else
                5e12
            end
            find_zero(y -> EoS_P_from_rho(y) - P, rho_guess)
        end
        return ρ
    end
    return EoS_P_from_rho, EoS_rho_from_P
end
