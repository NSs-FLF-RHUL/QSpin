using DocStringExtensions: TYPEDSIGNATURES
using ....PhysicalConstants: electron_volt, hbar, neutron_mass, speed_of_light_vacuum
using DataInterpolations: PCHIPInterpolation
using Roots: find_zero, Bisection
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
    ρ_drip = 4.e11, # in g/cm^3
    Ye = 0.4,
)
    ħ = hbar * 1e3 * 1e4 # convert to g * cm^2 / s
    mn = neutron_mass * 1e3 # convert to g
    c = speed_of_light_vacuum * 1e2 # convert to cm / s
    function EoS_P_from_rho(ρ)
        P = if ρ < 0
            0.0
        elseif ρ < ρ_drip
            ħ * c * (3 * π^2 * Ye * ρ / mn) .^ (4/3) / 12 / π^2
        else
            nb = ρ / (neutron_mass * 1e3) # in the unit of g/cm^3
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

        return @. P
    end

    P_drip = EoS_P_from_rho(ρ_drip)
    N_intp = 100
    ρ_beg = ρ_drip
    ρ_end = 1e16
    ρ_intp = exp.(range(log(ρ_beg), log(ρ_end), length = N_intp))
    P_intp = EoS_P_from_rho.(ρ_intp)
    EoS_rho_from_P_intp = PCHIPInterpolation(ρ_intp, P_intp)


    function EoS_rho_from_P(P)
        ρ = if P <= 0
            0.0
        elseif isapprox(P, P_drip; rtol = 1e-12)
            ρ_drip
        elseif P < P_drip
            (12 * π^2 * P / ħ / c)^(3/4) * mn / (3 * π^2 * Ye)
        elseif P > P_drip
            rho_guess = EoS_rho_from_P_intp.(P)
            func = y -> EoS_P_from_rho(y) - P

            try
                find_zero(
                    func,
                    (
                        maximum([ρ_drip, rho_guess*0.95]),
                        maximum([1.05*ρ_drip, 1.05*rho_guess]),
                    ),
                    Bisection(),
                )
            catch eval
                @show P rho_guess EoS_P_from_rho(rho_guess)
                rethrow(eval)
            end



        end
        return ρ
    end
    return EoS_P_from_rho, EoS_rho_from_P
end
