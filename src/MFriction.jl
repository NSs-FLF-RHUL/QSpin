module MFriction

using JSON: JSON
using DataInterpolations: ExtrapolationType, QuadraticSpline
using DocStringExtensions: TYPEDSIGNATURES
using ..PhysicalConstants: hbar, neutron_mass, electron_volt

"""
$(TYPEDSIGNATURES)

This function reads in JSON data for a given file path and calculates the mutual friction parameters for a neutron star crust based on the Graber et al. 2018 model.
According to Graber et al. 2018, the mutual friction coefficients are calculated based on the superfluid density and other physical parameters.
The function returns a tuple containing the input parameters and the calculated mutual friction parameters, including the qubic spline interpolations for the mutual friction coefficients as functions of the superfluid density (in kg * m^-3, while the coverted input is in kg fm^-3) in log-log space.

# Arguments
- `file_path`: A string representing the path to the JSON file containing the input parameters for the mutual friction calculations. The JSON file should contain an array of objects, each representing a different region of the neutron star crust with specific parameters such as baryon number density (nb), proton number (Z), neutron number (N), proton fraction (x), superfluid density (ns), lattice spacing (a), nuclear radius (RN), and pinning energy parameters (Es, E1, DE, xi, Ep).

# Returns
- 'output': A tuple containing the input parameters (in their original units from the input JSON file) and the calculated mutual friction parameters in array forms. The qubic spline interpolations for the mutual friction coefficients, B_EW and B_J, as functions of the superfluid density (in kg * m^-3, while the coverted input is in kg fm^-3) are included.
"""
function VNparaGraber2018(file_path)

    MeV = electron_volt * 1e6 # convert to kg * fm^2 / s^2
    δ = 1e-2 # dimensionless coefficient for the pinning energy reduction due to vortex tension
    κ = hbar / neutron_mass * π  # SI units

    data = JSON.parsefile(file_path)
    println("Successfully loaded JSON data!\n")
    nb = zeros(length(data)-1)
    Z = zeros(length(data)-1)
    N = zeros(length(data)-1)
    x = zeros(length(data)-1)
    ns = zeros(length(data)-1)
    a = zeros(length(data)-1)
    Rn = zeros(length(data)-1)
    Es = zeros(length(data)-1)
    E1 = zeros(length(data)-1)
    ΔE = zeros(length(data)-1)
    ξ = zeros(length(data)-1)
    Ep = zeros(length(data)-1)
    #    A = zeros(length(data)-1)
    #    RWS = zeros(length(data)-1)
    for dd = 2:length(data)
        nb[dd-1] = data[dd].nb
        Z[dd-1] = data[dd].Z
        N[dd-1] = data[dd].N
        x[dd-1] = data[dd].x
        ns[dd-1] = data[dd].ns
        a[dd-1] = data[dd].a
        Rn[dd-1] = data[dd].RN
        Es[dd-1] = data[dd].Es
        E1[dd-1] = data[dd].E1
        ΔE[dd-1] = data[dd].DE
        ξ[dd-1] = data[dd].xi
        Ep[dd-1] = data[dd].Ep
        #A[dd-1]   = data[dd].Z * (1 + 1/data[dd].x)
        #Rws[dd-1] = (3*(data[dd].N+data[dd].Z) / (4 * π * data[dd].nb*1e-4))^(1/3)
    end


    A = Z .* (1 .+ 1 ./ x)
    Rws = (3*(N .+ Z) ./ (4 * π * nb * 1e-4)) .^ (1/3)
    n1 = 3/4/π ./ Rws .^ 3*1e6
    ρs = ns * 1e-4 * neutron_mass * 1e45# in kg * m^-3
    Reb =
        2.8 * sqrt.(0.5 * neutron_mass/hbar) * sqrt.(abs.(Ep * MeV) * δ ./ ρs / κ) .* Rn ./
        a .^ (1.5) * 10^(7.5)
    Rj =
        0.5 / sqrt(π) *
        (0.5 * neutron_mass/hbar)^(1/2) *
        (abs.(Ep * MeV) * δ ./ ρs / κ) .^ (1/2) .* a .^ 0.5 ./ ξ * 10^(7.5)
    Beb = Reb ./ (1 .+ Reb .^ 2)
    Bj = Rj ./ (1 .+ Rj .^ 2)
    Beb_itp = QuadraticSpline(
        log10.(Beb),
        log10.(ρs);
        extrapolation_right = ExtrapolationType.Extension,
        extrapolation_left = ExtrapolationType.Constant,
    )
    Bj_itp = QuadraticSpline(
        log10.(Bj),
        log10.(ρs);
        extrapolation_right = ExtrapolationType.Extension,
        extrapolation_left = ExtrapolationType.Constant,
    )
    output = (;
        nb,
        Z,
        N,
        x,
        ns,
        n1,
        a,
        Rn,
        Es,
        E1,
        ΔE,
        ξ,
        Ep,
        A,
        Rws,
        ρs,
        Reb,
        Rj,
        Beb,
        Bj,
        Beb_itp,
        Bj_itp,
    )
    return output
end

function MutualFrictionCoefficients(Param, Beb_itp, Bj_itp; ρ_drip = 4e11*1e3, Rcci = 1e4)
    log_ρs = log10.(Param.ρs)
    Beb = exp10.(Beb_itp.(log_ρs))
    Bj = exp10.(Bj_itp.(log_ρs))
    Beb[Param.ρs .< ρ_drip] .= 0.0
    Bj[Param.ρs .< ρ_drip] .= 0.0
    Beb[Param.r .< Rcci] .= Param.Beb_core
    Bj[Param.r .< Rcci] .= Param.Bj_core

    Bs = (; Beb, Bj)
    return Bs

end

end
