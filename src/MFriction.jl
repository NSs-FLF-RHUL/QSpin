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

Here things are converted in the SI units so as MutualFrictionCoefficients

# Arguments
- `file_path`: A string representing the path to the JSON file containing the input parameters for the mutual friction calculations. The JSON file should contain an array of objects, each representing a different region of the neutron star crust with specific parameters such as baryon number density (nb), proton number (Z), neutron number (N), proton fraction (x), superfluid density (ns), lattice spacing (a), nuclear radius (RN), and pinning energy parameters (Es, E1, DE, xi, Ep).

# Returns
- 'output': A tuple containing the input parameters (in their original units from the input JSON file) and the calculated mutual friction parameters in array forms. The qubic spline interpolations for the mutual friction coefficients, B_EW and B_J, as functions of the superfluid density (in kg * m^-3, while the coverted input is in kg fm^-3) are included.
"""
function VNparaGraber2018(file_path; units = "CGS")

    MeV = electron_volt * 1e6 # convert to kg * fm^2 / s^2
    δ = 1e-2 # dimensionless coefficient for the pinning energy reduction due to vortex tension
    κ = hbar / neutron_mass * π  # SI units

    data = JSON.parsefile(file_path)
    @info "JSON data successfully loaded from $(file_path)"
    nb = zeros(length(data)-1)
    Z = zeros(length(data)-1)
    N = zeros(length(data)-1)
    x = zeros(length(data)-1)
    ns = zeros(length(data)-1)
    ρ = zeros(length(data)-1) # in g * cm^-3
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
        row = data[dd]
        nb[dd-1] = row["nb"]
        Z[dd-1] = row["Z"]
        N[dd-1] = row["N"]
        x[dd-1] = row["x"]
        ns[dd-1] = row["ns"]
        ρ[dd-1] = row["ρ"]
        a[dd-1] = row["a"]
        Rn[dd-1] = row["RN"]
        Es[dd-1] = row["Es"]
        E1[dd-1] = row["E1"]
        ΔE[dd-1] = row["DE"]
        ξ[dd-1] = row["xi"]
        Ep[dd-1] = row["Ep"]
        #A[dd-1]   = data[dd].Z * (1 + 1/data[dd].x)
        #Rws[dd-1] = (3*(data[dd].N+data[dd].Z) / (4 * π * data[dd].nb*1e-4))^(1/3)
    end


    A = Z .* (1 .+ 1 ./ x)
    Rws = (3*(N .+ Z) ./ (4 * π * nb * 1e-4)) .^ (1/3)
    n1 = 3/4/π ./ Rws .^ 3*1e6
    ρs = ns * 1e-4 * neutron_mass * 1e45# in kg * m^-3
    ρ = ρ * 1e12 * 1e3 # in kg * m^-3
    EpA = sqrt.(Es .^ 2 + Es .* E1 + 0.5 * E1 .^ 2)
    RA =
        2.8 * sqrt.(0.5 * neutron_mass/hbar) * sqrt.(abs.(EpA * MeV) * δ ./ ρs / κ) .* Rn ./
        a .^ (1.5) * 10^(7.5)
    Reb =
        2.8 * sqrt.(0.5 * neutron_mass/hbar) * sqrt.(abs.(Ep * MeV) * δ ./ ρs / κ) .* Rn ./
        a .^ (1.5) * 10^(7.5)
    Rj =
        0.5 / sqrt(π) *
        (0.5 * neutron_mass/hbar)^(1/2) *
        (abs.(Ep * MeV) * δ ./ ρs / κ) .^ (1/2) .* a .^ 0.5 ./ ξ * 10^(7.5)
    BA = RA ./ (1 .+ RA .^ 2)
    Beb = Reb ./ (1 .+ Reb .^ 2)
    Bj = Rj ./ (1 .+ Rj .^ 2)
    BA_itp = dlog10_fit(
        QuadraticSpline(
            log10.(BA),
            log10.(ρs);
            extrapolation_right = ExtrapolationType.Extension,
            extrapolation_left = ExtrapolationType.Constant,
        ),
    )
    Beb_itp = dlog10_fit(
        QuadraticSpline(
            log10.(Beb),
            log10.(ρs);
            extrapolation_right = ExtrapolationType.Extension,
            extrapolation_left = ExtrapolationType.Constant,
        ),
    )
    Bj_itp = dlog10_fit(
        QuadraticSpline(
            log10.(Bj),
            log10.(ρs);
            extrapolation_right = ExtrapolationType.Extension,
            extrapolation_left = ExtrapolationType.Constant,
        ),
    )
    R = (RA, Reb, Rj)
    B = (BA, Beb, Bj)
    B_itp = (BA_itp, Beb_itp, Bj_itp)
    output = (; nb, Z, N, x, ns, ρ, n1, a, Rn, Es, E1, ΔE, ξ, Ep, A, Rws, ρs, R, B, B_itp)
    return output
end
"""
$(TYPEDSIGNATURES)
# Arguments

Evaluate the mutual friction coefficent for an input DataInterpolation that fits the data points in double-log10 space.

# Arguments
- `B_log10_intepr`: The interpolation function for B in log-log space

# Returns
- `dlog10_fit_inner: Recovering the results back to linear space.
"""
function dlog10_fit(B_log10_intepr)
    function dlog10_fit_inner(ρ::Union{Float64,AbstractArray{Float64}})
        exp10.(B_log10_intepr(log10.(ρ)))
    end
    return dlog10_fit_inner
end

"""
    MutualFrictionCoefficients(Param, Beb_itp, Bj_itp; input_units="SI",
                               ρ_drip=nothing, Rcci=nothing)

Evaluate the electron-vortex and Jones mutual-friction coefficients along a
density/radius profile. `Param.ρs`, `Param.r`, and any explicit `ρ_drip` or
`Rcci` values use `input_units`. Supported systems are `"SI"` (kg/m³, m) and
`"CGS"` (g/cm³, cm). The interpolation tables are always evaluated in kg/m³.

The defaults represent `ρ_b = ρ_drip = 4e14 kg/m³` and `Rcci = 10 km`.
"""
function MutualFrictionCoefficients(
    Param,
    B_itp;
    input_units::String = "SI",
    ρ_b = nothing,
    R_cci = nothing,
)
    if input_units == "SI"
        density_to_si = 1.0
        ρ_b = something(ρ_b, 4e14)
        R_cci = something(R_cci, 1e4)
    elseif input_units == "CGS"
        density_to_si = 1e3
        ρ_b = something(ρ_b, 4e11)
        R_cci = something(R_cci, 1e6)
    else
        throw(ArgumentError("input_units must be \"SI\" or \"CGS\""))
    end

    BA = B_itp[1](Param.ρs .* density_to_si)
    Beb = B_itp[2](Param.ρs .* density_to_si) #exp10.(Beb_itp.(log_ρs))
    Bj = B_itp[3](Param.ρs .* density_to_si) #exp10.(Bj_itp.(log_ρs))
    BA[Param.ρs .< ρ_b] .= B_itp[1](ρ_b)
    Beb[Param.ρs .< ρ_b] .= B_itp[2](ρ_b)
    Bj[Param.ρs .< ρ_b] .= B_itp[3](ρ_b)
    Bs = (BA, Beb, Bj)
    return Bs

end

end
