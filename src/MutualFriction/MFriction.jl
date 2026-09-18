module MFriction

using JSON: JSON
using DataInterpolations: ExtrapolationType, QuadraticSpline
using DocStringExtensions: TYPEDSIGNATURES
using ..PhysicalConstants: hbar, neutron_mass, electron_volt
using ..Parameters: ParameterType

include("mfrictionGraber2016.jl")
include("mfrictionGraber2018.jl")


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

    BA = B_itp[1].(Param.ρ .* density_to_si)
    Beb = B_itp[2].(Param.ρ .* density_to_si) #exp10.(Beb_itp.(log_ρs))
    Bj = B_itp[3].(Param.ρ .* density_to_si) #exp10.(Bj_itp.(log_ρs))
    BA[Param.ρ .< ρ_b] .= B_itp[1].(ρ_b * density_to_si)
    Beb[Param.ρ .< ρ_b] .= B_itp[2].(ρ_b * density_to_si)
    Bj[Param.ρ .< ρ_b] .= B_itp[3].(ρ_b * density_to_si)
    Bs = (; BA, Beb, Bj)
    return Bs

end

end
