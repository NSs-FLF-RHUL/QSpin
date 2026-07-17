using DocStringExtensions: TYPEDSIGNATURES
using ....Parameters: ParameterType
using Roots: find_zero
using CSV: File
using DataFrames: DataFrame, select!
using DataInterpolations: ExtrapolationType, QuadraticSpline

"""
$(TYPEDSIGNATURES)
Loading the data from a pre-computed equation of state (EoS) and using quadraticspline interpolation to get the function of equation of state for TOV equation computations.

# Arguments
- `file_input::`: the file directory and file name of the pre-computed equation of state in dat format or lading them in an N-by-x array with a correct EoS_indices.
- `EoS_indices::tuple`: the indices in the dat for density, rho, and pressure, press.

# Returns
- `EoS_P_from_rho`: The interpolation for the equation of state.
- `EoS_rho_from_P`: The interpolation for the inverse equation of state.
"""
function EoS_LInterp(
    file_input::Union{String,AbstractArray,Matrix{Float64}},
    EoS_indices::Tuple{Int64,Int64},
)

    if typeof(file_input) == String
        df = DataFrame(File(file_input, delim = " "))
        select!(df, [k for (k, v) in pairs(eachcol(df)) if !all(ismissing, v)])
    elseif typeof(file_input) == AbstractArray || typeof(file_input) == Matrix{Float64}
        df = file_input
    else
        error(
            "Unsupported input -- file_input must be an N-by-2 array or a string for the directory of an .dat input.",
        )
    end
    i_rho, i_press = EoS_indices
    rho = sort(df[:, i_rho])
    press = sort(df[:, i_press])
    EoS_P_from_rho =
        QuadraticSpline(press, rho; extrapolation = ExtrapolationType.Extension);
    EoS_rho_from_P =
        QuadraticSpline(rho, press; extrapolation = ExtrapolationType.Extension);
    return EoS_P_from_rho, EoS_rho_from_P#, rho, press
end
