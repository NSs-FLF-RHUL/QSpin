using DocStringExtensions: TYPEDSIGNATURES
using ....Parameters: ParameterType
using Roots: find_zero
using CSV
using DataFrames
using DataInterpolations: ExtrapolationType, QuadraticSpline

"""
$(TYPEDSIGNATURES)
Loading the data from a pre-computed equation of state (EoS) and using quadraticspline interpolation to get the function of equation of state for TOV equation computations.

# Arguments
- `file_name::`: the file directory and file name of the pre-computed equation of state in dat format
- `EoS_indices::tuple`: the indices in the dat for density, rho, and pressure, press.

# Returns
- `EoS_P_from_rho`: The interpolation for the equation of state.
- `EoS_rho_from_P`: The interpolation for the inverse equation of state.
"""
function EoS_LInterp(file_name::String, EoS_indices::Tuple{Int64,Int64})
    df = DataFrame(CSV.File(file_name, delim = " "))
    select!(df, [k for (k, v) in pairs(eachcol(df)) if !all(ismissing, v)])
    #deleteat!(permutedims(df), all.(ismissing, eachcol(df))) |> permutedims # Seems to be equivalent.
    i_rho, i_press = EoS_indices
    rho = sort(df[:, i_rho])
    press = sort(df[:, i_press])
    EoS_P_from_rho =
        QuadraticSpline(press, rho; extrapolation = ExtrapolationType.Extension);
    EoS_rho_from_P =
        QuadraticSpline(rho, press; extrapolation = ExtrapolationType.Extension);
    return EoS_P_from_rho, EoS_rho_from_P#, rho, press
end
