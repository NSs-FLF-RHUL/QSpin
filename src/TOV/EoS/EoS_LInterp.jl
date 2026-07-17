using DocStringExtensions: TYPEDSIGNATURES
using ....Parameters: ParameterType
using Roots: find_zero
using CSV
using DataFrames
using DataInterpolations: ExtrapolationType, QuadraticSpline

function EoS_LInterp(file_name, EoS_indices)
    df = DataFrame(CSV.File(file_name, delim = " "))
    df = select!(df, [k for (k, v) in pairs(eachcol(df)) if !all(ismissing, v)])
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
