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
function EoS_LInterpSkyrme(
    input::Union{DataFrame,AbstractMatrix{Float64}};
    EoS_indices::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64} = (
        1,
        2,
        3,
        9,
        10,
        14,
        16,
    ),
)
    i_rho, i_press, i_nb, i_Yp, i_mp_ast, i_kFe, i_kFn = EoS_indices
    rho = sort(input[:, i_rho])
    press = sort(input[:, i_press])
    nb = sort(input[:, i_nb])
    Yp = sort(input[:, i_Yp])
    mp_ast = sort(input[:, i_mp_ast])
    kFe = sort(input[:, i_kFe])
    kFn = sort(input[:, i_kFn])
    EoS_P_from_rho =
        QuadraticSpline(press, rho; extrapolation = ExtrapolationType.Extension)
    EoS_rho_from_P =
        QuadraticSpline(rho, press; extrapolation = ExtrapolationType.Extension)
    EoS_rho_from_nb = QuadraticSpline(rho, nb; extrapolation = ExtrapolationType.Extension)
    EoS_rho_from_Yp = QuadraticSpline(rho, Yp; extrapolation = ExtrapolationType.Extension)
    EoS_rho_from_mp_ast =
        QuadraticSpline(rho, mp_ast; extrapolation = ExtrapolationType.Extension)
    EoS_rho_from_kFe =
        QuadraticSpline(rho, kFe; extrapolation = ExtrapolationType.Extension)
    EoS_rho_from_kFn =
        QuadraticSpline(rho, kFn; extrapolation = ExtrapolationType.Extension)
    return EoS_P_from_rho,
    EoS_rho_from_P,
    EoS_rho_from_nb,
    EoS_rho_from_Yp,
    EoS_rho_from_mp_ast,
    EoS_rho_from_kFe,
    EoS_rho_from_kFn
end

function EoS_LInterpSkyrme(
    file_input::String;
    EoS_indices::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64} = (
        1,
        2,
        3,
        9,
        10,
        14,
        16,
    ),
)
    df = DataFrame(File(file_input, delim = " "))
    select!(df, [k for (k, v) in pairs(eachcol(df)) if !all(ismissing, v)])
    return EoS_LInterpSkyrme(df; EoS_indices = EoS_indices)
end
