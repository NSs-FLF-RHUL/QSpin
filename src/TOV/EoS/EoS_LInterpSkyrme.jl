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
    input = Array(input)
    input = sortslices(input, dims = i_rho)
    rho = input[:, i_rho]
    press = input[:, i_press]
    nb = input[:, i_nb]
    Yp = input[:, i_Yp] # Do we really need this? It seems not used in the mutual fricition evaluations.
    mp_ast = input[:, i_mp_ast]
    kFe = input[:, i_kFe]
    kFn = input[:, i_kFn]
    EoS_P_from_rho =
        QuadraticSpline(press, rho; extrapolation = ExtrapolationType.Extension)
    EoS_rho_from_P =
        QuadraticSpline(rho, press; extrapolation = ExtrapolationType.Extension)
    EoS_nb_from_rho = QuadraticSpline(nb, rho; extrapolation = ExtrapolationType.Extension)
    EoS_Yp_from_rho = QuadraticSpline(Yp, rho; extrapolation = ExtrapolationType.Extension)
    EoS_mp_ast_from_rho =
        QuadraticSpline(mp_ast, rho; extrapolation = ExtrapolationType.Extension)
    EoS_kFe_from_rho =
        QuadraticSpline(kFe, rho; extrapolation = ExtrapolationType.Extension)
    EoS_kFn_from_rho =
        QuadraticSpline(kFn, rho; extrapolation = ExtrapolationType.Extension)
    return EoS_P_from_rho,
    EoS_rho_from_P,
    EoS_nb_from_rho,
    EoS_Yp_from_rho,
    EoS_mp_ast_from_rho,
    EoS_kFe_from_rho,
    EoS_kFn_from_rho
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
