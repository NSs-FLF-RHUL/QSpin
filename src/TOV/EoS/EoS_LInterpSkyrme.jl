using DocStringExtensions: TYPEDSIGNATURES
using ....Parameters: ParameterType
using Roots: find_zero
using CSV: File
using DataFrames: DataFrame, select!
using DataInterpolations: ExtrapolationType, PCHIPInterpolation
"""
$(TYPEDSIGNATURES)
Loading the data from a pre-computed equation of state (EoS) and using PCHIPInterpolation interpolation to get the function of equation of state for TOV equation computations.

# Arguments
- `file_input::Union{DataFrame,AbstractMatrix{Float64}}`: the file directory and file name of the pre-computed equation of state in dat format or lading them in an N-by-x array with a correct EoS_indices.
- `EoS_indices::tuple`: the indices in the dat for density - rho, and pressure - press, baryonic number density - nb, neutron density - rhoN, proton fraction - Yp, effecive proton mass - mp, and fermi momenta for electrons and neutrons - kFe and kFn respectively.

# Returns
- 'EoS': The NamedTuple contains EoS, EoS_inv, EoS_nb, EoS_rhoN, EoS_Yp, EoS_mp, EoS_kFe and EoS_kFn are in the format of ::DataInterpolations, and that maps ρ to different quantities listed above and EoS_inv is the inverse of EoS, that ρ as a function of press. ρ_cci is a Floating number about the core-crust interface density.
"""
function EoS_LInterpSkyrme(
    input::Union{DataFrame,AbstractMatrix{Float64}};
    EoS_indices::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64} = (
        1,
        2,
        3,
        4,
        9,
        10,
        14,
        16,
    ),
)
    i_rho, i_press, i_nb, i_rhoN, i_Yp, i_mp_ast, i_kFe, i_kFn = EoS_indices
    input = Array(input)
    input = sortslices(input, dims = i_rho)
    rho = input[:, i_rho]
    if minimum(rho) > 0.0
        @info "The minimum density in the input file is greater than zero. Adding a zero density point to the input data."
        rho = [0; input[:, i_rho]]
        press = [0; input[:, i_press]]
        nb = [0; input[:, i_nb]]
        rhoN = [0; input[:, i_rhoN]]
        Yp = [0; input[:, i_Yp]] # Do we really need this? It seems not used in the mutual fricition evaluations.
        mp_ast = [0; input[:, i_mp_ast]]
        kFe = [0; input[:, i_kFe]]
        kFn = [0; input[:, i_kFn]]
    elseif minimum(rho) < 0.0
        error("The density in the input file must be non-negative.")
    else
        press = input[:, i_press]
        nb = input[:, i_nb]
        rhoN = input[:, i_rhoN]
        Yp = input[:, i_Yp] # Do we really need this? It seems not used in the mutual fricition evaluations.
        mp_ast = input[:, i_mp_ast]
        kFe = input[:, i_kFe]
        kFn = input[:, i_kFn]
    end


    EoS_P_from_rho =
        PCHIPInterpolation(press, rho; extrapolation = ExtrapolationType.Extension)
    EoS_rho_from_P =
        PCHIPInterpolation(rho, press; extrapolation = ExtrapolationType.Extension)
    EoS_nb_from_rho =
        PCHIPInterpolation(nb, rho; extrapolation = ExtrapolationType.Extension)
    EoS_rhoN_from_rho =
        PCHIPInterpolation(rhoN, rho; extrapolation = ExtrapolationType.Extension)
    EoS_Yp_from_rho =
        PCHIPInterpolation(Yp, rho; extrapolation = ExtrapolationType.Extension)
    EoS_mp_ast_from_rho =
        PCHIPInterpolation(mp_ast, rho; extrapolation = ExtrapolationType.Extension)
    EoS_kFe_from_rho =
        PCHIPInterpolation(kFe, rho; extrapolation = ExtrapolationType.Extension)
    EoS_kFn_from_rho =
        PCHIPInterpolation(kFn, rho; extrapolation = ExtrapolationType.Extension)

    imp_zero = findall(iszero, EoS_mp_ast_from_rho.u)
    inb_cci = findlast(x->x<0.2, EoS_nb_from_rho.u[imp_zero]*1e-39)
    ρ_cci = EoS_rho_from_P.u[imp_zero[inb_cci]]

    return (
        EoS = EoS_P_from_rho,
        EoS_inv = EoS_rho_from_P,
        EoS_nb = EoS_nb_from_rho,
        EoS_rhoN = EoS_rhoN_from_rho,
        EoS_Yp = EoS_Yp_from_rho,
        EoS_mp = EoS_mp_ast_from_rho,
        EoS_kFe = EoS_kFe_from_rho,
        EoS_kFn = EoS_kFn_from_rho,
        ρ_cci = ρ_cci,
    )
end

function EoS_LInterpSkyrme(
    file_input::String;
    EoS_indices::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64} = (
        1,
        2,
        3,
        4,
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
