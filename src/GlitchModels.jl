module GlitchModels
using DocStringExtensions: TYPEDSIGNATURES
using ..Parameters: ParameterType

"""
$(TYPEDSIGNATURES)
"""
function ThreeCompSolid!(
    dΩ::AbstractArray,
    Ω::AbstractArray,
    Param::ParameterType,
    time::Float64,
)
    Ω_sf = Ω[2];
    Ω_core = Ω[3];
    Ω_crust = Ω[1];
    dΩ[2] = Param.B_sf * (Ω_crust - Ω_sf);
    dΩ[3] = Param.B_core * (Ω_crust - Ω_core);
    dΩ[1] =
        -Param.N_ext / Param.I_crust - Param.I_sf / Param.I_crust * dΩ[2] -
        Param.I_core / Param.I_crust * dΩ[3];
end

function ThreeCompGCA2018!(EoMSetup::ParameterType; value_check::Bool = false)

    Nr = length(EoMSetup.r);
    I_total = 0.35 * EoMSetup.M_NS * EOMSetup.R_NS^2
    dr = [diff(EoMSetup.r); diff(EoMSetup.r)[end]]
    spherical_momentum_of_innertia = EoMSetup.ρr .* EoMSetup.r .^ 4 .* dr
    unit_momenum_of_innertia = EvoMSetup.ρr .* EoMSetup.r .^ 3 .* dr
    I_crust_total = 8 * π * sum(spherical_momentum_of_innertia[i_Rcci:end])
    I_crust_unitheight = 2 * π * sum(nit_momenum_of_innertia[i_Rcci:end])
    h_crust = I_crust_unitheight / I_crust_total / 2.0

    I_core = EoMSetup.α_core * (I_total - I_crust_total)
    I_proton = EoMSetup.α_proton * (I_total - I_crust_total)
    i_R_drip
    I_sf = 4 * π * EvoMSetup.ρr .* EvoMSetup.r .^ 3 .* dr;
    I_sf = sum(I_sf(i_Rcci:iR_drip)) * h_crust
    if value_check
        println("")
    end

    function ThreeComMod_inner!(
        dΩ::AbstractArray,
        Ω::AbstractArray,
        Param::ParameterType,
        time::Float64,
    )
        # dΩ_sf/dt
        Ω_sf = Ω[2:(Nr+1)];
        #Bsf = Param.B_sf * Param.ρr ./ maximum(Param.ρr); # Scaling B_sf with the local density profile
        dΩ_sfdr = [diff(Ω_sf) ./ diff(Param.r); 0.0];
        dΩ[2:(Nr+1)] = Param.B_sf .* (2 * Ω_sf + EoMSetup.r .* dΩ_sfdr) .* (Ω[1] .- Ω_sf);
        dΩ_sf_net =
            4 .* π .* sum(
                ((Param.r .^ 2) .* EoMSetup.ρr .* dΩ[2:(Nr+1)]) .*
                [diff(EoMSetup.r); diff(EoMSetup.r)[end]],
            );
        # dΩ_core/dt
        dΩ[Nr+2] = 2 * Param.B_core * Ω[2] * (Ω[1] - Ω[2]);
        # dΩ_crust/dt
        dΩ[1] =
            -Param.N_ext/Param.I_crust - Param.I_core/Param.I_crust * dΩ[2] -
            dΩ_sf_net/Param.I_crust;
    end
    return ThreeComMod_inner!
end

function integral_moi_sph(
    ρ::AbstractArray,
    r::AbstractArray;
    r_range::Union{Tuple{Float64,Float64},Nothing} = nothing,
)
    if length(ρ) == length(r)
        r_low = isnothing(r_range) ? r[1] : r_range[1]
        r_up = isnothing(r_range) ? r[end] : r_range[2]
        if r_low < minimum(r)
            error("The lower bound of r_range is not in the ragne of r")
        end
        if r_up > maximum(r)
            error("The upper bound of r_range is not in the range of r")
        end
        i_low = argmin(abs.(r .- r_low))
        i_up = argmin(abs.(r .- r_up))

        i_low < 2 ? i_min = 1 : i_min = i_low-1
        i_up > length(r)-1 ? i_max = i_up : i_max = i_up + 1

        dV = ρ .* r .^ 4 .* [diff(r); diff(r)[end]];
        return 8 * π * sum(dV[i_low:i_up]) / 3
    else
        error("Input ρ and r are not in the same size.")
    end
end

function integral_moi_cyl(
    ρ::AbstractArray,
    r::AbstractArray;
    r_range::Union{Tuple{Float64,Float64},Nothing} = nothing,
)
    if length(ρ) == length(r)
        r_low = isnothing(r_range) ? r[1] : r_range[1]
        r_up = isnothing(r_range) ? r[end] : r_range[2]
        if r_low < minimum(r)
            error("The lower bound of r_range is not in the ragne of r")
        end
        if r_up > maximum(r)
            error("The upper bound of r_range is not in the range of r")
        end
        i_low = argmin(abs.(r .- r_low))
        i_up = argmin(abs.(r .- r_up))

        i_low < 2 ? i_min = 1 : i_min = i_low-1
        i_up > length(r)-1 ? i_max = i_up : i_max = i_up + 1

        dV = ρ .* r .^ 3 .* [diff(r); diff(r)[end]];
        return 2 * π * sum(dV[i_low:i_up])
    else
        error("Input ρ and r are not in the same size.")
    end
end

end
