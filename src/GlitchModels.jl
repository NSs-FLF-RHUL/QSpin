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

function ThreeCompGCA2018!(
    EoMSetup::ParameterType;
    ρ_drip::Float64 = 4e11,
    value_check::Bool = false,
)
    dr = [diff(r); diff(r)[end]]
    Nr = length(EoMSetup.r);
    I_total = 0.35 * EoMSetup.M_NS * EOMSetup.R_NS^2
    R_drip = r[argmin(abs.(EoMSetup.rho-ρ_drip))]
    I_sf = integral_moi_sph(EoMSetup.rho, EoMSetup.r; r_range = (EoMSetup.Rcci, R_drip))
    I_crust_total =
        integral_moi_sph(EoMSetup.rho, EoMSetup.r; r_range = (EoMSetup.Rcci, EoMSetup.R))
    I_core = 0.95 * (I_total - I_crust_total)
    I_crust = I_total - I_core - I_sf
    h_eff =
        0.5 * I_sf /
        integral_moi_cyl(EoMSetup.rho, EoMSetup.r; r_rnage = (EoMSetup.Rcci, EoMSetup.R))
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
        dΩ_sfdr = [diff(Ω_sf); 0.0] ./ dr;
        dΩ[2:(Nr+1)] = Param.B_sf .* (2 * Ω_sf + EoMSetup.r .* dΩ_sfdr) .* (Ω[1] .- Ω_sf);
        dΩ_sf_net =
            2 *
            h_eff *
            integral_moi_cyl(
                EoMSetup.rho .* dΩ[2:(Nr+1)],
                EoMSetup.r;
                r_range = (EoMSetp.Rcci, R_drip),
            )
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
    r_range::Union{Tuple{Float64,Float64},AbstractArray,Nothing} = nothing,
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

        dV = ρ .* r .^ 4 .* [diff(r); diff(r)[end]];
        return 8 * π * sum(dV[i_low:i_up]) / 3
    else
        error("Input ρ and r are not in the same size.")
    end
end

function integral_moi_cyl(
    ρ::AbstractArray,
    r::AbstractArray;
    r_range::Union{Tuple{Float64,Float64},AbstractArray,Nothing} = nothing,
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

        dV = ρ .* r .^ 3 .* [diff(r); diff(r)[end]];
        return 2 * π * sum(dV[i_low:i_up])
    else
        error("Input ρ and r are not in the same size.")
    end
end

end
