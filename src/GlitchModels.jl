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

function ThreeCompGCA2018!(EoMSetup::ParameterType)
    Nr = length(EoMSetup.r;value_check::Bool=false);
    I_total = 0.35 * M_NS * RNS^2
    dr = [diff(EoMSetup.r); diff(EoMSetup.r)[end]]
    spherical_momentum_of_innertia = EoMSetup.ρr .* EoMSetup.r.^4 .* dr
    unit_momenum_of_innertia = EvoMSetup.ρr .* EoMSetup.r.^3 .* dr
    I_crust_total = 8 * π * sum(spherical_momentum_of_innertia[i_Rcci:end])
    I_crust_unitheight = 2 * π * sum(nit_momenum_of_innertia[i_Rcci:end])
    h_crust = I_crust_unitheight / I_crust_total / 2.0

    I_core = EoMSetup.α_core * (I_total - I_crust_total)
    I_proton = EoMSetup.α_proton * (I_total - I_crust_total)
    i_R_drip
    I_sf = 4 * π * EvoMSetup.ρr .* EvoMSetup.r.^3 .* dr;
    I_sf = sum(I_sf(i_Rcci:iR_drip)) * h_crust
    if value_check
        print?
    end

    function ThreeComMod_inner!(
        dΩ::AbstractArray,
        Ω::AbstractArray,
        Param::ParameterType,
        time::Float64,
    )
        # dΩ_sf/dt
        Ω_sf = Ω[2:Nr+1];
        #Bsf = Param.B_sf * Param.ρr ./ maximum(Param.ρr); # Scaling B_sf with the local density profile
        dΩ_sfdr = [diff(Ω_sf) ./ diff(Param.r); 0.0];
        dΩ[2:Nr+1] = Param.B_sf .* (2 * Ω_sf + EoMSetup.r .* dΩ_sfdr) .* (Ω[1] .- Ω_sf);
        dΩ_sf_net =
            4 .* π .* sum(
                ((Param.r .^ 2) .* EoMSetup.ρr .* dΩ[2:Nr+1]) .*
                [diff(EoMSetup.r); diff(EoMSetup.r)[end]]
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

function integrand_sph(ρ,r;Rin::Float64=r[1],Rend::Float64=r[end])
    i_Rin
    8 * π * sum(spherical_momentum_of_innertia[i_Rcci:end])
end

function integrand_cyl(ρ,r;Rin::Float64,Rend::Float64)
    2 * π * sum(nit_momenum_of_innertia[i_Rcci:end])
end

end
