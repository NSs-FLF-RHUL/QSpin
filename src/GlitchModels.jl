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
    dΩ::AbstractArray,
    Ω::AbstractArray,
    Param::ParameterType,
    time::Float64,
)
    # dΩ_sf/dt
    Ω_sf = Ω[3:end];
    #Bsf = Param.B_sf * Param.ρr ./ maximum(Param.ρr); # Scaling B_sf with the local density profile
    dΩ_sfdr = [diff(Ω_sf) ./ diff(Param.r); 0.0];
    dΩ[3:end] = Param.B_sf .* (2 * Ω_sf + Param.r .* dΩ_sfdr) .* (Ω[1] .- Ω_sf);
    dΩ_sf_net =
        4 .* π .* sum(
            ((Param.r .^ 2) .* Param.ρr .* dΩ[3:end]) .*
            [diff(Param.r); diff(Param.r)[end]],
        );
    # dΩ_core/dt
    dΩ[2] = 2 * Param.B_core * Ω[2] * (Ω[1] - Ω[2]);
    # dΩ_crust/dt
    dΩ[1] =
        -Param.N_ext/Param.I_crust - Param.I_core/Param.I_crust * dΩ[2] -
        dΩ_sf_net/Param.I_crust;
end


end
