module GlitchModels
using DocStringExtensions: TYPEDSIGNATURES
using JSON: JSON
using ..Parameters: ParameterType

"""
$(TYPEDSIGNATURES)

Construct the set of equation of motion for a three-component pulsar glitch model, consisting of a crust, a superfluid component, and a core component.
Each component is represented by its angular velocity while they are regarded as rigid bodies.
The equations of motion are based on the mutual friction coefficients and the external torque applied to the crust.

# Arguments
- 'Param::ParameterType': A struct containing the parameters for the glitch model, including the mutual friction coefficients (B_sf and B_core), the external torque (N_ext), and the moments of inertia (I_crust, I_sf, I_core) for each component.

# Returns
- `ThreeCompSolid!::Function`: The function that evaluates the equations of motion for the three-component pulsar glitch model.

"""
function ThreeCompSolid!(
    dΩ::AbstractArray,
    Ω::AbstractArray,
    Param::ParameterType,
    time::Float64,
)
    Ω_sf = Ω[2]
    Ω_core = Ω[3]
    Ω_crust = Ω[1]
    dΩ[2] = Param.B_sf * (Ω_crust - Ω_sf)
    dΩ[3] = Param.B_core * (Ω_crust - Ω_core)
    dΩ[1] =
        -Param.N_ext / Param.I_crust - Param.I_sf / Param.I_crust * dΩ[2] -
        Param.I_core / Param.I_crust * dΩ[3]
end

"""
$(TYPEDSIGNATURES)

Construct the set of equation of motion for a three-component pulsar glitch model, consisting of a crust, a superfluid component, and a core component according to [Graber et al. (2018)](http://arxiv.org/abs/1804.02706).
Here the superfluid component has a radial dependence, so as the mutual friction coefficient for the sf part, with a self constent determination on the moment of inertia, while the other two regmain essentially treated as rigid bodies.
The equations of motion are based on the mutual friction coefficients and the external torque applied to the crust.

# Arguments
- 'EoMSetup::ParameterType': A struct containing the parameters for the glitch model

# Returns
- 'ThreeComMod_inner!::Function`: The function that evaluates the EoM.

"""
function ThreeCompGCA2018!(
    EoMSetup::ParameterType;
    ρ_drip::Float64 = 4e11,
    value_check::Bool = false,
)
    dr = [diff(EoMSetup.r); diff(EoMSetup.r)[end]]
    Nr = length(EoMSetup.r)
    I_total = 0.35 * EoMSetup.M_NS * EoMSetup.R_NS^2
    R_drip = EoMSetup.r[argmin(abs.(EoMSetup.rho .- ρ_drip))]
    I_sf = integral_moi_sph(EoMSetup.rho, EoMSetup.r; r_range = (EoMSetup.R_cci, R_drip))
    I_crust_total = integral_moi_sph(
        EoMSetup.rho,
        EoMSetup.r;
        r_range = (EoMSetup.R_cci, EoMSetup.R_NS),
    )
    I_core = 0.95 * (I_total - I_crust_total)
    I_crust = I_total - I_core - I_sf
    I_crust_unit = integral_moi_cyl(
        EoMSetup.rho,
        EoMSetup.r;
        r_range = (EoMSetup.R_cci, EoMSetup.R_NS),
    )
    h_eff = 0.5 * I_crust_total / I_crust_unit
    if value_check
        println(" * I_total = ", I_total)
        println(" * I_core = ", I_core)
        println(" * I_sf = ", I_sf, ", I_sf/_total = ", I_sf/I_total)
        println(
            " * I_crust_total = ",
            I_crust_total,
            ", I_crust_total/I_total = ",
            I_crust_total/I_total,
        )
        println(" * I_crust_unit = ", I_crust_unit)
        println(" * h_eff = ", h_eff)
        println(" * I_crust = ", I_crust)
        println(" * R_drip = ", R_drip)
    end

    function ThreeComMod_inner!(
        dΩ::AbstractArray,
        Ω::AbstractArray,
        Param::ParameterType,
        time::Float64,
    )
        # dΩ_sf/dt
        Ω_sf = Ω[2:(Nr+1)]
        #Bsf = Param.B_sf * Param.ρr ./ maximum(Param.ρr); # Scaling B_sf with the local density profile
        dΩ_sfdr = [diff(Ω_sf); 0.0] ./ dr
        dΩ[2:(Nr+1)] = EoMSetup.B_sf .* (2 * Ω_sf + EoMSetup.r .* dΩ_sfdr) .* (Ω[1] .- Ω_sf)
        dΩ_sf_net =
            2 *
            h_eff *
            integral_moi_cyl(
                EoMSetup.rho .* dΩ[2:(Nr+1)],
                EoMSetup.r;
                r_range = (EoMSetup.R_cci, R_drip),
            )
        # dΩ_core/dt
        dΩ[Nr+2] = 2 * EoMSetup.B_core * Ω[end] * (Ω[1] - Ω[end])
        # dΩ_crust/dt
        dΩ[1] = -EoMSetup.N_ext/I_crust - I_core/I_crust * dΩ[end] - dΩ_sf_net/I_crust
    end
    return ThreeComMod_inner!
end

function integral_range(
    r::AbstractArray,
    r_range::Union{Tuple{Float64,Float64},AbstractArray,Nothing},
)
    r_low = isnothing(r_range) ? r[1] : r_range[1]
    r_up = isnothing(r_range) ? r[end] : r_range[2]
    if r_low < minimum(r)
        error("The lower bound of r_range is not in the range of r")
    end
    if r_up > maximum(r)
        error("The upper bound of r_range is not in the range of r")
    end
    i_low = argmin(abs.(r .- r_low))
    i_up = argmin(abs.(r .- r_up))
    return i_low, i_up
end

"""
$(TYPEDSIGNATURES)

The integral for the moment of inertia in the spherical coordinate, given the density profile `ρ` and the radial profile `r`.

# Arguments
- 'ρ::AbstractArray': The density profile of the star.
- 'r::AbstractArray': The radial profile of the star.
- 'r_range::Union{Tuple{Float64,Float64},AbstractArray,Nothing}': The range of radii over which to integrate.

# Returns
- '8 * π * sum(dV[i_low:i_up]) / 3': The result of the integral for the moment of inertia in the spherical coordinate.

"""
function integral_moi_sph(
    ρ::AbstractArray,
    r::AbstractArray;
    r_range::Union{Tuple{Float64,Float64},AbstractArray,Nothing} = nothing,
)
    if length(ρ) == length(r)
        i_low, i_up = integral_range(r, r_range)

        dr = diff(r)
        dV = ρ .* r .^ 4 .* [dr; dr[end]]
        return 8 * π * sum(@view dV[i_low:i_up]) / 3
    else
        error("Input ρ and r are not in the same size.")
    end
end

"""
$(TYPEDSIGNATURES)

The integral for the moment of inertia elementin the cylindrical coordinate, given the density profile `ρ` and the radial profile `r`.

# Arguments
- 'ρ::AbstractArray': The density profile of the star.
- 'r::AbstractArray': The radial profile of the star.
- 'r_range::Union{Tuple{Float64,Float64},AbstractArray,Nothing}': The range of radii over which to integrate.

# Returns
- '2 * π * sum(dV[i_low:i_up])': The result of the integral for the moment of inertia in the cylindrical coordinate.

"""
function integral_moi_cyl(
    ρ::AbstractArray,
    r::AbstractArray;
    r_range::Union{Tuple{Float64,Float64},AbstractArray,Nothing} = nothing,
)
    if length(ρ) == length(r)
        i_low, i_up = integral_range(r, r_range)

        dr = diff(r)
        dV = ρ .* r .^ 3 .* [dr; dr[end]]
        return 2 * π * sum(@view dV[i_low:i_up])
    else
        error("Input ρ and r are not in the same size.")
    end
end

"""
$(TYPEDSIGNATURES)

Reading the input JSON file for the glitch model, and return a struct containing the parameters for the glitch model.

# Arguments
- 'file_path::String': The path to the input JSON file.

# Returns
- 'Sim_Input::ParameterType': A struct containing the parameters for the glitch model, including the mutual friction coefficients, initial conditions and the EoM solver setup.
"""
function gm_input(file_path::String)
    data = JSON.parsefile(file_path)
    @info "JSON data successfully loaded from $(file_path)"
    return Sim_Input = (
        B_core = data.EoM_Setup["B_core"],
        B_sf_type = data.EoM_Setup["B_sf_type"],
        Ω_crust = data.EoM_Initial_Conditions["Ω_crust"],
        Ω_sf = data.EoM_Initial_Conditions["Ω_sf"],
        Ω_core = data.EoM_Initial_Conditions["Ω_core"],
        N_ext = data.EoM_Setup["N_ext"],
        dt = data.EoM_Setup["dt"],
        Dt = data.EoM_Setup["Dt"],
        t_start = data.EoM_Setup["t_start"],
        t_end = data.EoM_Setup["t_end"],
        ρ0 = data.TOV["ρ0"],
        dr = data.TOV["dr"],
        Dr = data.TOV["Dr"],
        r_beg = data.TOV["r_beg"],
        R_cci = data.TOV["R_cci"],
        M_core = data.TOV["M_core"],
        tov_units = data.TOV["Units"],
    )
end
end
