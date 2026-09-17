"""
$(TYPEDSIGNATURES)

This function reads in JSON data for a given file path and calculates the mutual friction parameters for a neutron star crust based on the Graber et al. 2018 model.
According to Graber et al. 2018, the mutual friction coefficients are calculated based on the superfluid density and other physical parameters.
The function returns a tuple containing the input parameters and the calculated mutual friction parameters, including the qubic spline interpolations for the mutual friction coefficients as functions of the superfluid density (in kg * m^-3, while the coverted input is in kg fm^-3) in log-log space.

Here things are converted in the SI units so as MutualFrictionCoefficients

# Arguments
- `file_path`: A string representing the path to the JSON file containing the input parameters for the mutual friction calculations. The JSON file should contain an array of objects, each representing a different region of the neutron star crust with specific parameters such as baryon number density (nb), proton number (Z), neutron number (N), proton fraction (x), superfluid density (ns), lattice spacing (a), nuclear radius (RN), and pinning energy parameters (Es, E1, DE, xi, Ep).

# Returns
- 'output': A tuple containing the input parameters (in their original units from the input JSON file) and the calculated mutual friction parameters in array forms. The qubic spline interpolations for the mutual friction coefficients, B_EW and B_J, as functions of the superfluid density (in kg * m^-3, while the coverted input is in kg fm^-3) are included.
"""

function mfrictionGraber2016(type::String)
    if type == "s"
        Δ0 = 68.0
        g0 = 0.1
        g1 = 4.0
        g2 = 1.7
        g3 = 4.0
    elseif type == "p"
        Δ0 = 0.068
        g0 = 1.28
        g1 = 0.1
        g2 = 2.37
        g3 = 0.02
    else
        throw(
            ArgumentError(
                "type must be 's' or 'p' for single and triplet pairing SFs respectively",
            ),
        )
    end
    kFn = kFe = B_sf = 4e-4

    Δn = @. Δ0 * (kFn - g0) ^ 2 / ((kFn - g0) ^ 2 + g1) * (kFn - g2) ^ 2 /
       ((kFn-g2) .^ 2 + g3)
    Δn[kFn .< g0] .= 1e-9
    Δn[kFn .> g2] .= 1e-9
    β1 = @. 4.1 *
       sqrt((1/mnast) * (mnstar - 1 + mpstar)^(-1) * ρ * 1e-14 * (yp/0.05)) *
       (kFn / 2.0) *
       (0.05 / Δn)
    β2 = @. 8e2 * (1/mnast) * (kFe / 0.75) * (kFb / 2)
    B_core

end

function gap_n(kFn::Union{Float64,AbstractArray}, type::String)
    if type == "s"
        Δ0 = 68.0
        g0 = 0.1
        g1 = 4.0
        g2 = 1.7
        g3 = 4.0
    elseif type == "p"
        Δ0 = 0.068
        g0 = 1.28
        g1 = 0.1
        g2 = 2.37
        g3 = 0.02
    else
        throw(
            ArgumentError(
                "type must be 's' or 'p' for single and triplet pairing SFs respectively",
            ),
        )
    end

    Δn = @. Δ0 * (kFn - g0) ^ 2 / ((kFn - g0) ^ 2 + g1) * (kFn - g2) ^ 2 /
       ((kFn-g2) .^ 2 + g3)
    Δn[kFn .< g0] .= 1e-9
    Δn[kFn .> g2] .= 1e-9
    return Δn
end
