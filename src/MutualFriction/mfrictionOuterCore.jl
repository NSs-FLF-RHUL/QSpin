using SciMLBase: ContinuousCallback, terminate!, ODEProblem
using CommonSolve: CommonSolve
using OrdinaryDiffEqTsit5: Tsit5
using SpecialFunctions: besselj1

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

function mfrictionGraber2016(type::String, Params::ParameterType)
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
    kFn = Params.kFn
    kFb = Params.kFb
    kFe = Params.kFe
    B_sf = Params.B_sf

    Δn = @. Δ0 * (kFn - g0) ^ 2 / ((kFn - g0) ^ 2 + g1) * (kFn - g2) ^ 2 /
       ((kFn-g2) .^ 2 + g3)
    Δn[kFn .< g0] .= 1e-9
    Δn[kFn .> g2] .= 1e-9

    mn_ast, mp_ast = skyrme_effective_mass(Params.nb, Params.Yp, Params.a, Params.b)
    println(mn_ast)

    println(mp_ast)
    β1 = @. 4.1 *
       sqrt(
           (1.0 / mn_ast) *
           (mn_ast - 1.0 + mp_ast)^(-1) *
           Params.nb *
           1e-14 *
           (Params.Yp/0.05),
       ) *
       (kFn / 2.0) *
       (0.05 / Δn)
    β2 = @. 8e2 * (1.0 / mn_ast) * (kFe / 0.75) * (kFb / 2.0)

    B_core =
        @. 3 * π / 2 * Params.Yp / (1-Params.Yp) * (1/mn_ast)^2 * (1 - mp_ast)^2 * β1^4 /
           β2^3 * B_integral(Bcore_integrand, (β1, β2))
    return B_sf, B_core
end

function skyrme_effective_mass(
    nb::Union{Float64,AbstractArray},
    Yp::Union{Float64,AbstractArray},
    a::Float64,
    b::Float64,
)
    δ = 1 .- 2 * Yp

    mp_ast = @. 1 / (1 + a * nb + b * nb * δ)
    mn_ast = @. 1 / (1 + a * nb - b * nb * δ)

    return mn_ast, mp_ast
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

function jinc(x::Union{Float64,AbstractArray,AbstractVector})
    @. iszero(x) ? 0.5 : besselj1(x) / (x)
end

function Bcore_integrand(u, params, t)
    β1, β2 = params
    return (β2^2 + 0.5 * t .^ 2) ./ (β1 + t .^ 2)^2 * (jinc(t)) .^ 2
end

function B_integral(integrand::Function, params)
    u0 = 0.0
    problem = ODEProblem(integrand, u0, (0.0, params[2]), params)
    sol = CommonSolve.solve(
        problem,
        alg = Tsit5(),
        reltol = 1e-6,
        abstol = 1e-6,
        saveat = params[2],
    )
    return sol.u[end][1]
end
