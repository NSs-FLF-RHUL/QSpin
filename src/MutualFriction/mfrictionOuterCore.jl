using SciMLBase: ODEProblem
using CommonSolve: CommonSolve
using OrdinaryDiffEqTsit5: Tsit5
using SpecialFunctions: besselj1
using DataInterpolations: ExtrapolationType, QuadraticSpline

"""
$(TYPEDSIGNATURES)

This fucntion computes the mutual friction coefficient for the outer core region of a neutron star, B_core, using the model descibed in Graber's thesis in 2016 and reading that for the inner crust region, B_sf, that can be an arbitary input in the Params input arguement. The types of neutron superfluid have different gap parameters, controlled by "type" of the input argruments.

# Arguments
- 'type': A string for the type of neutron superfluid. It can be "s" or "p" for singlet and triplet pairings.
- `Params`: A parameter set contains B_sf, Yp, kFe, kFn and the Skyrme model parameters, a and b.

# Returns
- 'B_sf' and 'B_core' in NamedTuple format.
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
    kFe = Params.kFe
    B_sf = Params.B_sf

    Δn = @. Δ0 * (kFn - g0) ^ 2 / ((kFn - g0) ^ 2 + g1) * (kFn - g2) ^ 2 /
       ((kFn-g2) .^ 2 + g3)
    Δn[kFn .< g0] .= 1e-9
    Δn[kFn .> g2] .= 1e-9

    mn_ast, mp_ast = skyrme_effective_mass(Params.nb*1e-39, Params.Yp, Params.a, Params.b)

    β1 = @. 4.1 *
       sqrt(
           (1.0 / mn_ast) *
           abs(mn_ast - 1.0 + mp_ast)^(-1) *
           Params.ρ *
           1e-14 *
           (Params.Yp/0.05),
       ) *
       (kFn / 2.0) *
       (0.05 / Δn)
    β2 = @. 8e2 * (1.0 / mn_ast) * (kFe / 0.75) * (kFn / 2.0) * (0.05 / Δn)


    integral = ifelse.(
        (kFn .> g0) .& (kFn .< g2),
        B_integral.(Bcore_integrand, 0.0, β2, β1, β2),
        0.0,
    )

    B_core =
        @. 3 * π / 2 * Params.Yp / (1-Params.Yp) * (1/mn_ast)^2 * (1 - mp_ast)^2 * β1^4 /
           β2^3 * integral
    return (; B_sf, B_core)
end


"""
$(TYPEDSIGNATURES)

Computing the effective proton and neutron masses in the Skyrme model.

# Arguments
- 'nb::Union{Float64,AbstractArray}': bayronic number density.
- 'Yp::Union{Float64,AbstractArray};: Proton fraction.
- 'a::Float64': Skyrme model parameter.
- 'b::Float64': Skyrme model parameter.
# Returns
- 'mn_ast': effective neutron mass.
- 'mp_ast': effective proton mass.
"""
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
"""
$(TYPEDSIGNATURES)

Computing a and b in the Skyrme model parameters for given set of effective isoscalar and effective mass splitting ratio and baryonic number density.

# Arguments
- 'nb::Union{Float64,AbstractArray}': bayronic number density.
- 'ms::Union{Float64,AbstractArray}': isoscalar mass ratio.
- 'Δm::Union{Float64,AbstractArray': effective mass splitting ratio.

# Returns
- 'a::Float64': Skyrme model parameter.
- 'b::Float64': Skyrme model parameter.
"""
function skyrme_effective_mass_inv(
    nb::Union{Float64,AbstractArray},
    ms::Union{Float64,AbstractArray},
    Δm::Union{Float64,AbstractArray},
)
    a = @. (1 / ms - 1) / nb
    b = @. 0.5 / nb * Δm / ms^2
    return a, b
end

function jinc(x)
    iszero(x) ? 0.5 : besselj1(x) / (x)
end

"""
$(TYPEDSIGNATURES)

The integrand function in the computation of mutual friction coefficient in the core componont.

# Arguments
- 'u': general sytax for using ODE solver.
- 'params::Union{Float64,AbstractArray}': β1 and β2 set.
- 't': Integrated variable.

# Returns
- 'a::Float64': Skyrme model parameter.
- 'b::Float64': Skyrme model parameter.
"""
function Bcore_integrand(u, params, t)
    β1, β2 = params
    return @. (β2^2 + 0.5 * t ^ 2) / (β1^2 + t ^ 2) ^ 2 * jinc(t) ^ 2
end

"""
$(TYPEDSIGNATURES)

Computing the integral of mutual friction coefficient in the core componont.

# Arguments
- 'integrand::Function': general sytax for using ODE solver.
- 'lower_bound::Union{Float64,AbstractArray}': lower bound of the integral.
- 'upper_bound::Union{Float64,AbstractArray}': upper bound of the integral.
- 'β1::Union{Float64,AbstractArray}': β1 coefficient.
- 'β2::Union{Float64,AbstractArray}': β2 coefficient.

# Returns
- 'a::Float64': Skyrme model parameter.
- 'b::Float64': Skyrme model parameter.
"""
function B_integral(
    integrand::Function,
    lower_bound::Union{Float64,AbstractArray},
    upper_bound::Union{Float64,AbstractArray},
    β1::Union{Float64,AbstractArray},
    β2::Union{Float64,AbstractArray},
)
    u0 = 0.0
    if β2 == 0.0
        return 0.0
    end
    problem = ODEProblem(integrand, u0, (lower_bound, upper_bound), (β1, β2))
    sol =
        CommonSolve.solve(problem, alg = Tsit5(), reltol = 1e-6, abstol = 1e-6, saveat = β2)
    return sol.u[end][1]
end
