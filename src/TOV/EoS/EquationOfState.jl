"""
Submodule containing specific Equations of State that may be used in the TOV equations.
"""
module EquationOfState

include("GraberCummingAnderson2018.jl")
include("NegeleVautherin1973.jl")
include("TwoComponentPolytrope.jl")
include("EoS_LInterp.jl")
include("EoS_LInterpSkyrme.jl")

"""
$(TYPEDSIGNATURES)

Equation of motion control

# Arguments
- `EoSName`: A string representing the names of the EoS.
- 'Parameters': A struct containing the parameters for the EoS. This is optional and can be set to `nothing` if not needed.

# Returns
- 'EoS': The pressure-denisty relation for the specified EoS.
- 'EoS_inv': The density presure relation of the specific EoS.

"""
function EoS_Type(
    EoSName::String;
    Parameters::Union{ParameterType,Nothing} = nothing,
    FileInput = nothing,
)
    if EoSName == "LinterpSkyrme"
        EoS, EoS_inv, EoS_ρ2nb, EoS_ρ2Yp, EoS_ρ2mp_ast, EoS_ρ2kFe, EoS_ρ2kFn =
            EoS_LInterpSkyrme(FileInput)
        return EoS, EoS_inv, EoS_ρ2nb, EoS_ρ2Yp, EoS_ρ2mp_ast, EoS_ρ2kFe, EoS_ρ2kFn
    else
        EoS, EoS_inv = if EoSName == "GCA2018"
            EoS_GCA2018()
        elseif EoSName == "TwoCompPoly"
            EoS_two_component_polytrope(Parameters)
        elseif EoSName == "NV1973"
            EoS_NegeleVautherin1973()
        elseif EoSName == "Interp"
            EoS_LInterp(Parameters.file_name, Parameters.EoS_indices);
        else
            error(
                "EoS Type only supports specific types: GCA2018, TwoCompPoly, NV1973, and Interp",
            )
        end
        return EoS, EoS_inv
    end
end

end
