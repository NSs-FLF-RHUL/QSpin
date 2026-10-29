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

This is an equation of state wrapper giving simple inputs to apply different equation of states, named by

"GCA2018": The EoS used in V. Graber, A. Cumming & N. Anderson, ApJ 865, 23 (2018).
"TwoCompPoly": The EoS described by a two-component polytrope form.
"NV1973": EoS for the crust described in J. Negele & D. Vautherin, Nuclear Physics A 207, 298 (1973).
"Interp": Interpolating a given set of pressure - rho
"LinterpSkyrme": Loading a specific file based on a Skyrme model and using interpolation for the variables.

# Arguments
- `EoSName`: A string representing the names of the EoS.
- 'Parameters': A struct containing the parameters for the EoS. This is optional and can be set to `nothing` if not needed.
- 'FileInput': An optional input file for generating EoS based on specific input file.

# Returns
- 'EoS': The NamedTuple contains EoS, EoS_inv, EoS_nb, EoS_rhoN, EoS_Yp, EoS_mp, EoS_kFe and EoS_kFn are in the format of ::DataInterpolations, and that maps ρ to different quantities listed above and EoS_inv is the inverse of EoS, that ρ as a function of press. ρ_cci is a Floating number about the core-crust interface density.
"""
function EoS_Type(
    EoSName::String;
    Parameters::Union{ParameterType,Nothing} = nothing,
    FileInput = nothing,
)
    if EoSName == "LinterpSkyrme"
        EOS = EoS_LInterpSkyrme(FileInput)
        return EOS
    else
        EOS = if EoSName == "GCA2018"
            EoS_GCA2018()
        elseif EoSName == "TwoCompPoly"
            EoS_two_component_polytrope(Parameters)
        elseif EoSName == "NV1973"
            EoS_NegeleVautherin1973()
        elseif EoSName == "Interp"
            EoS_LInterp(Parameters.file_name, Parameters.EoS_indices)
        else
            error(
                "EoS Type only supports specific types: GCA2018, TwoCompPoly, NV1973, and Interp",
            )
        end
        return EOS
    end
end

end
