module MutualFriction

using Dierckx

function MutualFriction(ParamIn::ParameterType)

    if typeof(ParamIn.B_core) == Float64
        B_core = ParamIn.B_core;
    elseif typeof(ParamIn.B_core) == Matrix{Float64}
        println("B_core Interpolation")
        spl = Spline1D(ParamIn.ρr_core, ParamIn.B_core, k = 3, bc = "extrapolate")
        B_core = spl(ParamIn.ρr);
    else
        tyopeof(ParamIn.B_core)
        println("B_core is a function of density")
        B_core = ParamIn.B_core(ParamIn.ρr);
    end

    if typeof(ParamIn.B_sf) == Float64
        B_sf = ParamIn.B_sf;
    elseif typeof(ParamIn.B_sf) == Matrix{Float64}
        println("B_sf Interpolation")
        spl = Spline1D(ParamIn.ρr_sf, ParamIn.B_sf, k = 3, bc = "extrapolate")
        B_sf = spl(ParamIn.ρr);
    else
        tyopeof(ParamIn.B_sf)
        println("B_sf is a function of density")
        B_sf = ParamIn.B_sf(ParamIn.ρr);
    end

    B_sf[B_sf .< ParamIn.ρ_dripping] .= 0;
    B_core[B_core .< ParamIn.ρ_dripping] .= 0;
    return B_crust, B_core, B_sf
end

end
