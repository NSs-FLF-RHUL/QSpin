using QSpin
using QSpin.Parameters: ParameterType
using QSpin.TOV: TOV_Solve, rho_c_from_M
using QSpin.TOV.EquationOfState: EoS_Type
using QSpin.GlitchModels: ThreeCompGeneral!, gm_input, moment_of_inertia
using Plots, LaTeXStrings
using QSpin.MFriction: mfrictionGraber2016, skyrme_effective_mass_inv
using QSpin.OdeSolve: evolve
using OrdinaryDiffEqTsit5: Tsit5
using QSpin.PhysicalConstants: mass_sun
file_name = "src/TOV/EoS/WNewton/eos_SkXi450_31.0_60.00_-100.00_glitch.dat"

ms_all = [0.6, 0.8, 1.0]
Δm_all = [-0.3, 0.17, 0.64]

Sim_Input = (
    ρ0 = 1e15, # Initial central density in g/cm^3, above 2e16 seems to be unstable
    dr = 0.0001*1e6, # Radial step in cm
    Dr = 0.0005*1e6, # Radial interval for recording values in cm
    r_beg = 0.e5,
    units = "CGS", # optional
    B_sf = 1e-3, # Mutual Friction Parameter for the SF component
    N_ext = 0,
    Ω_vela = 2*π*11.195, # Initial angular velocity of the crust in rad/s
    dnu_vela = 16.0227,
    ρ_drip = 4e11, # cgs units
    t_end = 1e4, # Total simulation time in seconds
    Dt = 1e-1, # Time step for saving data in seconds
);

DataLenY = Int(length(ms_all) * length(Δm_all)) + 1
DataLenX = Int(Sim_Input.t_end / Sim_Input.Dt) + 1 + 4

Data = zeros(DataLenY, DataLenX)

EOS = EoS_Type("LinterpSkyrme", FileInput = file_name)
ρ0 = rho_c_from_M(EOS.EoS, EOS.EoS_inv, 1.4*1e3*mass_sun)
u0 = [EOS.EoS(ρ0); 0.0];

# Callback setup to terminate the integration when the pressure drops below zero
@time TOV_sol = TOV_Solve(
    EOS.EoS_inv,
    u0,
    Sim_Input.dr,
    Sim_Input.Dr,
    Sim_Input.r_beg;
    reltol = 1e-13,
    abstol = 1e-13,
    input_units = Sim_Input.units, # optional
    rho_ref = 2.8e14, # optional - nuclear saturation density in g/cm^3
);

mp = EOS.EoS_mp(TOV_sol.ρr)
mp[mp .< 0] .= 0.0;
kFn = EOS.EoS_kFn(TOV_sol.ρr)
kFn[kFn .< 0.0] .= 0.0;
kFe = EOS.EoS_kFe(TOV_sol.ρr)
kFe[kFe .< 0.0] .= 0.0;
Yp = EOS.EoS_Yp(TOV_sol.ρr)
Yp[Yp .< 0.0] .= 0.0;
nb = EOS.EoS_nb(TOV_sol.ρr)
nb[nb .< 0.0] .= 0.0;

MOI = moment_of_inertia((
    rho = TOV_sol.ρr,
    rhoN = EOS.EoS_rhoN(TOV_sol.ρr),
    r = TOV_sol.r,
    R_NS = TOV_sol.R,
))

i_drip = MOI.i_drip
i_cci = MOI.i_cci

Ω_c0 = Sim_Input.Ω_vela
Ω_sf0 = Sim_Input.Ω_vela + 2 * π * 1e-6 * Sim_Input.dnu_vela * MOI.I_total/MOI.I_sf
Ω_core0 = Sim_Input.Ω_vela

# Setting initial conditions for the glitch model
Ω_ini = [
    Ω_c0;
    Ω_sf0 * ones(size(TOV_sol.r[(i_cci+1):i_drip]));
    Ω_core0 * ones(size(TOV_sol.r[1:i_cci]))
]

cc = let cc = 1
    for aa = 1:length(ms_all)
        for bb = 1:length(Δm_all)
            a, b = skyrme_effective_mass_inv(0.16, ms_all[aa], Δm_all[bb])
            @info a, b, ms_all[aa], Δm_all[bb]
            Data[1+cc, 1:4] = [a, b, ms_all[aa], Δm_all[bb]]
            Params = (
                kFn = kFn,
                kFe = kFe,
                nb = nb,
                Yp = Yp,
                ρ = TOV_sol.ρr,
                B_sf = 4e-4,
                a = a,
                b = b,
            )
            Bs = mfrictionGraber2016("p", Params)

            # Wrap the parameters in the glith model
            EoMSetup = (
                rho = TOV_sol.ρr,
                rhoN = EOS.EoS_rhoN(TOV_sol.ρr),
                r = TOV_sol.r,
                M_NS = TOV_sol.M,
                R_NS = TOV_sol.R,
                ρ_cci = TOV_sol.ρr[i_cci], # 10 km in cm
                B_sf = Sim_Input.B_sf, # Mutual Friction Parameter
                B_core = Bs.B_core[1:i_cci], # Mutual Friction Parameter
                N_ext = Sim_Input.N_ext,
                ρ_drip = Sim_Input.ρ_drip, # Drip density in g/cm^3
            )

            EoM! = ThreeCompGeneral!(EoMSetup, MOI)
            @time sol = evolve(
                EoM!,
                Ω_ini,
                0.0,
                Sim_Input.t_end,
                Sim_Input;
                alg = Tsit5(),
                saveat = Sim_Input.Dt,
                reltol = 1e-14,
                abstol = 1e-14,
            )

            # Data Reading and Plotting
            Data[cc+1, 5:end] = Array(sol)[1, :];
            Data[1, 5:end] = sol.t;
            cc += 1

        end
    end
end

t = Data[1, 5:end]
Ωt = Data[2:end, 5:end]
pp = let pp = 1
    for aa = 1:length(ms_all)
        for bb = 1:length(Δm_all)
            @info pp, aa, bb
            if pp == 1
                plot(
                    t[2:end],
                    (Ωt[pp, 2:end] .- Ωt[pp, 1])*1e6/2/π,
                    xaxis = :log10,
                    yaxis = :log10,
                    xlabel = L"t\;(\textrm{s})",
                    ylabel = L"\Delta\nu\;(\mu\textrm{Hz})",
                    xticks = [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4],
                    yticks = [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4],
                    framestyle = :box,
                    xminorticks = 10,
                    yminorticks = 10,
                    lw = 2,
                    ms = 2,
                    label = "m_s= $(ms_all[aa]), δm/δ= $(Δm_all[bb])",
                    xlims = [0.5, 1e4],
                    ylims = [0.5, 1e3],
                    legend = :bottomleft,
                )
            else
                plot!(
                    t[2:end],
                    (Ωt[pp, 2:end] .- Ωt[pp, 1])*1e6/2/π,
                    lw = 2,
                    ms = 2,
                    label = "m_s= $(ms_all[aa]), δm/δ= $(Δm_all[bb])",
                )
            end
            pp += 1
        end
    end

end
gui()
