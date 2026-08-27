Ω_sf = Ωt[2:length(TOV_sol.r), :]
#t_plots = [0.12, 0.6, 3.0, 30]
t_plots = [0.3, 1.8, 6.0, 18.0];
idt = 1
plt1 = plot(
    TOV_sol.r[1:(end-1)]*1e-5,
    Ω_sf[:, idt],
    xflip = true,
    label = string(L"t=", 0, L"\;\textrm{s}"),
    linewidth = 3,
)
for pp = 1:length(t_plots)
    idt = Int64(t_plots[pp]/Sim_Input.Dt)+1
    plot!(
        plt1,
        TOV_sol.r[1:(end-1)]*1e-5,
        Ω_sf[:, idt],
        line = (2, :dot),
        label = string(L"t=", t_plots[pp], L"\;\textrm{s}"),
        xflip = true,
        xlabel = L"r (\mathrm{km})",
        ylabel = L"Ω_\mathrm{sf}\;(\textrm{rad/s})",
    )
end
plot!(
    plt1,
    TOV_sol.r[1:(end-1)]*1e-5,
    Ω_sf[:, end],
    linewidth = 3,
    label = string(L"t=", t[end], L"\;\textrm{s}"),
    xflip = true,
    xlims = (10.0, 10.45),
    legend = :outertopright,
    framestyle = :box,
)

plt2 = plot(t, Ωt[1, :], label = L"\Omega_\mathrm{crust}", linewidth = 2)
plot!(
    plt2,
    t,
    Ωt[end, :],
    linewidth = 2,
    label = L"\Omega_\mathrm{core}",
    xlabel = L"t\;(\mathrm{s})",
    ylabel = L"Ω\;(\textrm{rad/s})",
    title = string(L"\mathcal{B}_\mathrm{core}=", Sim_Input.B_core),
    framestyle = :box,
)

plt3 = heatmap(
    t,
    TOV_sol.r[1:(end-1)]*1e-5,
    Ω_sf,
    xlabel = L"t\quad(\textrm{s})",
    ylabel = L"r\quad(\textrm{km})",
    title = L"Ω_\mathrm{sf}(t,r)\;(\textrm{rad/s})",
    framestyle = :box,
)

l = @layout [a b; c]
plot(plt2, plt3, plt1, layout = l)
