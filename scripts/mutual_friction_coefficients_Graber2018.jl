using QSpin
using Plots, LaTeXStrings

file_path = "scripts/mutual_friction_input.json"
output = QSpin.MFriction.VNparaGraber2018(file_path)

ρc_scan = exp10.(range(13, stop = log10(2e17), length = 150)) # in kg * m^-3, which is equivalent to 1e-3 times the input in kg * fm^-3

yBew = exp10.(output.Beb_itp(log10.(ρc_scan)))
yBj = exp10.(output.Bj_itp(log10.(ρc_scan)))

plot(
    ρc_scan*1e-3,
    yBew,
    ls = :dashdot,
    lc = :blue,
    linewidth = 2,
    xaxis = (:log10, [3.5e11, 2e14]),
    yaxis = (:log10, [1e-5, 1e-1]),
    label = L"\textrm{(B)} \; \mathcal{B}_\mathrm{EB} \; \textrm{with }\; E_{p} ",
)
scatter!(output.ρs*1e-3, output.Beb, label = false, mc = :blue, lc = :blue)
plot!(
    ρc_scan*1e-3,
    yBj,
    ls = :dash,
    lc = RGB(0.94, 0.65, 0.25),
    linewidth = 2,
    xaxis = (:log10, [3.5e11, 2e14]),
    yaxis = (:log10, [1e-5, 1e-1]),
    label = L"\textrm{(C)} \; \mathcal{B}_\mathrm{J} \; \textrm{with }\; E_{p} ",
)
scatter!(
    output.ρs*1e-3,
    output.Bj,
    label = false,
    mc = RGB(0.94, 0.65, 0.25),
    lc = RGB(0.94, 0.65, 0.25),
)
xlabel!(L"\rho_s \; (\textrm{g} \; \textrm{cm}^{-3})")
ylabel!(L"\mathcal{B}")
