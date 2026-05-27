module MFriction
using JSON
using ..PhysicalConstants: hbar, neutron_mass, electron_volt, gravitational_constant


function ReadJSON(file_path)
    return JSON.parsefile(file_path)
end

function VNparaGraber2018(file_path)

    MeV = electron_volt * 1e6 # convert to kg * fm^2 / s^2
    δ = 1e-2 # dimensionless coefficient for the pinning energy reduction due to vortex tension
    κ = hbar / neutron_mass * π  # SI units

    data = ReadJSON(file_path)
    println("Successfully loaded JSON data!\n")
    nb = zeros(length(data)-1)
    Z = zeros(length(data)-1)
    N = zeros(length(data)-1)
    x = zeros(length(data)-1)
    ns = zeros(length(data)-1)
    a = zeros(length(data)-1)
    Rn = zeros(length(data)-1)
    Es = zeros(length(data)-1)
    E1 = zeros(length(data)-1)
    ΔE = zeros(length(data)-1)
    ξ = zeros(length(data)-1)
    Ep = zeros(length(data)-1)
    #    A = zeros(length(data)-1)
    #    RWS = zeros(length(data)-1)
    for dd = 2:length(data)
        nb[dd-1] = data[dd].nb
        Z[dd-1] = data[dd].Z
        N[dd-1] = data[dd].N
        x[dd-1] = data[dd].x
        ns[dd-1] = data[dd].ns
        a[dd-1] = data[dd].a
        Rn[dd-1] = data[dd].RN
        Es[dd-1] = data[dd].Es
        E1[dd-1] = data[dd].E1
        ΔE[dd-1] = data[dd].DE
        ξ[dd-1] = data[dd].xi
        Ep[dd-1] = data[dd].Ep
        #A[dd-1]   = data[dd].Z * (1 + 1/data[dd].x)
        #Rws[dd-1] = (3*(data[dd].N+data[dd].Z) / (4 * π * data[dd].nb*1e-4))^(1/3)
    end


    A = Z .* (1 .+ 1 ./ x)
    Rws = (3*(N .+ Z) ./ (4 * π * nb * 1e-4)) .^ (1/3)
    n1 = 3/4/π ./ Rws .^ 3*1e6
    ρs = ns * 1e-4 * neutron_mass * 1e45# in kg * fm^-3
    Reb =
        2.8 * sqrt.(0.5 * neutron_mass/hbar) * sqrt.(abs.(Ep * MeV) * δ ./ ρs / κ) .* Rn ./
        a .^ (1.5) * 10^(7.5)
    Rj =
        0.5 / sqrt(π) *
        (0.5 * neutron_mass/hbar)^(1/2) *
        (abs.(Ep * MeV) * δ ./ ρs / κ) .^ (1/2) .* a .^ 0.5 ./ ξ * 10^(7.5)


    output = (
        nb = nb,
        Z = Z,
        N = N,
        x = x,
        ns = ns,
        n1 = n1,
        a = a,
        Rn = Rn,
        Es = Es,
        E1 = E1,
        ΔE = ΔE,
        ξ = ξ,
        Ep = Ep,
        A = A,
        Rws = Rws,
        ρs = ρs,
        Reb = Reb,
        Rj = Rj,
    )



    return output
end

end
