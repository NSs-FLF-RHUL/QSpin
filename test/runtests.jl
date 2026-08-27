using QSpin
using Test

@testset "OdeSolve" begin
    include("odesolve/ode_rk4.jl")
end

@testset "Grids" begin
    include("grids/grid_setup.jl")

    include("grids/fft_ke.jl")
end

include("hamiltonian.jl")

@testset "TOV" begin
    include("tov/tov_eq.jl")
    include("tov/eos.jl")
    include("tov/glitch_model.jl")
end

@testset "MFriction" begin
    include("mfriction/mfriction.jl")
end
@testset "gltich model" begin
    include("tov/glitch_model.jl")
end

@testset "Quality assurance" begin
    include("quality_assurance.jl")
end
