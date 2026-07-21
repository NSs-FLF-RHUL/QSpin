using QSpin
using Test

@testset "OdeSolve" begin
    include("odesolve/ode_rk4.jl")
end

@testset "Grids" begin
    include("grids/grid_setup.jl")

    include("grids/fft_ke.jl")
end

@testset "TOV" begin
    include("tov/tov_eq.jl")
end

@testset "Quality assurance" begin
    include("quality_assurance.jl")
end
