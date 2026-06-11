using QSpin
using Test

@testset "OdeSolve" begin
    include("odesolve/ode_rk4.jl")
end

@testset "Grids" begin
    include("grids/grid_setup.jl")

    include("grids/fft_ke.jl")
end

@testset "Quality assurance" begin
    include("quality_assurance.jl")
end

@testset "TOV Module" begin
    include("tov/dimensionless_lengths.jl")
end
