using QSpin
using Test

@testset "QSpin.jl" begin
    # Write your tests here.
    @test QSpin.double(2) == 4
    @test QSpin.triple(2) == 6
end
