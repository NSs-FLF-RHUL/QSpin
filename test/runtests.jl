using MyJuliaPackage
using Test

@testset "MyJuliaPackage.jl" begin
    # Write your tests here.
    @test MyJuliaPackage.double(2) == 4
    @test MyJuliaPackage.triple(2) == 6
end
