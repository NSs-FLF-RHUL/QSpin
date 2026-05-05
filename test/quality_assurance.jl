using QSpin: QSpin
using Aqua: Aqua
using Test: @testset, @test

@testset "Aqua" begin
    Aqua.test_all(QSpin)
end
