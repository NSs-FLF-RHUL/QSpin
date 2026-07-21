@testset "Hamiltonian" begin
    kinetic(ψ, _, _) = 2 .* ψ
    potential(ψ, _, _) = 3 .* ψ
    ψ = ComplexF64[1, 2]

    dψ = similar(ψ)
    QSpin.Hamiltonian.hamiltonian!(kinetic, potential)(dψ, ψ, (), 0.0)
    @test dψ == -5im .* ψ

    imaginary_time = QSpin.Hamiltonian.hamiltonian!(kinetic, potential, -1)
    imaginary_time(dψ, ψ, (), 0.0)
    @test dψ == -5 .* ψ
end
