@testset "Non-dimensionalise EoS" begin
    using QSpin.TOV: characteristic_lengths_TOV
    using QSpin.TOV.EquationOfState: nondimensional_EoS

    function simple_builder(gradient::Number = 1.0; offset::Number = 0.0)
        # Forward "EoS" is y = gradient * x - offset.
        # Arguments deliberately one positional and one keyword for testing purposes.
        return x -> gradient * x - offset, y -> (y + offset) / gradient
    end

    r_eq_1_lengths = characteristic_lengths_TOV(; length = 1.0)
    n_pts = 1000
    dimensionless_range = LinRange(-5.0, 5.0, n_pts)

    @testset "Scaling matches" begin
        forward, backward = simple_builder()
        nd_forward, nd_backward = nondimensional_EoS(r_eq_1_lengths, simple_builder)

        # Check that expected scaling has been applied
        @test isapprox(
            nd_forward.(dimensionless_range),
            forward.(dimensionless_range * r_eq_1_lengths.Rho) / r_eq_1_lengths.Q,
        )
        @test isapprox(
            nd_backward.(dimensionless_range),
            backward.(dimensionless_range * r_eq_1_lengths.Q) / r_eq_1_lengths.Rho,
        )
    end

    @testset "Behaves nicely with optional args" begin
        offset = r_eq_1_lengths.Rho
        gradient = 2.0

        forward, backward = simple_builder(gradient; offset = offset)
        nd_forward, nd_backward =
            nondimensional_EoS(r_eq_1_lengths, simple_builder, gradient; offset = offset)

        @test isapprox(
            nd_forward.(dimensionless_range),
            forward.(dimensionless_range * r_eq_1_lengths.Rho) / r_eq_1_lengths.Q,
        )
        @test isapprox(
            nd_backward.(dimensionless_range),
            backward.(dimensionless_range * r_eq_1_lengths.Q) / r_eq_1_lengths.Rho,
        )
    end
end
