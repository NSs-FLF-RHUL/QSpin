@testset "characteristic_lengths_TOV" begin
    using QSpin.TOV: characteristic_lengths_TOV
    using QSpin.PhysicalConstants: gravitational_constant, speed_of_light_vacuum

    # Quick-reference for characteristic length relations.
    # R = R
    # M = GR / c^2
    # Q = M c^2 / R^3
    # Rho = M / R^3

    @testset "Sanity check: 0 scale" begin
        # If 0 is provided for a scale, we should get stupid answers.
        # M = 0, the other scales are infinite.
        output = characteristic_lengths_TOV(; length = 0.0)
        @test isapprox(output.R, 0.0)
        @test isapprox(output.M, 0.0)
        @test isinf(output.Q)
        @test isinf(output.Rho)
    end

    @testset "R = 1 scale" begin
        # If the characteristic length is set to R = 1, check that the
        # returned values for the other lengths match the substituted values.
        output = characteristic_lengths_TOV(; length = 1.0)
        @test isapprox(output.R, 1.0)
        @test isapprox(output.M, gravitational_constant / speed_of_light_vacuum^2)
        @test isapprox(output.Q, gravitational_constant)
        @test isapprox(output.Rho, gravitational_constant / speed_of_light_vacuum^2)
    end

    @testset "M = 1 scale" begin
        output = characteristic_lengths_TOV(; mass = 1.0)
        @test isapprox(output.R, speed_of_light_vacuum^2 / gravitational_constant)
        @test isapprox(output.M, 1.0)
        @test isapprox(output.Q, gravitational_constant^3 / speed_of_light_vacuum^4)
        @test isapprox(output.Rho, gravitational_constant^3 / speed_of_light_vacuum^6)
    end

    @testset "Expects exactly one argument" begin
        @test_throws "Multiple, or no, length scales provided." characteristic_lengths_TOV(;
            mass = 1.0,
            pressure = 1.0,
        )
        @test_throws "Multiple, or no, length scales provided." characteristic_lengths_TOV()
    end

    @testset "Using c = G = 1 units" begin
        for R in [1.0, 2.0]
            output = characteristic_lengths_TOV(; length = R, c = 1.0, G = 1.0)
            @test isapprox(output.R, R)
            @test isapprox(output.M, R)
            @test isapprox(output.Q, 1.0/R^2)
            @test isapprox(output.Rho, 1.0/R^2)
        end
    end

end
