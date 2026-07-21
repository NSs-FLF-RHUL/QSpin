@testset "ode_rk4" begin
    # This begins a set of tests for the ode_rk4 function

    @testset "ode_rk4 No Evolution" begin
        # This begins a subset of tests for the ode_rk4 function,
        # where we will use a time-independent field to confirm some simple properties.

        # A simple "EoM" that doesn't evolve the field (dψdt = 0)
        no_evolution(ψ, time) = zeros(size(ψ))
        # An arbitrary initial field, that we don't expect to evolve
        initial_field = [1.0; 2.0]

        # Confirm that there is no evolution according to RK4
        @test isapprox(
            QSpin.OdeSolve.ode_rk4(initial_field, 1.0, 0.0, no_evolution),
            initial_field,
        )
        # When there is no evolution, the timestep δt should not matter
        @test isapprox(
            QSpin.OdeSolve.ode_rk4(initial_field, 1.0, 0.0, no_evolution),
            QSpin.OdeSolve.ode_rk4(initial_field, 10.0, 0.0, no_evolution),
        )
    end

    @testset "evolve default algorithm" begin
        decay!(du, u, _, _) = (du .= -u)
        solution = QSpin.OdeSolve.evolve(decay!, [1.0], 0.0, 1.0)

        @test solution.t[end] == 1.0
        @test isapprox(solution.u[end][1], exp(-1); rtol = 2e-5)
    end

    @testset "evolve_rk4 save scheduling" begin
        constant_rate(u, _) = ones(size(u))
        fields, times = QSpin.OdeSolve.evolve_rk4([0.0], 0.1, 0.3, 0.6, constant_rate)

        @test times ≈ [0.0, 0.3, 0.6]
        @test vec(fields) ≈ times
        @test_throws ArgumentError QSpin.OdeSolve.evolve_rk4(
            [0.0],
            0.2,
            0.1,
            0.4,
            constant_rate,
        )
        @test_throws ArgumentError QSpin.OdeSolve.evolve_rk4(
            [0.0],
            0.1,
            0.25,
            0.5,
            constant_rate,
        )
    end

    @testset "evolve_rk4 NaN truncation" begin
        nan_rate(u, _) = fill(NaN, size(u))
        fields, times = QSpin.OdeSolve.evolve_rk4([1.0], 0.1, 0.2, 1.0, nan_rate)

        @test size(fields) == (1, 1)
        @test fields[:, 1] == [1.0]
        @test times == [0.0]
    end
end
