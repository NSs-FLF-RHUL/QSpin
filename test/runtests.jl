using QSpin
using Test

@testset "ode_rk4" begin
    # This begins a set of tests for the ode_rk4 function

    @testset "ode_rk4 No Evolution" begin
        # This begins a subset of tests for the ode_rk4 function,
        # where we will use a time-independent field to confirm some simple properties.

        # A simple "EoM" that doesn't evolve the field (dψdt = 0)
        no_evolution(ψ, time) = zeros(size(ψ))
        # An arbitrary initial field, that we don't expect to evolve
        initial_field = [1.; 2.]

        # Confirm that there is no evolution according to RK4
        @test isapprox(
            QSpin.OdeSolve.ode_rk4(initial_field, 1., 0., no_evolution),
            initial_field
        )
        # When there is no evolution, the timestep δt should not matter
        @test isapprox(
            QSpin.OdeSolve.ode_rk4(initial_field, 1., 0., no_evolution),
            QSpin.OdeSolve.ode_rk4(initial_field, 10., 0., no_evolution)
        )
    end
end