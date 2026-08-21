using QSpin
using QSpin.Parameters: ParameterType
using QSpin.OdeSolve: evolve
using QSpin.GlitchModels: ThreeCompSolid!, integral_moi_sph, integral_moi_cyl
using OrdinaryDiffEqTsit5: Tsit5

@testset "glith_riser" begin
    @testset "moment of inertia spherical integral" begin
        r = 0.0:5.0
        ρ = ones(size(r))
        @test isapprox(integral_moi_sph(ρ, r), sum(r .^ 4)*8*π/3)

    end
    @testset "moment of inertia spherical cyl" begin
        r = 0.0:5.0
        ρ = ones(size(r))
        @test isapprox(integral_moi_cyl(ρ, r), sum(r .^ 3)*2*π)

    end
    @testset "solid body" begin
        Sim_Input = (
            # Glitch Model Input
            Ω_crust = 70.34, # Initial angular velocity of the crust in rad/s
            Ω_sf = 70.34 + 6.3e-3, # Initial angular velocity of the superfluid in rad/s
            Ω_core = 70.34, # Initial angular velocity of the core in rad/s
            I_crust = 4.5e30, # Moment of inertia of the crust in kg m^2
            I_core = 0.8 * 4.5e30, # Moment of inertia of the core in kg m
            I_sf = 0.05*4.5e30,
            N_ext = 0.0,
            B_core = 0, # Mutual Friction Parameter
            B_sf = 0, # Mutual Friction Parameter
            # Glitch Model Solver Setup
            Dt = 0.1, # Time interval for recording values in the glitch model in seconds
            t_start = 0.0, # Start time for the glitch model simulation in seconds
            t_end = 120.0, # End time for the glitch model simulation in seconds
        )
        # A null propagation while the coulping is zero.
        Ω_ini = [Sim_Input.Ω_crust; Sim_Input.Ω_sf; Sim_Input.Ω_core]
        sol = evolve(
            ThreeCompSolid!,
            Ω_ini,
            Sim_Input.t_start,
            Sim_Input.t_end,
            Sim_Input;
            alg = Tsit5(),
            saveat = Sim_Input.Dt,
        )
        @test isapprox(sol[1, end], Ω_ini[1])
        @test isapprox(sol[2, end], Ω_ini[2])
        @test isapprox(sol[3, end], Ω_ini[3])

    end
end
