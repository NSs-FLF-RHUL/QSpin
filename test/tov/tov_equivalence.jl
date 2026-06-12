
@testset "TOV Equivalence" begin
    using QSpin.PhysicalConstants: hbar, mass_sun, neutron_mass
    using QSpin.TOV: solve_TOV, TOV_Solve, characteristic_lengths_TOV
    using QSpin.TOV.EquationOfState: EoS_two_component_polytrope

    import OrdinaryDiffEq as DE

    # Polytropic EoS parameters
    EoS_Param = (
        K_crust = (3*π^2)^(2/3) * hbar^2 / (5*neutron_mass^(8/3)),
        γ_crust = 5.0/3.0,
        γ_core = 3.0,
        ρ_b = 3e14*1e3,
    )
    # Simulation Input Parameters
    Sim_Input = (
        rho0 = 1e15*1e3, # Initial central density in kg/m^3 (1e15 g/cm^3)
        dr = 0.005*1e3, # Radial step in meters
        Dr = 0.01*1e3, # Radial interval for recording values in meters
        r_end = 20e3, # Maximum radius to solve up to in meters
    );
    # Function Setup for inverse EoS and TOV equation for the solver
    eos_forward, eos_backward = EoS_two_component_polytrope(EoS_Param);

    u0 = [eos_forward(Sim_Input.rho0); 0.0];

    # Solve the TOV equation using `QSpin` for the stiff and soft EoSs.
    TOV_sol_dimfull = TOV_Solve(
        u0,
        Sim_Input.dr,
        Sim_Input.Dr,
        Sim_Input.r_end,
        eos_backward;
        alg = DE.Tsit5(),
        reltol = 1e-8,
        abstol = [1.0, 1e15],
    )
    char_lengths = characteristic_lengths_TOV(; length = 1.0, c = 1, G = 1.0)
    # BUT THE EOS is not in C=1, G=1 units (the constants given are in whatever random units we choose!)
    # IE K_crust etc are not in C=1, G=1 units and the EoS translator does not account for this. Need to
    # be more careful in how the EoS are defined, probably best to ditch the non-dimensionaliser and just pass
    # the char lengths into it?
    TOV_sol_dimless = solve_TOV(
        u0,
        Sim_Input.dr,
        Sim_Input.Dr,
        Sim_Input.r_end,
        char_lengths,
        eos_backward;
        alg = DE.Tsit5(),
    )
    println(TOV_sol_dimfull.r)
    println(TOV_sol_dimless.r)


    # @test isclose(TOV_sol_dimfull.t, TOV_sol_dimless.t)
    # @test isclose(TOV_sol_dimfull.u, TOV_sol_dimless.u)
end
