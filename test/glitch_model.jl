using QSpin
using QSpin.Parameters: ParameterType
using QSpin.OdeSolve: evolve
using QSpin.GlitchModels:
    ThreeCompSolid!, integral_moi_sph, integral_moi_cyl, moment_of_inertia
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

@testset "Moment of Inertia" begin
    @testset "i_drip and i_cci failure checks" begin
        Nr = 101
        Rmax = 1e1
        r = range(0.0, stop = Rmax, length = Nr)

        Params = (rho = ones(Nr), rhoN = ones(Nr)*0.9, r = r, R_NS = Rmax)
        MOI = moment_of_inertia(Params; ρ_drip = 0.3, ρ_cci = 0.55)
        @test MOI.i_drip == Nr
        @test MOI.i_cci == Nr-1

        Params = (rho = ones(Nr)*0.2, rhoN = ones(Nr)*0.9, r = r, R_NS = Rmax)
        MOI = moment_of_inertia(Params; ρ_drip = 0.3, ρ_cci = 0.55)
        @test MOI.i_drip == 3
        @test MOI.i_cci == 2
    end
    @testset "Moment of Inertia for homogeneous density" begin
        Nr = 201
        Rmax = 1e1
        r = range(0.0, stop = Rmax, length = Nr)
        ρr = ones(Nr)
        ρsf = 0.5
        ρcrust = 0.25
        ρr[r .> 4.95] .= ρsf
        ρr[r .> 9.45] .= ρcrust
        Params = (rho = ρr, rhoN = ρr*0.9, r = r, R_NS = Rmax)
        MOI_T = moment_of_inertia(Params; ρ_drip = 0.3, ρ_cci = 0.55, mtd = :trapezoidal)
        MOI_R = moment_of_inertia(Params; ρ_drip = 0.3, ρ_cci = 0.55, mtd = :riemann)
        R_drip = r[MOI_T.i_drip]
        R_cci = r[MOI_T.i_cci]
        @test isequal(R_drip, r[MOI_R.i_drip])
        @test isequal(R_cci, r[MOI_R.i_cci])
        # Analytical moment of inertia for a uniform density sphere
        I_anal(a, b) = 8.0 * π / 15.0 * (b^5 - a^5)
        tol_riemann_sph(a, b, Nr) = maximum(b .^ 3 * 4) * (b-a)^2 / Nr/2 * 8 * π / 3
        tol_riemann_cyl(a, b, Nr) = maximum(b .^ 2 * 3) * (b-a)^2 / Nr/2 * 2 * π

        # Riemann sum maximum error estimate

        @test isapprox(MOI_T.I_total, 1.0e45; rtol = tol_riemann_sph(0, Rmax, Nr))
        @test isapprox(
            MOI_T.I_crust_total,
            I_anal(R_cci, R_drip)*ρsf+I_anal(R_drip, Rmax)*ρcrust;
            rtol = tol_riemann_sph(R_cci, R_drip, MOI_T.i_drip-MOI_T.i_cci)*ρsf +
                   tol_riemann_sph(R_drip, Rmax, Nr-MOI_T.i_drip)*ρcrust,
        )
        @test isapprox(
            MOI_T.I_total_n,
            I_anal(0, Rmax)*0.9;
            rtol = tol_riemann_sph(0, Rmax, Nr),
        )
        @test isapprox(
            MOI_T.I_sf,
            I_anal(R_cci, R_drip)*ρsf*0.9;
            rtol = tol_riemann_sph(R_cci, R_drip, MOI_T.i_drip-MOI_T.i_cci)*ρsf,
        )
        @test isapprox(
            MOI_T.I_core,
            I_anal(0, R_cci);
            atol = tol_riemann_sph(0, R_cci, MOI_T.i_cci),
        )

        @test isapprox(MOI_R.I_total, 1.0e45; rtol = tol_riemann_sph(0, Rmax, Nr))
        @test isapprox(
            MOI_R.I_crust_total,
            I_anal(R_cci, R_drip)*ρsf+I_anal(R_drip, Rmax)*ρcrust;
            rtol = tol_riemann_sph(R_cci, R_drip, MOI_R.i_drip-MOI_R.i_cci)*ρsf +
                   tol_riemann_sph(R_drip, Rmax, Nr-MOI_R.i_drip)*ρcrust,
        )
        @test isapprox(
            MOI_R.I_total_n,
            I_anal(0, Rmax)*0.9;
            rtol = tol_riemann_sph(0, Rmax, Nr),
        )
        @test isapprox(
            MOI_R.I_sf,
            I_anal(R_cci, R_drip)*ρsf*0.9;
            rtol = tol_riemann_sph(R_cci, R_drip, MOI_R.i_drip-MOI_R.i_cci)*ρsf,
        )
        @test isapprox(
            MOI_R.I_core,
            I_anal(0, R_cci);
            atol = tol_riemann_sph(0, R_cci, MOI_R.i_cci),
        )

    end
end
