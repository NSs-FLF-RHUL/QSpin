using QSpin.TOV: TOV_ref_units, tov_eq!, TOV_Solve
using QSpin.TOV.EquationOfState: EoS_two_component_polytrope, EoS_GCA2018
using QSpin.PhysicalConstants:
    gravitational_constant, speed_of_light_vacuum, neutron_mass, mass_sun, hbar

@testset "tov_eq" begin
    @testset "TOV_ref_units" begin
        rho_ref = 2.8e14
        G_cgs = gravitational_constant * 1e3
        c_cgs = speed_of_light_vacuum * 1e2
        length_ref_expected = c_cgs / sqrt(G_cgs * rho_ref)
        pressure_ref_expected = rho_ref * c_cgs^2
        mass_ref_expected = rho_ref * length_ref_expected^3

        units = TOV_ref_units(; input_units = "CGS", rho_ref = rho_ref)
        @test isapprox(units.length_ref, length_ref_expected)
        @test isapprox(units.pressure_ref, pressure_ref_expected)
        @test isapprox(units.mass_ref, mass_ref_expected)
        @test isapprox(units.rho_ref, rho_ref)
        @test units.G0 == 1.0
        @test units.c0 == 1.0

        for dim_units in ("CGS_dim", "SI_dim")
            dim = TOV_ref_units(; input_units = dim_units)
            @test dim.length_ref == 1.0
            @test dim.pressure_ref == 1.0
            @test dim.mass_ref == 1.0
            @test dim.rho_ref == 1.0
        end

        @test_throws Exception TOV_ref_units(; input_units = "unsupported")
    end

    @testset "tov_eq! RHS" begin
        rho0 = 2.8e14
        EoS_inv = Returns(rho0)
        tov! = tov_eq!(EoS_inv; units = "CGS", rho_ref = rho0)
        tovUnits = TOV_ref_units(; input_units = "CGS", rho_ref = rho0)

        # Dimensionless state: P = P_phys / P_ref, m = m_phys / m_ref
        P_dim = 1.0
        m_dim = 0.1
        u = [P_dim, m_dim]
        du = similar(u)

        tov!(du, u, nothing, 0.0)
        @test du[1] == 0.0
        @test du[2] == 0.0

        r = 0.5
        tov!(du, u, nothing, r)
        rho_dim = rho0 / tovUnits.rho_ref
        @test isapprox(du[2], 4 * pi * r^2 * rho_dim)
        @test du[1] < 0.0
    end

    @testset "EoS_two_component_polytrope roundtrip" begin
        K_crust = (3 * π^2)^(2 / 3) * (hbar * 1e7)^2 / (5 * (neutron_mass * 1e3)^(8 / 3))
        ρ_b = 3e14

        EoS_Param_Stiff = (K_crust = K_crust, γ_crust = 5.0 / 3.0, γ_core = 3.0, ρ_b = ρ_b)
        EoS_Param_Soft =
            (K_crust = K_crust, γ_crust = 5.0 / 3.0, γ_core = 5.0 / 2.0, ρ_b = ρ_b)

        for params in (EoS_Param_Stiff, EoS_Param_Soft)
            EoS, EoS_inv = EoS_two_component_polytrope(params)

            # Pressure continuity at the crust-core interface
            P_below = params.K_crust * ρ_b^params.γ_crust
            K_core = params.K_crust * ρ_b^(params.γ_crust - params.γ_core)
            P_above = K_core * ρ_b^params.γ_core
            @test isapprox(P_below, P_above)
            @test isapprox(EoS(ρ_b), P_below)

            for ρ in (1e13, ρ_b, 1e15)
                @test isapprox(EoS_inv(EoS(ρ)), ρ; rtol = 1e-10)
            end
        end
    end

    @testset "TOV_Solve smoke" begin
        K_crust = (3 * π^2)^(2 / 3) * (hbar * 1e7)^2 / (5 * (neutron_mass * 1e3)^(8 / 3))
        ρ_b = 3e14
        ρ0 = 1e15
        dr = 1.0
        Dr = 1e3
        r_beg = 0.0

        EoS_Stiff, EoS_inv_Stiff = EoS_two_component_polytrope((
            K_crust = K_crust,
            γ_crust = 5.0 / 3.0,
            γ_core = 3.0,
            ρ_b = ρ_b,
        ))
        EoS_Soft, EoS_inv_Soft = EoS_two_component_polytrope((
            K_crust = K_crust,
            γ_crust = 5.0 / 3.0,
            γ_core = 5.0 / 2.0,
            ρ_b = ρ_b,
        ))

        sol_stiff = TOV_Solve(EoS_inv_Stiff, [EoS_Stiff(ρ0); 0.0], dr, Dr, r_beg)
        sol_soft = TOV_Solve(EoS_inv_Soft, [EoS_Soft(ρ0); 0.0], dr, Dr, r_beg)

        for sol in (sol_stiff, sol_soft)
            @test all(isfinite, sol.r)
            @test all(isfinite, sol.Pr)
            @test all(isfinite, sol.mr)
            @test all(isfinite, sol.ρr)
            @test isfinite(sol.R) && sol.R > 0
            @test isfinite(sol.M) && sol.M > 0

            # Pressure non-negative and decreasing until the surface
            @test all(>=(0), sol.Pr)
            @test all(<=(0), diff(sol.Pr))

            R_km = sol.R / 1e5
            M_solar = sol.M * 1e-3 / mass_sun
            @test 5.0 <= R_km <= 20.0
            @test 0.1 <= M_solar <= 3.0
        end

        @test sol_stiff.M >= sol_soft.M

        # Central structure: density near ρ0, enclosed mass near zero, density nonincreasing
        @test isapprox(sol_stiff.ρr[1], ρ0)
        @test isapprox(sol_stiff.mr[1], 0.0; atol = eps(sol_stiff.M))
        @test all(>=(0), sol_stiff.Pr)
        @test all(>=(0), sol_stiff.ρr)
        @test all(<=(0), diff(sol_stiff.ρr))
    end

    @testset "CGS ↔ SI consistency" begin
        ρ0_cgs = 1e15
        ρ_b_cgs = 3e14
        K_crust_cgs =
            (3 * π^2)^(2 / 3) * (hbar * 1e7)^2 / (5 * (neutron_mass * 1e3)^(8 / 3))
        K_crust_si = (3 * π^2)^(2 / 3) * hbar^2 / (5 * neutron_mass^(8 / 3))

        EoS_cgs, EoS_inv_cgs = EoS_two_component_polytrope((
            K_crust = K_crust_cgs,
            γ_crust = 5.0 / 3.0,
            γ_core = 3.0,
            ρ_b = ρ_b_cgs,
        ))
        EoS_si, EoS_inv_si = EoS_two_component_polytrope((
            K_crust = K_crust_si,
            γ_crust = 5.0 / 3.0,
            γ_core = 3.0,
            ρ_b = ρ_b_cgs * 1e3,
        ))

        sol_cgs = TOV_Solve(
            EoS_inv_cgs,
            [EoS_cgs(ρ0_cgs); 0.0],
            1.0,
            1e3,
            0.0;
            input_units = "CGS",
            rho_ref = 2.8e14,
        )
        sol_si = TOV_Solve(
            EoS_inv_si,
            [EoS_si(ρ0_cgs * 1e3); 0.0],
            0.01,
            10.0,
            0.0;
            input_units = "SI",
            rho_ref = 2.8e14,
        )

        # CGS: cm, g → SI: m, kg
        @test isapprox(sol_cgs.R / 100, sol_si.R)
        @test isapprox(sol_cgs.M / 1000, sol_si.M)
    end

    @testset "EoS_GCA2018 roundtrip" begin
        EoS, EoS_inv = EoS_GCA2018()
        ρ_drip = 4.0e11

        # Outer crust (analytic inverse), drip interface, and inner crust / core (root finder)
        for ρ in (1e10, ρ_drip * 0.5, ρ_drip, 1e13, 1e14, 1e15)
            @test isapprox(EoS_inv(EoS(ρ)), ρ)
        end
    end

    @testset "TOV_Solve GCA2018 smoke" begin
        EoS, EoS_inv = EoS_GCA2018()
        ρ0 = 1e15
        sol = TOV_Solve(EoS_inv, [EoS(ρ0); 0.0], 1.0, 1e3, 0.0)

        @test isfinite(sol.R) && sol.R > 0
        @test isfinite(sol.M) && sol.M > 0
        @test isapprox(sol.ρr[1], ρ0)
        @test isapprox(sol.mr[1], 0.0; atol = eps(sol.M))

        @test all(>=(0), sol.Pr)
        @test all(<=(0), diff(sol.Pr))

        R_km = sol.R / 1e5
        M_solar = sol.M * 1e-3 / mass_sun
        @test 5.0 <= R_km <= 20.0
        @test 0.1 <= M_solar <= 3.0
    end

    @testset "constant-density analytic solution" begin
        G_cgs = gravitational_constant * 1e3
        c_cgs = speed_of_light_vacuum * 1e2
        ρ0 = 1e15
        R_exact = 1e6
        M_exact = 4 * π * ρ0 * R_exact^3 / 3
        surface_factor = sqrt(1 - 2 * G_cgs * M_exact / (c_cgs^2 * R_exact))
        P0 = ρ0 * c_cgs^2 * (1 - surface_factor) / (3 * surface_factor - 1)
        EoS_inv = Returns(ρ0)

        sol = TOV_Solve(
            EoS_inv,
            [P0, 0.0],
            100.0,
            1e4,
            0.0;
            input_units = "CGS_dim",
            r_max = 1.2e6,
            reltol = 1e-10,
            abstol = 1e-12,
        )

        interior = sol.r .<= 0.9 * R_exact
        r = sol.r[interior]
        radial_factor = sqrt.(1 .- 2 * G_cgs * M_exact .* r .^ 2 ./ (c_cgs^2 * R_exact^3))
        P_exact =
            ρ0 .* c_cgs^2 .* (radial_factor .- surface_factor) ./
            (3 * surface_factor .- radial_factor)
        m_exact = 4 * π * ρ0 .* r .^ 3 ./ 3

        @test isapprox(sol.Pr[interior], P_exact)
        @test isapprox(sol.mr[interior], m_exact)
        @test isapprox(sol.R, R_exact)
        @test isapprox(sol.M, M_exact)
    end

    @testset "solver keyword forwarding" begin
        ρ0 = 0.01
        EoS_inv = Returns(ρ0)
        sol = TOV_Solve(
            EoS_inv,
            [1e-3, 0.0],
            1e-3,
            1e-2,
            0.0;
            input_units = "CGS_dim",
            r_max = 0.2,
            saveat = 0.05,
        )

        @test all(isapprox(0.05), diff(sol.r))
        @test_throws Exception TOV_Solve(
            EoS_inv,
            [1e-3, 0.0],
            1e-3,
            1e-2,
            0.0;
            alg = :unsupported,
            input_units = "CGS_dim",
            r_max = 0.2,
        )
    end

    @testset "surface and r_max behaviour" begin
        ρ0 = 0.01
        R_surface = 1.0
        M_surface = 4 * π * ρ0 * R_surface^3 / 3
        surface_factor = sqrt(1 - 2 * M_surface / R_surface)
        P0 = ρ0 * (1 - surface_factor) / (3 * surface_factor - 1)
        EoS_inv = Returns(ρ0)

        truncated = TOV_Solve(
            EoS_inv,
            [P0, 0.0],
            1e-3,
            1e-2,
            0.0;
            input_units = "CGS_dim",
            r_max = 0.5,
        )
        @test truncated.R == truncated.r[end] == 0.5
        @test truncated.M == truncated.mr[end]
        @test truncated.Pr[end] > 0

        initially_negative = TOV_Solve(
            Returns(0.0),
            [-1.0, 0.0],
            1e-3,
            1e-2,
            0.0;
            input_units = "CGS_dim",
            r_max = 0.1,
        )
        @test initially_negative.R == initially_negative.r[1] == 0.0
        @test initially_negative.M == initially_negative.mr[1] == 0.0
    end
end
