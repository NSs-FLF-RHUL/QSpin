using QSpin.MFriction: VNparaGraber2018, MutualFrictionCoefficients

@testset "MFriction" begin
    input_file = joinpath(@__DIR__, "..", "..", "scripts", "mutual_friction_input.json")

    @testset "VNparaGraber2018" begin
        output = @test_logs (:info, r"JSON data successfully loaded") VNparaGraber2018(
            input_file,
        )

        @test length(output.nb) == 5
        @test output.nb == [8.8, 57.7, 204.0, 475.0, 789.0]
        @test output.ns == [4.8, 47.0, 184.0, 436.0, 737.0]
        @test output.ρ == [1.5, 9.6, 33.9, 78.9, 131.0]*1e15
        @test output.A ≈ output.Z .* (1 .+ 1 ./ output.x)

        for values in (
            output.n1,
            output.Rws,
            output.ρs,
            output.R.Reb,
            output.R.Rj,
            output.B.Beb,
            output.B.Bj,
        )
            @test all(isfinite, values)
            @test all(>(0), values)
        end
        @test all(<=(0.5), output.B.Beb)
        @test all(<=(0.5), output.B.Bj)
        @test all(>(0), diff(output.ρs))

        # The splines are built in log-log space and must reproduce their knots.
        @test output.B_itp.Beb_itp.(output.ρ) ≈ output.B.Beb
        @test output.B_itp.Bj_itp.(output.ρ) ≈ output.B.Bj

        # Left extrapolation is constant by construction.
        below_range = log10(first(output.ρ)) - 1
        @test output.B_itp.Beb_itp(below_range) ≈ first(output.B.Beb)
        @test output.B_itp.Bj_itp(below_range) ≈ first(output.B.Bj)

        @test_throws SystemError VNparaGraber2018(input_file * ".missing")
    end

    @testset "MutualFrictionCoefficients regions" begin
        ρ_drip = 4e14
        Rcci = 1e4
        Param = (
            ρ = [1e13, 1e15, 1e13, 1e15, ρ_drip, 1e15],
            r = [2e4, 2e4, 5e3, 5e3, 2e4, Rcci],
            Beb_core = 0.2,
            Bj_core = 0.3,
        )
        BA(ρ) = ρ > ρ_drip ? 0.1 : 0.05
        Beb(ρ) = ρ > ρ_drip ? 0.2 : 0.15
        Bj(ρ) = ρ > ρ_drip ? 0.3 : 0.25
        B_itp = (BA, Beb, Bj)

        Bs = MutualFrictionCoefficients(Param, B_itp; ρ_b = ρ_drip, R_cci = Rcci)
        @test Bs.BA ≈ [0.05, 0.1, 0.05, 0.1, 0.05, 0.1]
        @test Bs.Beb ≈ [0.15, 0.2, 0.15, 0.2, 0.15, 0.2]
        @test Bs.Bj ≈ [0.25, 0.3, 0.25, 0.3, 0.25, 0.3]
        crust_below_drip = (Param.ρ .< ρ_drip) .& (Param.r .>= Rcci)
        @test all(==(0.05), Bs.BA[crust_below_drip])
        @test all(==(0.15), Bs.Beb[crust_below_drip])
        @test all(==(0.25), Bs.Bj[crust_below_drip])
    end

    @testset "coefficient interpolation integration" begin
        output = @test_logs (:info, r"JSON data successfully loaded") VNparaGraber2018(
            input_file,
        )
        Param =
            (ρ = output.ρ, r = fill(2e4, length(output.ρs)), Beb_core = 0.2, Bj_core = 0.3)

        Bs = MutualFrictionCoefficients(Param, output.B_itp; ρ_b = 0.0, R_cci = 0.0)

        @test Bs.BA ≈ output.B.BA
        @test Bs.Beb ≈ output.B.Beb
        @test Bs.Bj ≈ output.B.Bj
    end

    @testset "CGS and SI profile equivalence" begin
        BA_itp(ρs) = -1.0 + 0.1 * (log10.(ρs)-15.0)
        Beb_itp(ρs) = -2.0 + 0.2 * (log10.(ρs)-15.0)
        Bj_itp(ρs) = -3.0 + 0.3 * (log10.(ρs)-15.0)

        B_itp = (BA_itp, Beb_itp, Bj_itp)
        common = (Beb_core = 0.2, Bj_core = 0.3)
        Param_si = (; common..., ρ = [1e13, 1e15], r = [2e4, 2e4])
        Param_cgs = (; common..., ρ = Param_si.ρ ./ 1e3, r = Param_si.r .* 1e2)

        Bs_si = MutualFrictionCoefficients(Param_si, B_itp; input_units = "SI")
        Bs_cgs = MutualFrictionCoefficients(Param_cgs, B_itp; input_units = "CGS")

        @test Bs_cgs.BA ≈ Bs_si.BA
        @test Bs_cgs.Beb ≈ Bs_si.Beb
        @test Bs_cgs.Bj == Bs_si.Bj
        @test_throws ArgumentError MutualFrictionCoefficients(
            Param_si,
            B_itp;
            input_units = "invalid",
        )
    end
end
