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
        @test output.A ≈ output.Z .* (1 .+ 1 ./ output.x)

        for values in (
            output.n1,
            output.Rws,
            output.ρs,
            output.R[2],
            output.R[3],
            output.B[2],
            output.B[3],
        )
            @test all(isfinite, values)
            @test all(>(0), values)
        end
        @test all(<=(0.5), output.B[2])
        @test all(<=(0.5), output.B[3])
        @test all(>(0), diff(output.ρs))

        # The splines are built in log-log space and must reproduce their knots.
        @test exp10.(output.B_itp[2].(log10.(output.ρs))) ≈ output.B[2]
        @test exp10.(output.B_itp[3].(log10.(output.ρs))) ≈ output.B[3]

        # Left extrapolation is constant by construction.
        below_range = log10(first(output.ρs)) - 1
        @test exp10(output.B_itp[2](below_range)) ≈ first(output.B[2])
        @test exp10(output.B_itp[3](below_range)) ≈ first(output.B[3])

        @test_throws SystemError VNparaGraber2018(input_file * ".missing")
    end

    @testset "MutualFrictionCoefficients regions" begin
        ρ_drip = 4e14
        Rcci = 1e4
        Param = (
            ρs = [1e13, 1e15, 1e13, 1e15, ρ_drip, 1e15],
            r = [2e4, 2e4, 5e3, 5e3, 2e4, Rcci],
            Beb_core = 0.2,
            Bj_core = 0.3,
        )
        Beb_itp = Returns(-2.0)
        Bj_itp = Returns(-3.0)

        Bs =
            MutualFrictionCoefficients(Param, Beb_itp, Bj_itp; ρ_drip = ρ_drip, Rcci = Rcci)

        @test Bs.Beb ≈ [0.0, 1e-2, 0.2, 0.2, 1e-2, 1e-2]
        @test Bs.Bj ≈ [0.0, 1e-3, 0.3, 0.3, 1e-3, 1e-3]
        crust_below_drip = (Param.ρs .< ρ_drip) .& (Param.r .>= Rcci)
        @test all(==(0.0), Bs.Beb[crust_below_drip])
        @test all(==(Param.Beb_core), Bs.Beb[Param.r .< Rcci])
        @test all(==(Param.Bj_core), Bs.Bj[Param.r .< Rcci])
    end

    @testset "coefficient interpolation integration" begin
        output = @test_logs (:info, r"JSON data successfully loaded") VNparaGraber2018(
            input_file,
        )
        Param = (
            ρs = output.ρs,
            r = fill(2e4, length(output.ρs)),
            Beb_core = 0.2,
            Bj_core = 0.3,
        )

        Bs = MutualFrictionCoefficients(
            Param,
            output.Beb_itp,
            output.Bj_itp;
            ρ_drip = 0.0,
            Rcci = 0.0,
        )

        @test Bs.Beb ≈ output.B[2]
        @test Bs.Bj ≈ output.B[3]
    end

    @testset "CGS and SI profile equivalence" begin
        Beb_itp(log_ρs) = -2 + 0.1 * (log_ρs - 15)
        Bj_itp(log_ρs) = -3 + 0.1 * (log_ρs - 15)
        common = (Beb_core = 0.2, Bj_core = 0.3)
        si = (; common..., ρs = [1e13, 1e15], r = [2e4, 2e4])
        cgs = (; common..., ρs = si.ρs ./ 1e3, r = si.r .* 1e2)

        Bs_si = MutualFrictionCoefficients(si, Beb_itp, Bj_itp; input_units = "SI")
        Bs_cgs = MutualFrictionCoefficients(cgs, Beb_itp, Bj_itp; input_units = "CGS")

        @test Bs_cgs.Beb ≈ Bs_si.Beb
        @test Bs_cgs.Bj ≈ Bs_si.Bj
        @test Bs_cgs.Beb[1] == Bs_si.Beb[1] == 0.0
        @test_throws ArgumentError MutualFrictionCoefficients(
            si,
            Beb_itp,
            Bj_itp;
            input_units = "invalid",
        )
    end
end
