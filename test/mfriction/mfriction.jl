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

        for values in
            (output.n1, output.Rws, output.ρs, output.Reb, output.Rj, output.Beb, output.Bj)
            @test all(isfinite, values)
            @test all(>(0), values)
        end
        @test all(<=(0.5), output.Beb)
        @test all(<=(0.5), output.Bj)
        @test all(>(0), diff(output.ρs))

        # The splines are built in log-log space and must reproduce their knots.
        @test exp10.(output.Beb_itp.(log10.(output.ρs))) ≈ output.Beb
        @test exp10.(output.Bj_itp.(log10.(output.ρs))) ≈ output.Bj

        # Left extrapolation is constant by construction.
        below_range = log10(first(output.ρs)) - 1
        @test exp10(output.Beb_itp(below_range)) ≈ first(output.Beb)
        @test exp10(output.Bj_itp(below_range)) ≈ first(output.Bj)

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

        @test Bs.Beb ≈ output.Beb
        @test Bs.Bj ≈ output.Bj
    end
end
