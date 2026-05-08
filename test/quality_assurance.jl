using QSpin: QSpin
using Aqua: Aqua
using ExplicitImports: ExplicitImports
using Test: @testset, @test

@testset "Aqua" begin
    Aqua.test_all(QSpin)
end

@testset "ExplicitImports" begin
    @testset "Explicit Imports" begin
        @test ExplicitImports.check_no_implicit_imports(QSpin) === nothing
    end

    @testset "Import via Owner" begin
        @test ExplicitImports.check_all_explicit_imports_via_owners(QSpin) === nothing
    end

    @testset "Stale Explicit Imports" begin
        @test ExplicitImports.check_no_stale_explicit_imports(QSpin) === nothing
    end

    @testset "Qualified Accesses" begin
        @test ExplicitImports.check_all_qualified_accesses_via_owners(QSpin) === nothing
    end

    @testset "Self Qualified Accesses" begin
        @test ExplicitImports.check_no_self_qualified_accesses(QSpin) === nothing
    end
end
