@testset "EoS" begin
    # This begin a set of tests for the EoS
    @testset "EoS_LInterp - Matrix input" begin
        # Making density array (unsorted to check the sorting part of the scripts)
        rho = [1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 9.1, 8.8, 10.0]
        # Quadratic EoS
        press = rho .^ 2
        # Making an input array for the test
        input = [rho press]
        # Creating equation of state and its inverse function
        EoS, EoS_inv = QSpin.TOV.EquationOfState.EoS_LInterp(input, (1, 2))
        # Comparing
        @test isapprox(EoS(rho), press)
        @test isapprox(EoS_inv(press), rho)
    end
    @testset "EoS_LInterp - Path input" begin
        file_name = joinpath(@__DIR__, "test.dat")
        EoS, EoS_inv = QSpin.TOV.EquationOfState.EoS_LInterp(file_name, (1, 2));
        rho = [1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 9.1, 8.8, 10.0]*1e12
        press = EoS(rho)
        @test isapprox(EoS_inv(press), rho, rtol = 0.02)
    end
end
