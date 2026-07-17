@testset "EoS" begin
    # This begin a set of tests for the EoS
    @testset "EoS_LInterp" begin
        # Making density array (unsorted to check the sorting part of the scripts)
        rho = [1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 9.1, 8.8, 10.0]
        # Quadratic EoS
        press = rho .^ 2
        # Making an input array for the test
        file_input = [rho press]
        # Creating equation of state and its inverse function
        EoS, EoS_inv = QSpin.TOV.EquationOfState.EoS_LInterp(file_input, (1, 2))
        # Comparing
        @test isapprox(EoS(rho), press)
        @test isapprox(EoS_inv(press), rho)
    end
end
