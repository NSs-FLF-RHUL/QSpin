@testset "fft_ke" begin
    @testset "1D fft_ke test" begin
        # Creating the grids
        x, kx, facx = QSpin.Grids.CartGrid([10.0], [4]);
        f = cos.(2*facx .* x); # A test function with a known kinetic energy
        KE_f = QSpin.Grids.fft_ke(0.5 .* kx .^ 2)(f, 0.0); # The kinetic energy of the test function
        KE_f_analytical = 0.5 * (2*facx)^2 * cos.(2*facx .* x); # The analytical kinetic energy of the test function
        # Checking the x coordingates, the kx coordinates, and the facx value
        @test isapprox(KE_f, KE_f_analytical);
    end
    @testset "2D fft_ke test" begin
        # Creating the grids
        x, y, X, Y, kx, ky, Kx, Ky, facx, facy = QSpin.Grids.CartGrid([10.0, 5.0], [6, 10]);
        f = cos.(3*facx .* X) .* cos.(2*facy .* Y);# A test function with a known kinetic energy
        KE_f = QSpin.Grids.fft_ke(0.5 .* (Kx .^ 2 + Ky .^ 2))(f, 0.0); # The kinetic energy of the test function
        KE_f_analytical = 0.5 * ((3*facx)^2+(2*facy)^2) * f; # The analytical kinetic energy of the test function
        # Checking the x coordingates, the kx coordinates, and the facx value
        @test isapprox(real(KE_f[:]), KE_f_analytical[:]);
    end
    @testset "3D fft_ke test" begin
        # Creating the grids
        x, y, z, X, Y, Z, kx, ky, kz, Kx, Ky, Kz, facx, facy, facz =
            QSpin.Grids.CartGrid([10.0, 5.0, 20.0], [6, 10, 20]);
        f = cos.(3*facx .* X) .* cos.(2*facy .* Y) .* cos.(facz .* Z);# A test function with a known kinetic energy
        KE_f = QSpin.Grids.fft_ke(0.5 .* (Kx .^ 2 + Ky .^ 2 + Kz .^ 2))(f, 0.0); # The kinetic energy of the test function
        KE_f_analytical = 0.5 * ((3*facx)^2+(2*facy)^2+facz^2) * f; # The analytical kinetic energy of the test function
        # Checking the x coordingates, the kx coordinates, and the facx value
        @test isapprox(real(KE_f[:]), KE_f_analytical[:]);
    end
end
