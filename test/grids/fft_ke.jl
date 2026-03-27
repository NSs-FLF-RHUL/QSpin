@testset "fft_ke" begin

    # KE functions take parameters arguments, but in these tests we will simply define
    # non-parametrised functions for ease.
    parameters = (;)
    @testset "1D fft_ke test" begin
        # Creating the grids
        x, kx, facx = QSpin.Grids.CartGrid([10.0], [4]);

        # A test function with a known kinetic energy
        f = cos.(2*facx .* x);

        # The kinetic energy of the test function
        KE_f = QSpin.Grids.fft_ke(0.5 .* kx .^ 2)(f, parameters, 0.0);

        # The analytical kinetic energy of the test function
        KE_f_analytical = 0.5 * (2*facx)^2 * cos.(2*facx .* x);

        @test isapprox(KE_f, KE_f_analytical);
    end
    @testset "2D fft_ke test" begin
        x, y, X, Y, kx, ky, Kx, Ky, facx, facy = QSpin.Grids.CartGrid([10.0, 5.0], [6, 10]);

        f = cos.(3*facx .* X) .* cos.(2*facy .* Y);
        KE_f = QSpin.Grids.fft_ke(0.5 .* (Kx .^ 2 + Ky .^ 2))(f, parameters, 0.0);
        KE_f_analytical = 0.5 * ((3*facx)^2+(2*facy)^2) * f;

        @test isapprox(real(KE_f[:]), KE_f_analytical[:]);
    end
    @testset "3D fft_ke test" begin
        x, y, z, X, Y, Z, kx, ky, kz, Kx, Ky, Kz, facx, facy, facz =
            QSpin.Grids.CartGrid([10.0, 5.0, 20.0], [6, 10, 20]);

        f = cos.(3*facx .* X) .* cos.(2*facy .* Y) .* cos.(facz .* Z);
        KE_f = QSpin.Grids.fft_ke(0.5 .* (Kx .^ 2 + Ky .^ 2 + Kz .^ 2))(f, parameters, 0.0);
        KE_f_analytical = 0.5 * ((3*facx)^2+(2*facy)^2+facz^2) * f;

        @test isapprox(real(KE_f[:]), KE_f_analytical[:]);
    end
end
