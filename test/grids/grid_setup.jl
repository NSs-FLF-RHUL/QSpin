@testset "CartGrid" begin
    # This begins a set of tests for the CartGrid function
    @testset "1D grid setup" begin
        # This begins a subset of test for the CartGrid function,
        # where we will test the 1D grid setup

        # Creating the grids
        x, kx, facx = QSpin.Grids.CartGrid([10.0], [4]);

        # Checking the x coordingates, the kx coordinates, and the facx value
        @test isapprox(x, [-10.0, -5.0, 0.0, 5.0])
        @test isapprox(kx, [0, 1, -2.0, -1.0] .* (pi/10))
        @test isapprox(facx, kx[2])
    end
    @testset "2D grid setup" begin
        # This begins a subset of test for the CartGrid function,
        # where we will test the 2D grid setup

        # Creating the grids
        x, y, X, Y, kx, ky, Kx, Ky, facx, facy = QSpin.Grids.CartGrid([10.0, 5.0], [4, 10]);

        # Checking the x and y coordingates, the kx and ky coordinates, the Kx and Ky meshgrid, and the facx and facy values
        @test isapprox(x, [-10.0, -5.0, 0.0, 5.0])
        @test isapprox(y, range(-5.0, stop = 4, length = 10))
        @test isapprox(X[1, :], X[2, :]) # X array should be the same along y direction.
        @test isapprox(Y[:, 1], Y[:, 2]) # Y array should be the same along x direction.
        @test isapprox(Y[:, 1], y) # Y array should be the same as y along y direction.
        @test isapprox(X[1, :], x) # X array should be the same as x along x direction.
        @test isapprox(kx, [0, 1, -2.0, -1.0] .* (pi/10))
        @test isapprox(
            ky,
            [range(0, stop = 4, length = 5); range(-5, stop = -1, length = 5)] .* (pi/5),
        )
        @test isapprox(Kx[1, :], Kx[2, :]) # Kx array should be the same along y direction.
        @test isapprox(Ky[:, 1], Ky[:, 2]) # Ky array should be the same along x direction.
        @test isapprox(Ky[:, 1], ky) # Ky array should be the same as ky along y direction.
        @test isapprox(Kx[1, :], kx) # Kx array should be the same as kx along x direction.
        @test isapprox(facx, kx[2])
        @test isapprox(facy, ky[2])
        @test size(X)==size(Y)
        @test size(Kx)==size(Ky)
        @test size(X)==(10, 4)
    end
    @testset "3D grid setup" begin
        # This begins a subset of test for the CartGrid function,
        # where we will test the 2D grid setup

        # Creating the grids
        x, y, z, X, Y, Z, kx, ky, kz, Kx, Ky, Kz, facx, facy, facz =
            QSpin.Grids.CartGrid([10.0, 5.0, 20.0], [4, 10, 20]);
        # Checking the x, y, and z coordingates, the kx, ky, and kz coordinates, the Kx and Ky meshgrid, and the facx, facy, and facz values
        @test isapprox(x, [-10.0, -5.0, 0.0, 5.0])
        @test isapprox(y, range(-5.0, stop = 4, length = 10))
        @test isapprox(z, range(-20.0, stop = 18, length = 20))
        @test isapprox(kx, [0, 1, -2.0, -1.0] .* (pi/10))
        @test isapprox(
            ky,
            [range(0, stop = 4, length = 5); range(-5, stop = -1, length = 5)] .* (pi/5),
        )
        @test isapprox(
            kz,
            [range(0, stop = 9, length = 10); range(-10, stop = -1, length = 10)] .*
            (pi/20),
        )
        @test isapprox(Kx[1, :, 1], Kx[2, :, 5]); # Kx array should be the same along y and z directions.
        @test isapprox(Ky[:, 1, 1], Ky[:, 2, 3]); # Ky array should be the same along x and z directions.
        @test isapprox(Kz[1, 2, :], Kz[3, 4, :]); # Kz array should be the same along x and y directions.
        @test isapprox(Kx[1, :, 1], kx); # Kx array should be the same as kx along x direction.
        @test isapprox(Ky[:, 1, 1], ky); # Ky array should be the same as ky along y direction.
        @test isapprox(Kz[1, 2, :], kz); # Kz array should be the same as kz along z direction.
        @test isapprox(facx, kx[2]);
        @test isapprox(facy, ky[2]);
        @test isapprox(facz, kz[2]);
        @test size(X)==size(Y)
        @test size(X)==size(Z)
        @test size(Kx)==size(Ky)
        @test size(Kx)==size(Kz)
        @test size(X)==(10, 4, 20)
        @test size(Y)==(10, 4, 20)
        @test size(Z)==(10, 4, 20)
        @test size(Kx)==(10, 4, 20)
        @test size(Ky)==(10, 4, 20)
        @test size(Kz)==(10, 4, 20)
    end
end
