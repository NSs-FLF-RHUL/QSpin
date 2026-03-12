using QSpin
using Test

@testset "ode_rk4" begin
    # This begins a set of tests for the ode_rk4 function

    @testset "ode_rk4 No Evolution" begin
        # This begins a subset of tests for the ode_rk4 function,
        # where we will use a time-independent field to confirm some simple properties.

        # A simple "EoM" that doesn't evolve the field (dψdt = 0)
        no_evolution(ψ, time) = zeros(size(ψ))
        # An arbitrary initial field, that we don't expect to evolve
        initial_field = [1.0; 2.0]

        # Confirm that there is no evolution according to RK4
        @test isapprox(
            QSpin.OdeSolve.ode_rk4(initial_field, 1.0, 0.0, no_evolution),
            initial_field,
        )
        # When there is no evolution, the timestep δt should not matter
        @test isapprox(
            QSpin.OdeSolve.ode_rk4(initial_field, 1.0, 0.0, no_evolution),
            QSpin.OdeSolve.ode_rk4(initial_field, 10.0, 0.0, no_evolution),
        )
    end
end

@testset "CartGrid" begin
    # This begins a set of tests for the CartGrid function
    @testset "1D grid setup" begin
        x, kx, facx = QSpin.Grids.CartGrid([10.0], [4]);
        @test isapprox(x, [-10.0, -5.0, 0.0, 5.0])
        @test isapprox(kx, [0, 1, -2., -1.].*(pi/10))
        @test isapprox(facx, pi / 10)
    end
    @testset "2D grid setup" begin
        x, y, X, Y, kx, ky, Kx, Ky, facx, facy = QSpin.Grids.CartGrid([10.0, 5.0], [4,10]);
        @test isapprox(x, [-10.0, -5.0, 0.0, 5.0])
        @test isapprox(y, range(-5.0, stop=4, length=10))
        @test isapprox(X[1,:],X[2,:]) # Check that X is constant along the y direction
        @test isapprox(Y[:,1],Y[:,2]) # Check that Y is constant
        @test isapprox(kx, [0, 1, -2., -1.].*(pi/10))
        @test isapprox(ky, [range(0,stop=4,length=5) ; range(-5,stop=-1,length=5)].*(pi/5))
        @test isapprox(Kx[1,:],Kx[2,:]) # Check that Kx is constant along the y direction
        @test isapprox(Ky[:,1],Ky[:,2]) # Check that Ky is constant along the x direction
        @test isapprox(Ky[:,1],ky) # Check that Ky is consistent with ky
        @test isapprox(Kx[1,:],kx) # Check that Kx is consistent with kx
        @test isapprox(facx, kx[2])
        @test isapprox(facy, ky[2])
    end
    @testset "3D grid setup" begin
        x, y, z, X, Y, Z, kx, ky, kz, Kx, Ky, Kz, facx, facy, facz = QSpin.Grids.CartGrid([10.0, 5.0, 20.], [4,10, 20]);
        @test isapprox(x, [-10.0, -5.0, 0.0, 5.0])
        @test isapprox(y, range(-5.0, stop=4, length=10))
        @test isapprox(z, range(-20.0, stop=18, length=20))
        @test isapprox(kx, [0, 1, -2., -1.].*(pi/10))
        @test isapprox(ky, [range(0,stop=4,length=5) ; range(-5,stop=-1,length=5)].*(pi/5))
        @test isapprox(kz, [range(0,stop=9,length=10) ; range(-10,stop=-1,length=10)].*(pi/20))
         @test isapprox(Kx[1,:,1],Kx[2,:,5]) # Check that Kx is constant along the y direction
        @test isapprox(Ky[:,1,1],Ky[:,2,3]) # Check that Ky is constant along the x direction
        @test isapprox(Ky[:,1,1],ky) # Check that Ky is consistent with ky
        @test isapprox(Kx[1,:,1],kx) # Check that Kx is consistent with kx
        @test isapprox(facx, kx[2])
        @test isapprox(facy, ky[2])
        @test isapprox(facz, kz[2])
    end
end
 