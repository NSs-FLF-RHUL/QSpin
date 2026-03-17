module Grids
using LinearAlgebra
using FFTW
using MAT
using ParallelStencil

"""
    CartGrid(CompDomain::Array{Float64},GridSize::Array{Int64})

    Setting up uniform Cartesian grids and applicable for Fourier spectral method.

        :param CompDomain: The half computational domain size. Input as an array for [Lx,Ly,Lz] and up to 3D.
        :param GridSize: The number of grid points in each dimension, in the form of [Nx,Ny,Nz] and up to 3D.
"""
function CartGrid(CompDomain::Array{Float64}, GridSize::Array{Int64})
    dims = length(GridSize)
    if dims <= 3
        Nx = GridSize[1]
        dx = 2 * CompDomain[1] / Nx
        x = range(-CompDomain[1], stop = CompDomain[1]-dx, length = Int(Nx))
        kx = [
            range(0, stop = Nx/2-1, length = Int(Nx/2));
            range(-Nx/2, stop = -1, length = Int(Nx/2))
        ]
        facx = pi / CompDomain[1]
        kx = kx .* facx
    end
    if dims >= 2
        Ny = GridSize[2]
        dy = 2 * CompDomain[2] / Ny
        y = range(-CompDomain[2], stop = CompDomain[2]-dy, length = Int(Ny))
        ky = [
            range(0, stop = Ny/2-1, length = Int(Ny/2));
            range(-Ny/2, stop = -1, length = Int(Ny/2))
        ]
        facy = pi / CompDomain[2]
        ky = ky .* facy
    end
    if dims == 3 # In progress, not tested yet
        Nz = GridSize[3]
        dz = 2 * CompDomain[3] / Nz
        z = range(-CompDomain[3], stop = CompDomain[3]-dz, length = Int(Nz))
        kz = [
            range(0, stop = Nz/2-1, length = Int(Nz/2));
            range(-Nz/2, stop = -1, length = Int(Nz/2))
        ]
        facz = pi / CompDomain[3]
        kz = kz .* facz
    end
    if dims == 1
        println("    Creating 1D Cartesian grid with ", Nx, " points.")
        return x, kx, facx
    elseif dims == 2
        X = repeat(x', Ny, 1)
        Y = repeat(y, 1, Nx)
        Kx = repeat(kx', Ny, 1);
        Ky = repeat(ky, 1, Nx);
        println("    Creating 2D Cartesian grid with ", Nx, " x ", Ny, " points.")
        return x, y, X, Y, kx, ky, Kx, Ky, facx, facy
    elseif dims == 3
        X = repeat(x', GridSize[2], 1, GridSize[3]);
        Y = repeat(y, 1, GridSize[1], GridSize[3]);
        Z = permutedims(repeat(z, 1, GridSize[1], GridSize[2]), [3 2 1]);
        Kx = repeat((kx)', GridSize[2], 1, GridSize[3])
        Ky = repeat((ky), 1, GridSize[1], GridSize[3])
        Kz = permutedims(repeat((kz), 1, GridSize[1], GridSize[2]), [3 2 1]); # Meshgrid K_sq
        println(
            "    Creating 3D Cartesian grid with ",
            Nx,
            " x ",
            Ny,
            " x ",
            Nz,
            " points.",
        )
        return x, y, z, X, Y, Z, kx, ky, kz, Kx, Ky, Kz, facx, facy, facz
    end
end


"""
    fft_ke(KE_mtx::Array{Float64})(ψ::Union{AbstractArray,Array{Float64},Array{ComplexF64}})

    Computing the quantum kinetic energy term in Schrodinger-type equations, namely, -∇^2ψ, using the Fourier spectral method.

    :param KE_mtx: the k-square matrix for computing kinetic energy in momentum space, which can be obtained by using the k-matrices in CartGrid function in Grids.jl.
    :param ψ: the field for computing.

"""
function fft_ke(KE_mtx::Array{Float64})
    function kinetic_energy(
        ψ::Union{AbstractArray,Array{Float64},Array{ComplexF64}},
        time::Float64,
    )
        return ifft(KE_mtx .* fft(ψ))
    end
    return kinetic_energy
end

"""
    fft_Lzψ(X::Array{Float64},Y::Array{Float64},Kx::Array{Float64}, Ky::Array{Float64})(ψ::Union{AbstractArray,Array{Float64},Array{ComplexF64}})

    Computing the quantum angular momentum term, namely, Lz ψ, using the Fourier spectral method.
    The angular momentum operator along the z-axis is given by Lz = -i (x ∂/∂y - y ∂/∂x).

    :param X: the x-coordinate matrix, which can be obtained by using the x-matrix in CartGrid function in Grids.jl.
    :param Y: the y-coordinate matrix, which can be obtained by using the y-matrix in CartGrid function in Grids.jl.
    :param Kx: the kx-coordinate matrix, which can be obtained by using the kx-matrix in CartGrid function in Grids.jl.
    :param Ky: the ky-coordinate matrix, which can be obtained by using the ky-matrix in CartGrid function in Grids.jl.
    :param ψ: the field for computing.

"""
function fft_Lzψ(
    X::Array{Float64},
    Y::Array{Float64},
    Kx::Array{Float64},
    Ky::Array{Float64},
)
    function angular_momentum_z(
        ψ::Union{AbstractArray,Array{Float64},Array{ComplexF64}},
        time::Float64,
    )
        ψk = fft(ψ);
        return im .* (Y .* ifft(im .* Kx .* ψk) - X .* ifft(im .* Ky .* ψk))
    end
    return angular_momentum_z
end
end
