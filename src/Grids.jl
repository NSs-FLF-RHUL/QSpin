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
function CartGrid(CompDomain::Array{Float64},GridSize::Array{Int64})
    dims = length(GridSize)
    if dims <=2
        Nx = GridSize[1]
        dx   = 2 * CompDomain[1]/ Nx
        x    = range(-CompDomain[1],stop=CompDomain[1]-dx,length=Int(Nx))
        kx   = [range(0,stop=Nx/2-1,length=Int(Nx/2)) ;range(-Nx/2,stop=-1,length=Int(Nx/2))]
        facx = pi / CompDomain[1]
    end
    if dims <=3
        Ny = GridSize[2]
        dy   = 2 * CompDomain[2] / Ny
        y    = range(-CompDomain[2],stop=CompDomain[2]-dy,length=Int(Ny))
        ky   = [range(0,stop=Ny/2-1,length=Int(Ny/2)) ;range(-Ny/2,stop=-1,length=Int(Ny/2))]
        facy = pi / CompDomain[2]

        X   = repeat(x',Ny,1)
        Y   = repeat(y,1,Nx)
        Kx  = repeat(kx',Ny,1)
        Ky  = repeat(ky,1,Nx);
    end
    if dims == 3 # In progress, not tested yet
        Nz = GridSize[3]
        dz   = 2 * CompDomain[3] / Nz
        z    = range(-CompDomain[2],stop=CompDomain[2]-dz,length=Int(Nz))
        kz   = [range(0,stop=Nz/2-1,length=Int(Nz/2)) ;range(-Nz/2,stop=-1,length=Int(Nz/2))]
        facz = pi / CompDomain[2]
    end 
    if dims == 1
        println("Creating 1D Cartesian grid with ", Nx, " points.")
        return x, kx, facx
    elseif dims == 2
        println("Creating 2D Cartesian grid with ", Nx, " x ", Ny, " points.")
        return x, y, X, Y, kx, ky,Kx, Ky, facx, facy
    elseif dims == 3
        println("Creating 3D Cartesian grid with ", Nx, " x ", Ny, " x ", Nz, " points.")
        return x, y, z, X, Y, Z, Kx, Ky, Kz, facx, facy, facz
    end
end 



function fft_ke(KE_mtx::Array{Float64})
    function kinetic_energy(ψ::Union{AbstractArray,Array{Float64},Array{ComplexF64}})
        return ifft(KE_mtx.*fft(ψ))
    end
    return kinetic_energy
end

function fft_Lzψ(X::Array{Float64},Y::Array{Float64},Kx::Array{Float64}, Ky::Array{Float64})
    function angular_momentum_z(ψ::Union{AbstractArray,Array{Float64},Array{ComplexF64}})
        return im .* (Y .* ifft(im .* Kx .*fft(ψ)) - X .* ifft(im .* Ky .*fft(ψ)))
    end
    return angular_momentum_z
end
end