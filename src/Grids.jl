module Grids
using FFTW: plan_fft
using DocStringExtensions: TYPEDSIGNATURES

using ..Parameters: ParameterType

"""
$(TYPEDSIGNATURES)

Set up a uniform Cartesian grid, applicable for the Fourier spectral method.

# Arguments
- `CompDomain::AbstractVector{<:Real}`: The half computational domain size.
    Input as an array for [Lx,Ly,Lz], and up to 3D.
- `GridSize::AbstractVector{<:Integer}`: The number of grid points in each dimension.
    Input as an array of the form [Nx,Ny,Nz], and up to 3D.
"""
function CartGrid(CompDomain::AbstractVector{<:Real}, GridSize::AbstractVector{<:Integer})
    dims = length(GridSize)
    dims in 1:3 || throw(ArgumentError("CartGrid supports one to three dimensions"))
    length(CompDomain) == dims ||
        throw(DimensionMismatch("CompDomain and GridSize must have equal lengths"))
    all(>(0), CompDomain) || throw(ArgumentError("domain sizes must be positive"))
    all(>(0), GridSize) || throw(ArgumentError("grid sizes must be positive"))

    if dims <= 3
        Nx = GridSize[1]
        dx = 2 * CompDomain[1] / Nx
        x = range(-CompDomain[1], stop = CompDomain[1]-dx, length = Int(Nx))
        kx = [
            0:fld(Nx-1, 2);
            (-fld(Nx, 2)):-1
        ]
        facx = pi / CompDomain[1]
        kx = kx .* facx
    end
    if dims >= 2
        Ny = GridSize[2]
        dy = 2 * CompDomain[2] / Ny
        y = range(-CompDomain[2], stop = CompDomain[2]-dy, length = Int(Ny))
        ky = [
            0:fld(Ny-1, 2);
            (-fld(Ny, 2)):-1
        ]
        facy = pi / CompDomain[2]
        ky = ky .* facy
    end
    if dims == 3 # In progress, not tested yet
        Nz = GridSize[3]
        dz = 2 * CompDomain[3] / Nz
        z = range(-CompDomain[3], stop = CompDomain[3]-dz, length = Int(Nz))
        kz = [
            0:fld(Nz-1, 2);
            (-fld(Nz, 2)):-1
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
$(TYPEDSIGNATURES)

Return a function that computes the quantum kinetic energy term in Schrodinger-type equations.

Namely, return a function that computes ``-\\nabla^2 \\psi``, given ``\\psi``, using the Fourier spectral method.
The returned function has signature `kinetic_energy(ψ, parameters, time)`.

The k-square matrix can be obtained from the `CartGrid` function in `Grids.jl`.
The plans for the forward and inverse Fourier transforms can be generated using the `plan_{i}fft` function(s) provided by `FFTW`.

# Arguments
- `KE_mtx::AbstractArray`: k-square matrix for computing kinetic energy in momentum space.
    This can be obtained by using the k-matrices in `CartGrid` function in `Grids.jl`.
- `PFFT::`: Plan for forward FFT.
- `PiFFT::`: Plan for inverse FFT.

# Returns
- `kinetic_energy::Function`: Function that evaluates the kinetic energy.
"""
function Pfft_ke(KE_mtx::AbstractArray{Float64}, PFFT, PiFFT)
    function kinetic_energy(ψ::AbstractArray, parameters::ParameterType, time::Float64)
        return PiFFT * (KE_mtx .* (PFFT * ψ))
    end

    return kinetic_energy
end

"""
$(TYPEDSIGNATURES)

Return a function that computes the quantum kinetic energy term in Schrodinger-type equations.

Namely, return a function that computes ``-\\nabla^2 \\psi``, given ``\\psi``, using the Fourier spectral method.
The returned function has signature `kinetic_energy(ψ, parameters, time)`.

The k-square matrix can be obtained from the `CartGrid` function in `Grids.jl`.

# Arguments
- `KE_mtx::AbstractArray`: k-square matrix for computing kinetic energy in momentum space.
    This can be obtained by using the k-matrices in `CartGrid` function in `Grids.jl`.

# Returns
- `kinetic_energy::Function`: Function that evaluates the kinetic energy.
"""
function fft_ke(KE_mtx::AbstractArray{Float64})
    PFFT = plan_fft(KE_mtx)
    PiFFT = inv(PFFT)

    return Pfft_ke(KE_mtx, PFFT, PiFFT)
end

"""
$(TYPEDSIGNATURES)

Return a function that computes the quantum angular momentum using the Fourier spectral method.

The angular momentum operator along the z-axis is given by

``
\\mathcal{L}_z = -\\mathrm{i} \\left(x \\frac{\\partial}{\\partial y} - y \\frac{\\partial}{\\partial x} \\right).
``

The returned function has signature `Lz(ψ, parameters, time)`.

The various coordinate matrices can be obtained via the `CartGrid` function in `Grids.jl`.

# Arguments
- `X::AbstractMatrix{Float64}`: x-coordinate matrix.
- `Y::AbstractMatrix{Float64}`: y-coordinate matrix.
- `Kx::AbstractMatrix{Float64}`: kx-coordinate matrix.
- `Ky::AbstractMatrix{Float64}`: ky-coordinate matrix.
"""
function fft_Lzψ(
    X::AbstractMatrix{Float64},
    Y::AbstractMatrix{Float64},
    Kx::AbstractMatrix{Float64},
    Ky::AbstractMatrix{Float64},
)

    PFFT = plan_fft(X)
    PiFFT = inv(PFFT)
    return Pfft_Lzψ(X, Y, Kx, Ky, PFFT, PiFFT)
end

"""
$(TYPEDSIGNATURES)

Return a function that computes the quantum angular momentum using the Fourier spectral method.

The angular momentum operator along the z-axis is given by

``
\\mathcal{L}_z = -\\mathrm{i} \\left(x \\frac{\\partial}{\\partial y} - y \\frac{\\partial}{\\partial x} \\right).
``

The returned function has signature `Lz(ψ, parameters, time)`.

The various coordinate matrices can be obtained via the `CartGrid` function in `Grids.jl`.
The plans for the forward and inverse Fourier transforms can be generated using the `plan_{i}fft` function(s) provided by `FFTW`.

# Arguments
- `X::AbstractMatrix{Float64}`: x-coordinate matrix.
- `Y::AbstractMatrix{Float64}`: y-coordinate matrix.
- `Kx::AbstractMatrix{Float64}`: kx-coordinate matrix.
- `Ky::AbstractMatrix{Float64}`: ky-coordinate matrix.
- `PFFT::`: Plan for forward FFT.
- `PiFFT::`: Plan for inverse FFT.
"""
function Pfft_Lzψ(
    X::AbstractMatrix{Float64},
    Y::AbstractMatrix{Float64},
    Kx::AbstractMatrix{Float64},
    Ky::AbstractMatrix{Float64},
    PFFT,
    PiFFT,
)
    function angular_momentum_z(
        ψ::Union{AbstractArray,Array{Float64},Array{ComplexF64}},
        parameters::ParameterType,
        time::Float64,
    )
        ψk = PFFT * ψ;
        return -(Y .* (PiFFT * (Kx .* ψk)) - X .* (PiFFT * (Ky .* ψk)));
    end
    return angular_momentum_z
end

end
