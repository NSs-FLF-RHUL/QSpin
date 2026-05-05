module OdeSolve
using LinearAlgebra
using FFTW
using MAT
using ParallelStencil
using DocStringExtensions: TYPEDSIGNATURES

import OrdinaryDiffEq as DE

using ..Parameters: ParameterType

"""
$(TYPEDSIGNATURES)

Time-evolve an equation of motion using `OrdinaryDiffEq` (DE).

This is a thin wrapper around `DE.ODEProblem` and `DE.solve`. The keyword arguments in
`solver_options` are handed directly to `DE.solve`. See
[their documentation](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#solver_options)
for the full range of options.

The equation of motion provided is expected to have signature `(dψ, ψ, parameters, t)`, where
`dψ` is the array to which the output is to be written in-place. `ψ`, `parameters`, and `t`
are the current field, parameters of the problem, and current time respectively.

# Arguments
- `eom!::Function`: Equation of motion.
-  `ψ0::AbstractArray`: Initial field value.
- `t_start::Float64`: Start time for system evolution.
- `t_end::Float64`: End time for system evolution.
- `p::NamedTuple`: Problem parameters, to pass to `DE.ODEProblem`.
- `solver_options`: Keyword arguments that will be passed to `DE.solve`.

# Returns
- `solution`: [Solution object](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#CommonSolve.solve-Tuple{SciMLBase.AbstractDEProblem,%20Vararg{Any}}) for the problem that was time-evolved.
"""
function evolve(
    eom!::Function,
    ψ0::AbstractArray,
    t_start::Float64 = 0.0,
    t_end::Float64 = 1.0,
    p::Union{ParameterType,Nothing} = nothing;
    solver_options...,
)
    if p === nothing
        p = ()
    end

    problem = DE.ODEProblem(eom!, ψ0, (t_start, t_end), p)
    return DE.solve(problem; solver_options...)
end

"""
$(TYPEDSIGNATURES)

Integrate an equation of motion using the [Runge-Kutta 4-th order method](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods).

# Arguments
- `u::AbstractArray`: The target solution of the equation of motion.
- `δt::Float64`: The integrating time step.
- `time::Float64`: Current time.
- `eom::Function`: The equation of motion of the problem.

# Returns
- `un::AbstractArray`: Field value at time `time + δt`.
"""
function ode_rk4(u::AbstractArray, δt::Float64, time::Float64, eom::Function)
    k1 = eom(u, time);
    k2 = eom(u+0.5*k1*δt, time+0.5*δt);
    k3 = eom(u+0.5*k2*δt, time+0.5*δt);
    k4 = eom(u+k3*δt, time+δt);
    un = u + δt * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    return un
end

"""
$(TYPEDSIGNATURES)

Time-evolve an equation of motion using the RK4 Runge-Kutta 4-th order method.

# Arguments
- `ψ0::AbstractArray`: Initial value for the field at time 0.
- `dt::Float64`: Timestep interval (used as the integral timestep in RK4).
- `Dt::Float64`: Time interval between recorded field values.
- `t_end::Float64`: End time for equation evolution.
- `eom::Function`: The equation of motion of the problem.

# Returns
- `ψall::AbstractArray`: Field values at recorded timestamps.
    `ψall[..., i]` is the field value at time `tspan[i]`.
- `tspan::Vector{Float64}`: Timestamps at which field values were recorded.
"""
function evolve_rk4(
    ψ0::Union{AbstractArray,Array{Float64},Array{ComplexF64}},
    dt::Float64,
    Dt::Float64,
    t_end::Float64,
    eom::Function,
)
    dims = ndims(ψ0)
    time_dimension_index = dims + 1
    println(
        "。  Solving ",
        dims,
        "D dimensional EOM. Time slices will be along dimension ",
        time_dimension_index,
        ".",
    )

    ΔNt = floor(Int, Dt / dt)
    Nt = floor(Int, t_end / Dt)
    ψall = zeros(eltype(ψ0), size(ψ0)..., Nt + 1)

    selectdim(ψall, time_dimension_index, 1) .= ψ0
    tspan = zeros(Nt + 1)

    t = 0.0
    ψcurrent = ψ0
    save_number = 1
    step_number = 0
    println("。  Simulation Begins。")
    println(" t = ", t)
    @inbounds while t < t_end
        ψcurrent = ode_rk4(ψcurrent, dt, t, eom)
        t += dt
        step_number += 1
        if sum(isnan.(ψcurrent[:]))>0
            println(
                "NaN detected in the field at time ",
                t,
                ". Time Step could be too big.",
            )
            break
        end
        if mod(step_number, ΔNt) == 0
            println(" t = ", t)
            selectdim(ψall, time_dimension_index, save_number + 1) .= ψcurrent
            tspan[save_number+1] = t
            save_number += 1
        end
    end
    return ψall, tspan
end

end
