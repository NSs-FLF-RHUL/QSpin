module OdeSolve

using CommonSolve: solve
using DocStringExtensions: TYPEDSIGNATURES
using OrdinaryDiffEqLowOrderRK: OrdinaryDiffEqLowOrderRK
using SciMLBase: ODEProblem

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
- `alg`: ODE algorithm; defaults to `DP5()`.
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
    alg = OrdinaryDiffEqLowOrderRK.DP5(),
    solver_options...,
)
    if p === nothing
        p = ()
    end

    problem = ODEProblem(eom!, ψ0, (t_start, t_end), p)
    return solve(problem, alg; solver_options...)
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

    dt > 0 || throw(ArgumentError("dt must be positive"))
    Dt > 0 || throw(ArgumentError("Dt must be positive"))
    t_end >= 0 || throw(ArgumentError("t_end must be non-negative"))

    save_every = round(Int, Dt / dt)
    save_every >= 1 || throw(ArgumentError("Dt must be at least dt"))
    isapprox(save_every * dt, Dt; rtol = 1e-12, atol = 1e-14) ||
        throw(ArgumentError("Dt must be an integer multiple of dt"))

    step_count = floor(Int, t_end / dt + 1e-12)
    save_count = fld(step_count, save_every) + 1
    ψall = zeros(eltype(ψ0), size(ψ0)..., save_count)

    selectdim(ψall, time_dimension_index, 1) .= ψ0
    tspan = zeros(save_count)

    t = 0.0
    ψcurrent = ψ0
    save_number = 1
    println("。  Simulation Begins。")
    println(" t = ", t)
    @inbounds for step_number = 1:step_count
        ψcurrent = ode_rk4(ψcurrent, dt, t, eom)
        t = step_number * dt
        if any(isnan, ψcurrent)
            println(
                "NaN detected in the field at time ",
                t,
                ". Time Step could be too big.",
            )
            valid_fields = copy(selectdim(ψall, time_dimension_index, 1:save_number))
            return valid_fields, tspan[1:save_number]
        end
        if mod(step_number, save_every) == 0
            println(" t = ", t)
            save_number += 1
            selectdim(ψall, time_dimension_index, save_number) .= ψcurrent
            tspan[save_number] = t
        end
    end
    return ψall, tspan
end

end
