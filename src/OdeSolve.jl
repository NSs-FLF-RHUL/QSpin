module OdeSolve
using LinearAlgebra
using FFTW
using MAT
using ParallelStencil

"""
Integrate an equation of motion using the Runge-Kutta 4-th order method.

:param u: The target solution of the equation of motion.
:param δt: The integrating time step.
:param time: Current time.
:param eom: The equation of motion of the problem.
"""
function ode_rk4(
    u::AbstractArray, 
    δt::Float64, 
    time::Float64, 
    eom::Function
)
    k1 = eom(u, time);
    k2 = eom(u+0.5*k1*δt, time+0.5*δt);
    k3 = eom(u+0.5*k2*δt, time+0.5*δt);
    k4 = eom(u+k3*δt, time+δt);
    un = u + δt * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    return un
end


"""
Time-evolve an equation of motion using the RK4 Runge-Kutta 4-th order method.

:param ψ0: Initial value for the field at time 0.
:param dt: Timestep interval (used as the integral timestep in RK4).
:param Dt: Time interval between recorded field values.
:param t_end: End time for equation evolution.
:param eom: The equation of motion of the problem.
:returns ψall: Field values at recorded timestamps.
    `ψall[.., i]` is the field value at time `tspan[i]`.
:returns tspan: Timestamps at which field values were recorded.
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
        "Field is ",
        dims,
        "D dimensional. Time slices will be along dimension ",
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

    while t < t_end
        ψcurrent = ode_rk4(ψcurrent, dt, t, eom)
        t += dt
        step_number += 1

        if mod(step_number, ΔNt) == 0
            println("t=", t)
            selectdim(ψall, time_dimension_index, save_number + 1) .= ψcurrent
            tspan[save_number+1] = t
            save_number += 1
        end
    end
    return ψall, tspan
end

end
