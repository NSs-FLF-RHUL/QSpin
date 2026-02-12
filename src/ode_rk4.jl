"""
Integrating an equation of motion usin the Runge-Kutta 4-th order method

:param u: the targed solution of the equation of motion.
:param δt: the Integrating time step.
:param time: time.
:param eom: the equation of motion of the problem.
"""
function ode_rk4(u::Array{Float64},δt::Float64,time::Float64,eom::Function)
    k1 = eom(u,time);
    k2 = eom(u+0.5*k1*δt,time+0.5*δt); 
    k3 = eom(u+0.5*k2*δt,time+0.5*δt); 
    k4 = eom(u+k3*δt,time+δt);
    un = u + δt * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    return un
end