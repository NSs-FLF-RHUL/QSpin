"""
Integrating an equation of motion usin the Runge-Kutta 4-th order method

:param u: the targed solution of the equation of motion.
:param δt: the Integrating time step.
:param time: time.
:param eom: the equation of motion of the problem.
"""
function ode_rk4(u::Array{Float64},dt::Float64,time::Float64,eom::function)   
    k1 = eom(u,time);
    k2 = eom(u+0.5*k1*dt,time+0.5*dt); 
    k3 = eom(u+0.5*k2*dt,time+0.5*dt); 
    k4 = eom(u+k3*dt,time+dt);
    un = u + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    return un
end