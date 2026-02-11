using QSpin
using Plots

# Setting the EoM problem
   
function eom(ψ::Array{Float64},time::Float64)
    A = [-2 1;1 -2]
    dψdt = A * ψ
end

# Running and Plotting
@time ut, t = QSpin.OdeSolve.evolve1d([.1;.2],1e-3,1e-2,10.,eom)
plot(t,ut[1,:])
plot!(t,ut[2,:],xlabel = "time",ylabel="displacement",title="Solving a set of coupled ODE")
