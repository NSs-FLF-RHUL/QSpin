using Plots

include("../src/OdeSolve.jl")

# Running and Plotting
@time ut, t = OdeSolve.evolve1d([.1;.2],1e-3,1e-2,10.)
plot(t,ut[1,:])
plot!(t,ut[2,:],xlabel = "time",ylabel="displacement",title="Solving a set of coupled ODE")
