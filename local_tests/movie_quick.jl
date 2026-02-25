using Plots

# Define data
x = range(0, 2π, length=100)
y = sin.(x)

# Create animation
anim = @animate for i in 1:length(x)
    plot(x[1:i], y[1:i], 
         title="Sine Wave", 
         xlims=(0, 2π), 
         ylims=(-1, 1), 
         legend=false, 
         lw=2)
end

# Save as GIF
gif(anim, "local_tests/sine_wave.gif", fps=30)
# To save as mp4 instead:
# gif(anim, "sine_wave.mp4", fps=30)