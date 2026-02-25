a = [1 2 3; 4 5 6]

typeof(a)

if typeof(a) == Matrix{Int64} 
    format = println(typeof(a))
end

# data defined beforehand...
using GLMakie

points = Observable(data[:, 1])

fig, ax = scatter(points)
limits!(ax, 1, size(data)[1], 0, max(data...))

frames = 2:(size(data)[2])

record(fig, "animation.mp4", frames;
        framerate = 30) do frame
    points[] = data[:, frame]
end