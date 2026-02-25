a = [1 2 3; 4 5 6]

typeof(a)

if typeof(a) == Matrix{Int64} 
    format = println(typeof(a))
end