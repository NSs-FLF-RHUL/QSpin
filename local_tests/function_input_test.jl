"""
    my_print(x::Union{Number, AbstractArray{<:Number}})

TBW
"""
function my_print(x::Union{Number, Array, AbstractArray{<:Number}})
  println("value = ", typeof(x))
  
end

