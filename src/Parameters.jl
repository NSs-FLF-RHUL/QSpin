"""
Defines the `ParameterType` used by `parameters` arguments to QSpin equation-of-motion functions.
"""
module Parameters

export ParameterType

"Static type for `parameter` arguments that QSpin functions take."
const ParameterType = Union{NamedTuple,Tuple}

end
