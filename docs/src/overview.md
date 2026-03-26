# [Quick Overview](@id using-desolve)

`QSpin` relies on [`OrdinaryDiffEq`](https://github.com/SciML/OrdinaryDiffEq.jl) to time-evolve equations of motion.
The `QSpin` package itself provides some helper methods for constructing frequently-used terms (or functions, or sub-functions) that appear in the equations of motion, as well as some thin wrappers around the functionality we need from `OrdinaryDiffEq`.

One important thing to keep in mind is that the functions `QSpin` provides need to be written in a format that `OrdinaryDiffEq.ODEProblem` is expecting.
The full documentation for these problems is available `OrdinaryDiffEq`'s website, but the main impact on how we write functions in `QSpin` is the functional form of the equation of motion $f$.

The function $f$ that defines the equation of motion should have one of two signatures:

- The simplest signature is `f(ψ, parameters, t)`.
  - These functions should explicitly `return` an array that corresponds to the evaluation of the equation of motion $f$.
  - `ψ` and `t` are the field value and current time respectively.
  - `parameters` is a container that includes all constant parameter values that are needed to evaluate the equation of motion.
- A more memory-efficient format for the functions is `f!(dψ, ψ, parameters, t)`.
  - `dψ` **should be overwritten** with the result of evaluating $f$, rather than having the result of the evaluation explicitly returned via `return`.
  - `ψ`, `t`, and `parameters` have the same interpretations as above.

## QSpin Conventions

### Helper Functions for the Equations of Motion

There are several functions within `QSpin` that can generate an appropriate function `f` from the component "pieces" of the equation of motion.
The `QSpin.Hamiltonian.hamiltonian!` function, for example, constructs a function `H!(dψ, ψ, parameters, t)` from two other functions that compute the kinetic and potential energy.
`QSpin` uses the convention that the names of these "constructor" functions have an `!` at the end if **the function they return** has the `f!(dψ, ψ, parameters, t)` signature - that is, when the returned function is expected to mutate its input in-place.
Conversely, constructor functions without the `!` at the end of their name return functions with the `f(ψ, parameters, t)` signature.

### `parameters` Data Type

`OrdinaryDiffEq` recommends using immutable containers for the `parameters` variable.
As such `QSpin` elects to store these values in `NamedTuple`s:

```julia
my_parameters = (
  parameter_name_1 = 1.0,
  parameters_name_2 = 3.14,
)
```

## Time-Evolving Equations of Motion

`QSpin.OdeSolve` provides the `evolve` method which is a thin wrapper that creates an appropriate [`OrdinaryDiffEq.ODEProblem`](https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/#SciMLBase.ODEProblem) and then passes it to [`OrdinaryDiffEq.solve`](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#CommonSolve.solve-Tuple{SciMLBase.AbstractDEProblem,%20Vararg{Any}}).
Solver options, usually taken by `OrdinaryDiffEq.solve`, can be passed via keyword arguments to `solver_options`.
The positional arguments to [`QSpin.OdeSolve.evolve`](@ref) are essentially analogues of the arguments that are passed to `OrdinaryDiffEq.ODEProblem` to create the problem.

### Example: Coupled System of ODEs

As a simple example, let's solve the system of equations

```math
\begin{aligned}
\frac{\mathrm{d}\psi}{\mathrm{d}t} = M \psi,
&\qquad
\psi_0 = \begin{pmatrix} 0.1 \\ 0.2 \end{pmatrix},
\end{aligned}
```

over the time interval $t\in[0,1]$, with

```math
M = \begin{pmatrix} -2 & 1 \\ 1 & -2 \end{pmatrix}.
```

Naively, we could just define our equation of motion function $f$ as

```julia
function f(ψ)
  return [-2 1; 1 -2] * ψ
end
```

however, this is not the format that [`evolve`](@ref QSpin.OdeSolve.evolve) is expecting.
It also `return`s the evaluation explicitly, rather than updating an argument in place.
To use this equation of motion with `evolve`, we need to give it the appropriate signature:

```@setup coupled-odes
using LaTeXStrings
using Plots
using OrdinaryDiffEq
using QSpin

default(show = false)

function create_plot!(ψ, plot_ref = nothing, label_1 = L"$\psi_1$", label_2 = L"$\psi_2$"; fmt_args...)
  if plot_ref === nothing
    plot_ref = plot();
  end

  plot!(
    plot_ref, ψ.t, ψ[1, :], label = label_1;
    fmt_args...,
  );
  plot!(
    plot_ref, ψ.t, ψ[2, :], label = label_2; 
    fmt_args...,
  );

  return plot_ref
end

default_plot_args = (
  xlabel = L"$t$",
  ylabel = L"$\psi$",
)
```

```@example coupled-odes; continued = true
function f!(dψ, ψ, parameters, t)
  # Perform an in-place overwrite of the input array,
  # replacing the current values with the result of evaluation.
  dψ .= [-2 1; 1 -2] * ψ
end
```

Notice that:

- Even though the equation of motion does not use the `parameters` or `t` arguments, they must be present in the function signature.
- The final line performs an in-place update (`.=`) to one of the input arrays.
- The exclamation mark at the end of the function name (`f!`) is consistent with Julia's naming practices for functions that edit their input arguments.

We could then time-evolve $f$ as follows:

```@example coupled-odes
ψ0 = [0.1; 0.2]
timestep = 1e-3
save_interval = 1e-1

# Note that the default time-range for evolve is [0, 1].
ψ = QSpin.OdeSolve.evolve(
  f!, ψ0;
  alg=OrdinaryDiffEq.Tsit5(),
  dt=timestep,
  saveat=save_interval,
)

create_plot!(ψ; default_plot_args...); # hide
```

### Example: Parametrised Coupling Matrix

In the example above, we hard-coded the value we wanted the matrix $M$ to take into our equation of motion.
This is inefficient, as every time we want to solve an equation of motion of the same form as $f$, we'd need to re-write (or define a new function) `f!`.
We can instead make use of `OrdinaryDiffEq`'s parametrisation functionality, and the `parameters` argument to `evolve`, to avoid having multiple "versions" of the same function floating around.

```@example coupled-odes
function f_parametrised!(dψ, ψ, parameters, t)
  dψ .= parameters.M * ψ
end

ψ0 = [0.1; 0.2]
timestep = 1e-3
save_interval = 1e-1
original_parameters = (
  M = [-2 1; 1 -2],
)
alternative_parameters = (
  M = [-4 1; 1 -4],
)

ψ_original = QSpin.OdeSolve.evolve(
  f_parametrised!, ψ0, p=original_parameters;
  alg=OrdinaryDiffEq.Tsit5(),
  dt=timestep,
  saveat=save_interval,
)
ψ_alternative = QSpin.OdeSolve.evolve(
  f_parametrised!, ψ0, p=alternative_parameters;
  alg=OrdinaryDiffEq.Tsit5(),
  dt=timestep,
  saveat=save_interval,
)

plot_ref = create_plot!( # hide
  ψ_original, # hide
  nothing, # hide
  L"$(\psi_{\mathrm{original}})_1$", # hide
  L"$(\psi_{\mathrm{original}})_2$"; # hide
  default_plot_args... # hide
); # hide
plot_ref = create_plot!( # hide
  ψ_alternative, # hide
  plot_ref, # hide
  L"$(\psi_{\mathrm{alternative}})_1$", # hide
  L"$(\psi_{\mathrm{alternative}})_2$"; # hide
  default_plot_args... # hide
); # hide
```
