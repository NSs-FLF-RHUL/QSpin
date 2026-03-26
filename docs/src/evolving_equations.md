# [Time-Evolving Equations](@id using-desolve)

`QSpin` relies on [`OrdinaryDiffEq`](https://github.com/SciML/OrdinaryDiffEq.jl) to time-evolve equations of motion.
The `QSpin` package itself provides some helper methods for constructing frequently-used terms (or functions, or sub-functions) that appear in the equations of motion, as well as some thin wrappers around the functionality we need from `OrdinaryDiffEq`.

One important thing to keep in mind is that the functions `QSpin` provides need to be written in a format that `OrdinaryDiffEq.ODEProblem` is expecting.
The full documentation for these problems is available `OrdinaryDiffEq`'s website, but the main impact on how we write functions in `QSpin` is the functional form of the equation of motion $f$.

The function $f$ that defines the equation of motion should have Julia signature `f(du, u, parameters, t)`.
`du` **should be overwritten** with the result of evaluating `f`, rather than having the result of the evaluation explicitly returned via `return`.
`u` and `t` are the field value and current time respectively.
`parameters` is a container that includes all constant parameter values that are needed to evaluate the equation of motion.
`OrdinaryDiffEq` recommends using immutable containers for the `parameters` variable, and as such `QSpin` elects to store these values in `NamedTuple`s.

## Time-Evolving Equations of Motion

`QSpin.OdeSolve` provides the `evolve` method which is a thin wrapper that creates an appropriate [`OrdinaryDiffEq.ODEProblem`](https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/#SciMLBase.ODEProblem) and then passes it to [`OrdinaryDiffEq.solve`](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#CommonSolve.solve-Tuple{SciMLBase.AbstractDEProblem,%20Vararg{Any}}).
Solver options, usually taken by `OrdinaryDiffEq.solve`, can be passed via keyword arguments to `solver_options`.
The positional arguments to `QSpin.OdeSolve.evolve` are essentially analogues of the arguments that are passed to `OrdinaryDiffEq.ODEProblem` to create the problem.

```@docs
QSpin.OdeSolve.evolve
```

## Example: Coupled System of ODEs

As a simple example, let's solve the system of equations

```math
\begin{aligned}
\frac{\mathrm{d}\psi}{\mathrm{d}t} = \begin{pmatrix} -2 & 1 \\ 1 & -2 \end{pmatrix} \psi,
&\qquad
\psi_0 = \begin{pmatrix} 0.1 \\ 0.2 \end{pmatrix},
\end{aligned}
```

over the time interval $t\in[0,1]$.
