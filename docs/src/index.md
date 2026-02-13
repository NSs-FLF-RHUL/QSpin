# My Package Documentation

```@docs
QSpin.OdeSolve.ode_rk4
```

## Building the Documentation

From the top-level package directory, run the `docs/make.jl` script using the `docs` project;

```bash
julia --project=docs docs/make.jl
```

This will create the `docs/build` folder, containing the built documentation.
The syntax for automatic documentation can be found on `Documenter.jl`'s website: <https://documenter.juliadocs.org/stable/man/syntax/>

## Linting

The codebase [uses `JuliaFormatter.jl` for its code linting style](https://github.com/domluna/JuliaFormatter.jl).
This is the default linter that is packaged with the VSCode Julia extension.

The project is configured with a `pre-commit` hook that will format (and attempt to fix in place) any `.jl` files that have been added to git's staging area.
To setup the `pre-commit` checks, you will first need to [install `pre-commit`](https://`pre-commit`.com/#installation), which will require a suitable Python environment (unless you want to install `pre-commit` system wide).
Once `pre-commit` is installed, you can run the usual

```bash
pre-commit install
```

within the project's root directory to configure the hooks.
