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

The project is configured with a `pre-commit` hook that will format (and attempt to fix in place) any `.jl` files that have been added to git's staging area with `git add`.

To setup the `pre-commit` checks, you will first need to [install `pre-commit`](https://pre-commit.com/#installation), which will require a suitable Python environment (unless you want to install `pre-commit` system wide).
Once you have created a suitable Python environment, activate it, and then install pre-commit:

```bash
pip install pre-commit
```

You will then need to set up `pre-commit` for the `QSpin` project.
Navigate to the repository root for your local clone of `QSpin`, and run

```bash
pre-commit install
```

to configure the hooks.
Once this has run, `pre-commit` should now be configured to run whenever you attempt to make a new commit.
You can also tell `pre-commit` to lint all files in the repository at any time using

```bash
pre-commit run -a
```

When run, `pre-commit` will (try to) fix the offending files, which you can then `git add` and (re)try the `git commit`.
If you get `files were modified by this hook` messages, have a look at the changes and commit them if they make sense.
There are some things that pre-commit can't fix.
But you'll get some helpful output, and the rule names that are failing, which you can then search and find what it needs you to fix.

**Note:** You will need to be working within the Python environment in which you installed pre-commit in subsequent sessions, if you want `pre-commit` to continue linting.
If you are making commits and notice that `pre-commit` isn't running, check that you have activated the corresponding environment.
