# QSpin

[![Build Status](https://github.com/willGraham01/MyJuliaPackage.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/willGraham01/MyJuliaPackage.jl/actions/workflows/CI.yml?query=branch%3Amain)

## About

:warning: This package is currently under construction and in pre-release.
The API and features may change suddenly without warning.

### Project Team

Vanessa Graber ([vanessa.graber@rhul.ac.uk](mailto:vanessa.graber@rhul.ac.uk)),

Gary Liu ([Gary.Liu@rhul.ac.uk](mailto:Gary.Liu@rhul.ac.uk)),

Will Graham ([william.graham@ucl.ac.uk](mailto:william.graham@ucl.ac.uk))

<!-- TODO: how do we have an array of collaborators - steal from s2fft -->

## Description

QSpin contains multi-functional PDE solvers for quantum fluid systems.
It is specifically designed for solving [Gross-Pitaevskii equation](https://iopscience.iop.org/article/10.3847/1538-4357/adc383), [Ginzburg-Landau model](https://www.mdpi.com/2218-1997/8/4/228), and similar non-linear Schrödinger equations for their stationary/equilibrium states and dynamics.
The former is approached by the steepest descent and imaginary-time propagation methods.
The time integration here is using Runge-Kutta methods, and adaptive time step is enable.

Here we also include an example script for a neutron star glitch simulator, based on a three-component model, [Graber et al. (2018)](http://arxiv.org/abs/1804.02706), and a python version can be accessed at [glitchsim](https://github.com/NSs-FLF-RHUL/glitchsim.git), which is based on the Jupyter notebook in [glitchraiser](https://github.com/vanessagraber/glitchrises.git).

## Getting Started

1. Download [Julia](https://julialang.org/downloads/) (version 1.12 or later).

2. Install essentially required packages, including [FFTW](https://github.com/JuliaMath/FFTW.jl), [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/), [MAT](https://juliaio.github.io/MAT.jl/stable/), [OrdinaryDiffEq](https://docs.sciml.ai/OrdinaryDiffEq/stable/):

- Julia version 1.12 or higher.
- Install required Julia packages:

  ```julia
  using Pkg
  Pkg.add(["LinearAlgebra","OrdinaryDiffEq","FFTW","MAT"])

If you are going to use GPU feature, pleae furtehr include

```julia
Pkg.add("ParallelStencil")

For data visualization, we primarily use [Plots](https://docs.juliaplots.org/stable/) and [LaTeXStrings](https://github.com/JuliaStrings/LaTeXStrings.jl) to plot results in our scirpts.

3. Download [QSpin](https://github.com/NSs-FLF-RHUL/QSpin/tree/gl_glitch_raiser_script?tab=readme-ov-file) package.

## Functions and Features in QSpin

There are 4 modules, OdeSolve, Grids, Hamiltonian, Parameters, in QSpin package. They are built to solve the targeted differential equation problems.

OdeSolve:

Grids: Geranting the spatial and momentum grids (CartGrid), computing quantum kinetic energy (fft_ke & Pfft_ke) and angular momenta (fft_Lz).

Hamiltonian: Combining kinetic and potential energy functions in Schrdinger-like equations.

Parameters:

## Usage



## Running your first PDE/ODE
