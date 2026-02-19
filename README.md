# QSpin

[![Build Status](https://github.com/willGraham01/MyJuliaPackage.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/willGraham01/MyJuliaPackage.jl/actions/workflows/CI.yml?query=branch%3Amain)

## About

:warning: This package is currently under construction and in pre-release.
The API and features may change suddenly without warning.

### Project Team

Vanessa Graber ([vanessa.graber@rhul.ac.uk](mailto:vanessa.graber@rhul.ac.uk)),
Gary Liu ([Gary.Liu@rhul.ac.uk](mailto:Gary.Liu@rhul.ac.uk))
<!-- TODO: how do we have an array of collaborators - steal from s2fft -->

## Description

QSpin package is an multi-functional PDE solvers for quantum fluid systems. It is specfically designed for solving [Gross-Pitaevskii equation](https://iopscience.iop.org/article/10.3847/1538-4357/adc383), [Ginzburg-Landau model](https://www.mdpi.com/2218-1997/8/4/228), and similar nonlinear Schrödinger equations for their stationary/equilibrium states and dynamics. The former is approached by the steepest descent and imaginary-time propgation methods. The time integration here is using Runge-Kutta methods, and adaptive time step is enable. 

Here we also include an example script for a neutron star glitch simulator, based on a three-component model, [Graber et al. (2018)](http://arxiv.org/abs/1804.02706), and a python version can be accessed at [glitchsim](https://github.com/NSs-FLF-RHUL/glitchsim.git), which is based on the Jupyter notebook in [glitchraiser](https://github.com/vanessagraber/glitchrises.git).

## Getting Started