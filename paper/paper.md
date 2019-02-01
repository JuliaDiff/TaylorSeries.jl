---
title: 'TaylorSeries.jl: Taylor expansions in one and several variables in Julia'
tags:
  - Taylor series
  - Automatic differentiation
  - Julia
authors:
 - name: Luis Benet
   orcid: 0000-0002-8470-9054
   affiliation: 1
 - name: David P. Sanders
   orcid: 0000-0001-5593-1564
   affiliation: 2
affiliations:
 - name: Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)
   index: 1
 - name: Departamento de Física, Facultad de Ciencias, Universidad Nacional Autónoma de México (UNAM)
   index: 2
date: September 20, 2018
bibliography: paper.bib
---

# Summary

`TaylorSeries.jl`
provides a framework to use and manipulate
Taylor polynomials in one and more variables
in the [Julia programming language](https://julialang.org) [@julia]. It allows to compute
elementary functions, e.g. `exp`, of
polynomials (`Taylor1`- or `TaylorN`-type objects), where
techniques of automatic differentiation are used
[@Tucker:ValidatedNumerics; @HaroEtAl:ParameterizMeth]. Differentiation and
integration are also implemented.

The package allows to work with different `Number` formats
as coefficients of the series, including complex numbers,
the arbitrary precision `BigFloat`s [@MPFR],
`Interval`s [@ValidatedNumerics.jl], `ArbFloat`s [@ArbFloats.jl],
as well as `Taylor1` and `TaylorN` objects.

`TaylorSeries.jl` is a core component of
[`TaylorIntegration.jl`](https://github.com/PerezHz/TaylorIntegration.jl)
[@TaylorIntegration.jl], whose aim is to perform accurate integration
of ODEs using the Taylor method, including jet transport techniques,
and of
[`TaylorModels.jl`](https://github.com/JuliaIntervals/TaylorModels.jl)
[@TaylorModels.jl],
which constructs rigorous polynomial approximations of functions.

## Acknowledgements

We are thankful for the additions of
[all contributors](https://github.com/JuliaDiff/TaylorSeries.jl/graphs/contributors)
to this project. We acknowledge financial support from PAPIME grants
PE-105911 and PE-107114, and PAPIIT grants IG-101113, IG-100616
and IN-117117. LB acknowledges support through a Cátedra Marcos Moshinsky (2013).

# References
