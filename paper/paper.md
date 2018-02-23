---
title: 'TaylorSeries.jl: Taylor expansions in one and many variables in Julia'
tags:
  - Taylor series
  - Automatic differentiation
  - Julia
authors:
 - name: Luis Benet
   orcid: 0000-0002-8470-9054
   affiliation: 1
 - name: David P. Sanders
   affiliation: 2
affiliations:
 - name: Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)
   index: 1
 - name: Departamento de Física, Facultad de Ciencias, Universidad Nacional Autónoma de México (UNAM)
   index: 2
date: February 23, 2018
bibliography: paper.bib
---

# Summary

[`TaylorSeries.jl`](https://github.com/JuliaDiff/TaylorSeries.jl)
[@TaylorSeries] provides a framework to use and manipulate
Taylor series polynomials in one and more independent variables
in [Julia](https://julialang.org) [@Julia]. It allows to compute
elementary functions of
polynomials (`Taylor1`- or `TaylorN`-type objects), where
techniques of automatic differentiation are exploited
[@Tucker:ValidatedNumerics, @HaroEtAl:ParameterizMeth]. Derivation,
integration is also implemented.

The package allows to work with a different `Number` formats
as coefficients of the series, including complex numbers,
the extended precision `BigFloat`s [@MPFR],
`Intervals` [@ValidatedNumerics], `ArbFloats` [@ArbFloats],
as well as `Taylor1` and `TaylorN` objects.

`TaylorSeries.jl` is a core component of
[`TaylorIntegration.jl`](https://github.com/PerezHz/TaylorIntegration.jl)
[@TaylorIntegration], whose aim is to perform accurate integration
of ODEs using Taylor's method, including Jet transport techniques.

## Acknowledgements

We are thankful for the additions of
[all contributors](https://github.com/JuliaDiff/TaylorSeries.jl/graphs/contributors)
to this project. We acknowledge financial support from PAPIME grants
PE-105911 and PE-107114, and PAPIIT grants IG-101113 and IG-100616. LB
acknowledges support through a Cátedra Marcos Moshinsky (2013).

# References
