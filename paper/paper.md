---
title: 'TaylorSeries.jl: Taylor expansions in Julia'
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
date: February 10, 2018
bibliography: paper.bib
---

# Summary

[`TaylorSeries.jl`](https://github.com/JuliaDiff/TaylorSeries.jl)
[@TaylorSeries] provides a framework to use and manipulate
Taylor series polynomials in one and more independent variables
in [Julia](https://julialang.org) [@Julia]. In order to compute
elementary functions of
polynomials (`Taylor1`- or `TaylorN`-type objects) we
exploit techniques of automatic differentiation
[@Tucker:ValidatedNumerics, @HaroEtAl:ParameterizMeth].

The package allows to work with a large number of `Number` formats
as coefficients of the series, including complex numbers,
the extended precision format `BigFloat` [@MPFR],
`Intervals` [@ValidatedNumerics], `ArbFloats` [@ArbFloats],
as well as `Taylor1` and `TaylorN` objects, which are defined
by the package.

This package is a core component of
[`TaylorIntegration.jl`](https://github.com/PerezHz/TaylorIntegration.jl)
[@TaylorIntegration], whose aim is to perform precise integration
of ODEs using Taylor's method.

## Acknowledgements

We are thankful for the additions of
[all contributors](https://github.com/JuliaDiff/TaylorSeries.jl/graphs/contributors)
to this project. We acknowledge financial support from PAPIME grants
PE-105911 and PE-107114, and PAPIIT grants IG-101113 and IG-100616. LB
acknowledges support through a Cátedra Marcos Moshinsky (2013).

# References
