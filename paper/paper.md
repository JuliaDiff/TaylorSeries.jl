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

The purpose of `TaylorSeries.jl` is to provide a framework to exploit
Taylor polynomials in one and more variables
in the [Julia programming language](https://julialang.org) [@julia]. In
certain cases, it provides a very primitive CAS (computer algebra system),
though it does the work numerically, and not on a symbolic level;
some examples can be found in the
[documentation](http://www.juliadiff.org/TaylorSeries.jl/stable/).
The package allows to manipulate polynomials of an specified maximum
degree, including power and composition, as well as series expansions
of some elementary functions of polynomials, e.g. `exp`,
where techniques of automatic differentiation are used
[@Tucker:ValidatedNumerics; @HaroEtAl:ParameterizMeth]. Differentiation and
integration are also implemented.

Two basic immutable types are defined, `Taylor1` and `TaylorN`,
which represent the series expansions in one of more variables,
respectively. These structures are essentially vectors of coefficients
ordered increasingly by its degree. In the case of `TaylorN`, the
coefficients are `HomogeneousPolynomials`, which in turn are vectors
of coefficients representing all monomials of given number of variables
and order (total degree), in some lexicographical ordering.
The package allows to work with different `Number` formats
as coefficients of the series, including complex numbers,
the arbitrary precision `BigFloat`s [@MPFR],
`Interval`s [@ValidatedNumerics.jl], `ArbFloat`s [@ArbFloats.jl],
as well as `Taylor1` and `TaylorN` objects.

`TaylorSeries.jl` is a core component of
[`TaylorIntegration.jl`](https://github.com/PerezHz/TaylorIntegration.jl)
[@TaylorIntegration.jl], whose aim is to perform accurate integration
of ODEs using the Taylor method, including jet transport techniques,
where a small region around an initial condition is integrated.
It is also an important component of
[`TaylorModels.jl`](https://github.com/JuliaIntervals/TaylorModels.jl)
[@TaylorModels.jl], whose aim is to construct rigorous polynomial
approximations of functions.

## Acknowledgements

We are thankful for the additions of
[all contributors](https://github.com/JuliaDiff/TaylorSeries.jl/graphs/contributors)
to this project. We acknowledge financial support from PAPIME grants
PE-105911 and PE-107114, and PAPIIT grants IG-101113, IG-100616
and IN-117117. LB acknowledges support through a Cátedra Marcos Moshinsky (2013).

# References
