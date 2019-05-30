# TaylorSeries.jl

A [Julia](http://julialang.org) package for Taylor polynomial expansions in one or more
independent variables.

[![Build Status](https://api.travis-ci.org/JuliaDiff/TaylorSeries.jl.svg?branch=master)](https://travis-ci.org/JuliaDiff/TaylorSeries.jl)
[![Coverage Status](https://coveralls.io/repos/JuliaDiff/TaylorSeries.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaDiff/TaylorSeries.jl?branch=master)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](http://www.juliadiff.org/TaylorSeries.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](http://www.juliadiff.org/TaylorSeries.jl/latest)

[![DOI](http://joss.theoj.org/papers/10.21105/joss.01043/status.svg)](https://doi.org/10.21105/joss.01043)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2601941.svg)](https://zenodo.org/record/2601941)

#### Authors
- [Luis Benet](http://www.cicc.unam.mx/~benet/), Instituto de Ciencias Físicas,
Universidad Nacional Autónoma de México (UNAM)
- [David P. Sanders](http://sistemas.fciencias.unam.mx/~dsanders/), Facultad
de Ciencias, Universidad Nacional Autónoma de México (UNAM)

Comments, suggestions and improvements are welcome and appreciated.

#### Examples
Taylor series in one varaible
```julia
julia> using TaylorSeries

julia> t = Taylor1(Float64, 5)
 1.0 t + 𝒪(t⁶)

julia> exp(t)
 1.0 + 1.0 t + 0.5 t² + 0.16666666666666666 t³ + 0.041666666666666664 t⁴ + 0.008333333333333333 t⁵ + 𝒪(t⁶)
 
 julia> log(1 + t)
 1.0 t - 0.5 t² + 0.3333333333333333 t³ - 0.25 t⁴ + 0.2 t⁵ + 𝒪(t⁶)
 ```
Multivariate Taylor series
 ```julia
julia> x, y = set_variables("x y", order=2);

julia> exp(x + y)
 1.0 + 1.0 x + 1.0 y + 0.5 x² + 1.0 x y + 0.5 y² + 𝒪(‖x‖³)
 
```
Differential and integral calculus on Taylor series:
```julia
julia> x, y = set_variables("x y", order=4);

julia> p = x^3 + 2x^2 * y - 7x + 2
 2.0 - 7.0 x + 1.0 x³ + 2.0 x² y + 𝒪(‖x‖⁵)

julia> ∇(p)
2-element Array{TaylorN{Float64},1}:
  - 7.0 + 3.0 x² + 4.0 x y + 𝒪(‖x‖⁵)
                    2.0 x² + 𝒪(‖x‖⁵)

julia> integrate(p, 1)
 2.0 x - 3.5 x² + 0.25 x⁴ + 0.6666666666666666 x³ y + 𝒪(‖x‖⁵)

julia> integrate(p, 2)
 2.0 y - 7.0 x y + 1.0 x³ y + 1.0 x² y² + 𝒪(‖x‖⁵)
```

For more details, please see the [docs](http://www.juliadiff.org/TaylorSeries.jl/stable).

#### License

`TaylorSeries` is licensed under the [MIT "Expat" license](./LICENSE.md).

#### Installation

`TaylorSeries` can be installed simply with `using Pkg; Pkg.add("TaylorSeries")`.

#### Contributing

There are many ways to contribute to this package:

- Report an issue if you encounter some odd behavior, or if you have suggestions to improve the package.
- Contribute with code addressing some open issues, that add new functionality or that improve the performance.
- When contributing with code, add docstrings and comments, so others may understand the methods implemented.
- Contribute by updating and improving the documentation.

#### References

- W. Tucker, Validated numerics: A short introduction to rigorous
computations, Princeton University Press (2011).
-  A. Haro, Automatic differentiation methods in computational dynamical
systems: Invariant manifolds and normal forms of vector fields at fixed points,
[preprint](http://www.maia.ub.es/~alex/admcds/admcds.pdf).

#### Acknowledgments
This project began (using `python`) during a Masters' course in the postgraduate
programs in Physics and in Mathematics at UNAM, during the second half of 2013.
We thank the participants of the course for putting up with the half-baked
material and contributing energy and ideas.

We acknowledge financial support from DGAPA-UNAM PAPIME grants
PE-105911 and PE-107114, and DGAPA-PAPIIT grants IG-101113
and IG-100616.
LB acknowledges support through a *Cátedra Moshinsky* (2013).
