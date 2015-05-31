<script src='https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML'>
</script>

# TaylorSeries.jl

A [Julia](http://julialang.org) package for Taylor expansions in one or more independent variables.

---

### Authors

- [Luis Benet](http://www.cicc.unam.mx/~benet/), Instituto de Ciencias Físicas,
Universidad Nacional Autónoma de México (UNAM)
- [David P. Sanders](http://sistemas.fciencias.unam.mx/~dsanders/), Facultad de Ciencias,
Universidad Nacional Autónoma de México (UNAM)

Comments, suggestions and additions are welcome and appreciated.

### Installation

TaylorSeries.jl is a [registered package](http://pkg.julialang.org), and is
simply installed by running

```julia
julia> Pkg.add("TaylorSeries")
```

### License

TaylorSeries is licensed under the MIT "Expat" license; see
[LICENSE](https://github.com/lbenet/TaylorSeries.jl/blob/master/LICENSE.md) for
the full license text.

### Related packages

- [Polynomials.jl](https://github.com/keno/Polynomials.jl): Polynomial
manipulations
- [PowerSeries.jl](https://github.com/jwmerrill/PowerSeries.jl): Truncated
power series for Julia
- [MultiPoly.jl](https://github.com/daviddelaat/MultiPoly.jl) Sparse
multivariate polynomials in Julia

### References

- W. Tucker, Validated numerics: A short introduction to rigorous
computations, Princeton University Press (2011).
-  A. Haro, Automatic differentiation methods in computational dynamical
systems: Invariant manifolds and normal forms of vector fields at fixed points,
[preprint](http://www.maia.ub.es/~alex/admcds/admcds.pdf).

### Acknowledgments

This project began (using `python`) during a Masters' course in the postgraduate
programs in Physics and in Mathematics at UNAM, during the second half of 2013.
We thank the participants of the course for putting up with the half-baked
material and contributing energy and ideas.

We acknowledge financial support from DGAPA-UNAM PAPIME grants PE-105911 and
PE-107114, and PAPIIT grant IG-101113. LB acknowledges support through a
*Cátedra Moshinsky* (2013).
