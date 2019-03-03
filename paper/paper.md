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
which works numerically and not symbolically.
The package allows to manipulate polynomials of a specified maximum
degree, including power and composition, as well as series expansions
of some elementary functions on polynomials, e.g. `exp`,
where techniques of automatic differentiation are used
[@Tucker:ValidatedNumerics; @HaroEtAl:ParameterizMeth]. Differentiation and
integration are also implemented.

Two basic immutable types are defined, `Taylor1` and `TaylorN`,
which represent the series expansions in one or several variables,
respectively. These structures are essentially vectors of coefficients
ordered increasingly by its degree. In the case of `TaylorN`, the
coefficients are `HomogeneousPolynomials`, which in turn are vectors
of coefficients representing all monomials of given number of variables
and order (total degree), in some lexicographical order.
The package allows to work with different `Number` formats
for the coefficients of the series, including complex numbers,
arbitrary precision `BigFloat`s [@MPFR],
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

# Examples

We present three examples to show the use of `TaylorSeries.jl`. Other
examples as well as a detailed user guide can be found in the
[documentation](http://www.juliadiff.org/TaylorSeries.jl/stable/).

As a first example we describe how to generate the [Hermite polynomials][@HermitePols]
("physicists" version), up to a maximum order. We begin exploiting directly
the recurrence relation satisfied by the polynomials.

```julia
julia> using TaylorSeries

julia> displayBigO(false);

julia> function generate_hermite_polynomials(::Type{T}, nmax::Int) where {T<:Integer}
  Hn = Vector{Taylor1{T}}(undef, nmax+1)
  Hn[1] = one(T)+Taylor1(T, 0) # order 0
  Hn[2] = 2*Taylor1(T, 1)      # order 1
  for n in 2:nmax
      x = Taylor1(Int128, n) # Taylor variable of n-th order
      # Recursion formula for degree n
      Hn[n+1] = 2*x*Hn[n] - 2*(n-1)*Hn[n-1]
  end
  return Hn
end;

julia> generate_hermite_polynomials(n) = generate_hermite_polynomials(Int, n);

julia> Hn = generate_hermite_polynomials(10);

julia> function HermitePol(n::Int)
    @assert 0 ≤ n ≤ length(Hn) "Not enough Hermite polynomials generated"
    return Hn[n+1]
end;

julia> HermitePol(6)
- 120 + 720 t² - 480 t⁴ + 64 t⁶

```

The example above can be slightly modified to compute the 100th Hermite polynomial.
In this case, the coefficients will be larger than `2^63-1`, so the modular
`Int64` arithmetic will not suffice. In this case, the polynomials should
be generated with `generate_hermite_polynomials(BigInt, 100)` to ensure
using extended integer precision.

As a second example, we shall describe a *numeric* form of obtaining the
Hermite polynomials from the generating function. The n-th Hermite polynomial
corresponds to the n-th derivative of the Taylor expansion with respect to `t`
of the function `exp(2t*x-t^2)`.

```julia
julia> 𝒢(x,t) = exp(2*t*x-t^2); # generating function; 𝒢 is typed as \scrG<TAB>

julia> xn = set_variables("x", numvars=1, order=10);

julia> x = xn[1];

julia> t = Taylor1([zero(x),one(x)], 10); # Taylor1{TaylorN{Float64}}

julia> gf = 𝒢(0*t+x, t+0*x); # Taylor expansion of 𝒢

julia> Hnn(n::Int) = derivative(n, gf); # n-th derivative of `gf`

julia> Hnn(6)
- 120.0 + 720.0 x₁² - 480.0 x₁⁴ + 63.99999999999999 x₁⁶
```

This example shows that the calculations are performed numerically and not
simbolically, manifested by the fact that the last coefficient is not an
integer. This example is aimed to describe the possibility offered
by `TaylorSeries.jl` as a polynomial manipulator.

--
## Acknowledgements

We are thankful for the additions of
[all contributors](https://github.com/JuliaDiff/TaylorSeries.jl/graphs/contributors)
to this project. We acknowledge financial support from PAPIME grants
PE-105911 and PE-107114, and PAPIIT grants IG-101113, IG-100616
and IN-117117. LB acknowledges support through a Cátedra Marcos Moshinsky (2013).

# References
