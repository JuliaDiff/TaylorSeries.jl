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
certain cases, it provides a primitive CAS (computer algebra system),
which works numerically and not symbolically.
The package allows to manipulate polynomials of a specified maximum
degree, including power and composition, as well as series expansions
of some elementary functions on polynomials, e.g. `exp`,
where techniques of automatic differentiation are used
[@Tucker:ValidatedNumerics; @HaroEtAl:ParameterizMeth]. Differentiation and
integration are also implemented.

Two basic immutable types are defined, `Taylor1` and `TaylorN`,
which represent the series expansions in one or several variables,
respectively. These are essentially vectors of coefficients,
ordered by increasing degree. In the case of `TaylorN`, the
coefficients are `HomogeneousPolynomials`, which in turn are vectors
of coefficients representing all monomials with a given number of variables
and order (total degree), in some lexicographical order.

Julia's parametrized type system allows the construction of Taylor series whose coefficient type is any subtype of the `Number` abstract type. Use cases include complex numbers,
arbitrary precision `BigFloat`s [@MPFR],
`Interval`s [@IntervalArithmetic.jl], `ArbFloat`s [@ArbFloats.jl],
as well as `Taylor1` and `TaylorN` objects themselves.

`TaylorSeries.jl` is the main component of
[`TaylorIntegration.jl`](https://github.com/PerezHz/TaylorIntegration.jl)
[@TaylorIntegration.jl], whose aim is to perform accurate integration
of ODEs using the Taylor method, including jet transport techniques,
where a small region around an initial condition is integrated.
It is also a key component of
[`TaylorModels.jl`](https://github.com/JuliaIntervals/TaylorModels.jl)
[@TaylorModels.jl], whose aim is to construct rigorous polynomial
approximations of functions.

# Examples

We present three examples to illustrate the use of `TaylorSeries.jl`. Other
examples, as well as a detailed user guide, can be found in the
[documentation](http://www.juliadiff.org/TaylorSeries.jl/stable).

## Hermite polynomials
As a first example we describe how to generate the [Hermite polynomials][@HermitePols_wikipedia]
("physicist's" version) up to a given maximum order. Firstly we directly exploit the recurrence relation satisfied by the polynomials.

```julia
julia> using TaylorSeries

julia> displayBigO(false)

julia> function hermite_polynomials(::Type{T}, nmax::Int) where {T <: Integer}

           x = Taylor1(T, nmax)    # Taylor variable
           H = fill(x, nmax + 1)   # vector of Taylor series to be overwritten

           H[1] = 1   # order 0
           H[2] = 2x  # order 1

           for n in 2:nmax
               # recursion relation for order n:
               H[n+1] = 2x * H[n] - 2(n-1) * H[n-1]
           end

           return H
     end

julia> hermite_polynomials(n) = hermite_polynomials(Int, n);

julia> H = hermite_polynomials(10);

julia> function hermite_polynomial(n::Int)
    @assert 0 ≤ n ≤ length(H) "Not enough Hermite polynomials generated"
    return H[n+1]
end

julia> hermite_polynomial(6)  # degree 6
- 120 + 720 t² - 480 t⁴ + 64 t⁶

```

The example above can be slightly modified to compute, for example, the 100th Hermite polynomial.
In this case, the coefficients will be larger than $2^63-1$, so the modular
behavior under overflow of the standard `Int64` type will not suffice. Rather, the polynomials should
be generated with `hermite_polynomials(BigInt, 100)` to ensure
the use of arbitrary length integers.

## Using a generating function
As a second example, we describe a numerical way of obtaining the
Hermite polynomials from their generating function: the $n$-th Hermite polynomial
corresponds to the $n$-th derivative of the function $\exp(2t \, x - t^2)$.

```julia
julia> 𝒢(x,t) = exp(2t * x - t^2)  # generating function; 𝒢 is typed as \scrG<TAB>

julia> xn = set_variables("x", numvars=1, order=10)

julia> x = xn[1]

julia> t = Taylor1([zero(x), one(x)], 10)  # Taylor1{TaylorN{Float64}}

julia> gf = 𝒢(x, t)  # Taylor1 expansion of 𝒢

julia> HH(n::Int) = derivative(n, gf)  # n-th derivative of `gf`

julia> HH(6)
- 120.0 + 720.0 x₁² - 480.0 x₁⁴ + 63.99999999999999 x₁⁶
```

This example shows that the calculations are performed numerically and not
symbolically, using `TaylorSeries.jl` as a polynomial manipulator; this
 is manifested by the fact that the last coefficient of `HH(6)` is not
 identical to an integer.

## Taylor method for integrating ordinary differential equations
As a final example, we give a simple implementation of Picard
iteration to integrate an ordinary differential equation, which is equivalent to
the Taylor method.

We consider the initial-value problem $\dot{x} = x$,
with initial condition $x(0) = 1$. One step of the integration corresponds
to constructing the Taylor series of the solution $x(t)$ in powers of $t$:

 ```julia

julia> ∫⬩dt(u::Taylor1) = integrate(u) # the symbol ∫ is obtained as \int<TAB>

julia> function taylor_step(f, u0)

           u = copy(u0)
           unew = u0 + ∫⬩dt(f(u))

           while unew != u
               u = unew
               unew = u0 + ∫⬩dt(f(u))   # Picard iteration
           end

           return u
       end

julia> f(x) = x  # Differential equation

julia> order = 20  # maximum order of the Taylor expansion for the solution

julia> u0 = Taylor1([1.0], order)  # initial condition given as a Taylor expansion

julia> solution = taylor_step(f, u0);   # solution

julia> solution(1.0) - exp(1.0) # compare the solution evaluated at t=1 with the exact value
0.0

```

Thus this Taylor expansion of order 20 around $t_0=0$
suffices to obtain the exact solution at $t=1$, while the error at time $t=2$
from the same expansion is $4.53 \times 10^{-14}$.
This indicates that a proper treatment should try to estimate what size step
should be taken as a function of the behavior of the solution.

--
## Acknowledgements

We are thankful for the additions of
[all contributors](https://github.com/JuliaDiff/TaylorSeries.jl/graphs/contributors)
to this project. We acknowledge financial support from PAPIME grants
PE-105911 and PE-107114, and PAPIIT grants IG-101113, IG-100616
and IN-117117. LB and DPS acknowledge support via the Cátedra Marcos Moshinsky (2013 and 2018, respectively).

# References
