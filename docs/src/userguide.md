# User guide

---

```@meta
CurrentModule = TaylorSeries
```

[TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl)
is a basic polynomial algebraic manipulator in one or more
variables; these two cases are treated separately.  Three new types are defined,
[`Taylor1`](@ref), [`HomogeneousPolynomial`](@ref) and [`TaylorN`](@ref),
which correspond to
expansions in one independent variable, homogeneous polynomials of various
variables, and the polynomial
series in many independent variables, respectively. These types are subtypes
of `AbstractSeries`, which in turn is a subtype of `Number`, and are defined
parametrically.

The package is loaded as usual:

```@repl userguide
using TaylorSeries
```

## One independent variable

Taylor expansions in one variable are represented by the [`Taylor1`](@ref) type,
which consists of a vector of coefficients (fieldname `coeffs`) and the maximum
order considered for the expansion (fieldname `order`). The
coefficients are arranged in ascending order with respect to the degree of the
monomial, so that
`coeffs[1]` is the constant term, `coeffs[2]` gives the first order term (`t^1`),
etc. This is a dense representation of the polynomial.
The order of the polynomial can be
omitted in the constructor, which is then fixed by the length of the
vector of coefficients. If the length of the vector does not correspond with
the `order`, `order` is used, which effectively truncates polynomial to degree `order`.

```@repl userguide
Taylor1([1, 2, 3],4) # Polynomial of order 4 with coefficients 1, 2, 3
Taylor1([0.0, 1im]) # Also works with complex numbers
Taylor1(ones(8), 2) # Polynomial truncated to order 2
shift_taylor(a) = a + Taylor1(typeof(a),5)  ## a + taylor-polynomial of order 5
t = shift_taylor(0.0) # Independent variable `t`
```

Note that the information about the maximum order considered is displayed
using a big-𝒪 notation.

The definition of `shift_taylor(a)` uses the method
[`Taylor1([::Type{Float64}], [order::Int64=1])`](@ref), which is a
shortcut to define the independent variable of a Taylor expansion,
of given type and order (defaults are `Float64` and `order=1`).
This is one of the easiest ways to work with the package.

The usual arithmetic operators (`+`, `-`, `*`, `/`, `^`, `==`) have been
extended to work with the [`Taylor1`](@ref) type, including promotions that
involve `Number`s. The operations return a valid Taylor expansion of
maximum order. This is apparent in the last example below, where
the answer is beyond the order of the expansion.

```@repl userguide
t*(3t+2.5)
1/(1-t)
t*(t^2-4)/(t+2)
tI = im*t
(1-t)^3.2
(1+t)^t
t^6  # t is of order 5
```

If no valid Taylor expansion can be computed, an error is thrown, for instance
when a derivative is not defined (or simply diverges):

```@repl userguide
1/t
t^3.2
abs(t)
```

Several elementary functions have been implemented; their coefficients
are computed recursively. At the moment of this writing, these functions
are `exp`, `log`, `sqrt`, the trigonometric functions
`sin`, `cos` and `tan`, their inverses, as well as the hyperbolic functions
`sinh`, `cosh` and `tanh` and their inverses;
more will be added in the future. Note that this way of obtaining the
Taylor coefficients is not the *laziest* way, in particular for many independent
variables. Yet, it is quite efficient, especially for the integration of
ordinary differential equations, which is among the applications we have in
mind (see also
[TaylorIntegration.jl](https://github.com/PerezHz/TaylorIntegration.jl)).

```@repl userguide
exp(t)
log(1-t)
sqrt(1 + t)
imag(exp(tI)')
real(exp(Taylor1([0.0,1im],17))) - cos(Taylor1([0.0,1.0],17)) == 0.0
convert(Taylor1{Rational{Int64}}, exp(t))
```

Again, errors are thrown whenever it is necessary.

```@repl userguide
sqrt(t)
log(t)
```

Differentiating and integrating is straightforward for polynomial expansions in
one variable, using [`derivative`](@ref) and [`integrate`](@ref). These
functions return the corresponding [`Taylor1`](@ref) expansions.
The last coefficient of a derivative is set to zero to keep the
same order as the original polynomial; for the integral, an
integration constant may be set by the user (the default is zero). The
order of the resulting polynomial is not changed. The value of the
``n``-th (``n \ge 0``)
derivative is obtained using `derivative(n,a)`, where `a` is a Taylor series.

```@repl userguide
derivative(exp(t))
integrate(exp(t))
integrate( exp(t), 1.0)
integrate( derivative( exp(-t)), 1.0 ) == exp(-t)
derivative(1, exp(shift_taylor(1.0))) == exp(1.0)
derivative(5, exp(shift_taylor(1.0))) == exp(1.0) # 5-th derivative of `exp(1+t)`
```

To evaluate a Taylor series at a given point, Horner's rule is used via the
function `evaluate(a, dt)`. Here, `dt` is the increment from
the point ``t_0`` around which the Taylor expansion of `a` is calculated,
i.e., the series
is evaluated at ``t = t_0 + dt``. Omitting `dt` corresponds to ``dt = 0``;
see [`evaluate`](@ref).

```@repl userguide
evaluate(exp(shift_taylor(1.0))) - e # exp(t) around t0=1 (order 5), evaluated there (dt=0)
evaluate(exp(t), 1) - e # exp(t) around t0=0 (order 5), evaluated at t=1
evaluate(exp( Taylor1(17) ), 1) - e # exp(t) around t0=0, order 17
tBig = Taylor1(BigFloat, 50) # Independent variable with BigFloats, order 50
eBig = evaluate( exp(tBig), one(BigFloat) )
e - eBig
```

Another way to obtain the value of a `Taylor1` polynomial `p` at a given value `x`, is to call `p` as if it were a function, i.e., `p(x)`:

```@repl userguide
t = Taylor1(15)
p = sin(t)
evaluate(p, pi/2) # value of p at pi/2 using `evaluate`
p(pi/2) # value of p at pi/2 by evaluating p as a function
p(pi/2) == evaluate(p, pi/2)
p(0.0)
p() == p(0.0) # p() is a shortcut to obtain the 0-th order coefficient of `p`
```

Note that the syntax `p(x)` is equivalent to `evaluate(p, x)`, whereas `p()` is
equivalent to `evaluate(p)`. For more details about function-like behavior for a
given type in Julia, see the [Function-like objects](https://docs.julialang.org/en/stable/manual/methods/#Function-like-objects-1)
section of the Julia manual.

Useful shortcuts are `taylor_expand` are `update!`. The former returns
the expansion of a function around a given value `t0`. In turn, `update!`
provides an in-place update of a given Taylor polynomial, that is, it shifts
it further by the provided amount.

```@repl userguide
p = taylor_expand( x -> sin(x), pi/2, order=16) # 16-th order expansion of sin(t) around pi/2
update!(p, 0.025) # updates the expansion given by p, by shifting it further by 0.025
p
```


## Many variables

A polynomial in ``N>1`` variables can be represented in (at least) two ways:
As a vector whose coefficients are homogeneous polynomials of fixed degree, or
as a vector whose coefficients are polynomials in ``N-1`` variables. The
current implementation of `TaylorSeries.jl` corresponds to the first option,
though some infrastructure has been built that permits to develop the second
one. An elegant (lazy) implementation of the second representation
was discussed  [here](https://groups.google.com/forum/#!msg/julia-users/AkK_UdST3Ig/sNrtyRJHK0AJ).

The structure [`TaylorN`](@ref) is constructed as a vector of parameterized
homogeneous polynomials
defined by the type [`HomogeneousPolynomial`](@ref), which in turn is a vector of
coefficients of given order (degree). This implementation imposes the user
to specify the (maximum) order considered and the number of independent
variables at the beginning, which can be conveniently done using
[`set_variables`](@ref). A vector of the resulting Taylor variables is returned:

```@repl userguide
x, y = set_variables("x y")
typeof(x)
x.order
x.coeffs
```

As shown, the resulting objects are of `TaylorN{Float64}` type.
There is an optional `order` keyword argument in [`set_variables`](@ref),
used to specify the maximum order of the `TaylorN` polynomials.

```@repl userguide
set_variables("x y", order=10)
```

Similarly, numbered variables are also available by specifying a single
variable name and the optional keyword argument `numvars`:

```@repl userguide
set_variables("α", numvars=3)
```

Alternatively to `set_variables`, [`get_variables`](@ref) can be used if one
doesn't want to change internal dictionaries. `get_variables` returns a vector
of independent variables  of a desired `order` (lesser than `get_order` so
internals doesn't have to change) with the length and variable names defined
by `set_variables` initially.

```@repl userguide
get_variables(order=2) #vector of independent variables of order 2
```

The function [`show_params_TaylorN`](@ref) displays the current values of the
parameters, in an info block.

```@repl userguide
show_params_TaylorN()
```

Internally, changing the parameters (maximum order and number of variables)
redefines the hash-tables that
translate the index of the coefficients of a [`HomogeneousPolynomial`](@ref)
of given order into the corresponding
multi-variable monomials, or the other way around.
Fixing these values from the start is imperative; the initial (default) values
are `order = 6` and `num_vars=2`.

The easiest way to construct a [`TaylorN`](@ref) object is by defining
the independent variables. This can be done using `set_variables` as above,
or through the method [`TaylorN{T<:Number}(::Type{T}, nv::Int)`](@ref)
for the `nv` independent `TaylorN{T}` variable;
the order can be also specified using the optional keyword argument `order`.

```@repl userguide
x, y = set_variables("x y", numvars=2, order=6);
x
TaylorN(1, order=4) # variable 1 of order 4
TaylorN(Int, 2)    # variable 2, type Int, order=get_order()=6
```

Other ways of constructing [`TaylorN`](@ref) polynomials involve
using [`HomogeneousPolynomial`](@ref)
objects directly, which is uncomfortable.

```@repl userguide
set_variables("x", numvars=2);
HomogeneousPolynomial([1,-1])
TaylorN([HomogeneousPolynomial([1,0]), HomogeneousPolynomial([1,2,3])],4)
```

The Taylor expansions are implemented around 0 for all variables; if the
expansion
is needed around a different value, the trick is a simple translation of
the corresponding independent variable, i.e. ``x \to x+a``.

As before, the usual arithmetic operators (`+`, `-`, `*`, `/`, `^`, `==`)
have been extended to work with [`TaylorN`](@ref) objects, including the
appropriate
promotions to deal with numbers. (Some of the arithmetic operations have
been extended for
[`HomogeneousPolynomial`](@ref), whenever the result is a
[`HomogeneousPolynomial`](@ref); division, for instance, is not extended.)

Also, the elementary functions have been
implemented, again by computing their coefficients recursively:

```@repl userguide
x, y = set_variables("x y", order=10);
exy = exp(x+y)
```

The function [`get_coeff`](@ref)
gives the normalized coefficient of the polynomial that corresponds to the
monomial specified by a vector `v` containing the powers. For instance, for
the polynomial `exy` above, the coefficient of the monomial ``x^3 y^5`` is

```@repl userguide
get_coeff(exy, [3,5])
rationalize(ans)
```

Partial differentiation is also implemented for [`TaylorN`](@ref) objects,
through the function [`derivative`](@ref), specifying the number
of the variable as the second argument; integration is yet to be implemented.

```@repl userguide
p = x^3 + 2x^2 * y - 7x + 2
q = y - x^4
derivative( p, 1 )   # partial derivative with respect to 1st variable
derivative( q, 2 )
```

If we ask for the partial derivative with respect to a non-defined variable,
an error is thrown.

```@repl userguide
derivative( q, 3 )   # error, since we are dealing with 2 variables
```

[`evaluate`](@ref) can also be used for [`TaylorN`](@ref) objects, using
it on vectors of
numbers (`Real` or `Complex`); the length of the vector must coincide with the
number of independent variables.

```@repl userguide
evaluate(exy, [.1,.02]) == e^0.12
```

Analogously to `Taylor1`, another way to obtain the value of a `TaylorN` polynomial `p` at a given point `x`, is to call `p` as if it were a function:

```@repl userguide
exy([.1,.02])
exy([.1,.02]) == e^0.12
```

Again, the syntax `p(x)` for `p::TaylorN` is equivalent to `evaluate(p,x)`, and
`p()` is equivalent to `evaluate(p)`.

The functions `taylor_expand` and `update!` work as well for `TaylorN`.

```@repl userguide
xysq = x^2 + y^2
update!(xysq, [1.0, -2.0]) # expand around (1,-2)
xysq
update!(xysq, [-1.0, 2.0]) # shift-back
xysq == x^2 + y^2
```

Functions to compute the gradient, Jacobian and
Hessian have also been implemented. Using the
polynomials ``p = x^3 + 2x^2 y - 7 x + 2`` and ``q = y-x^4`` defined above,
we may use [`gradient`](@ref) (or `∇`); the results are of
type `Array{TaylorN{T},1}`. To compute the Jacobian and Hessian of a vector field
evaluated at a point, we use respectively [`jacobian`](@ref) and
[`hessian`](@ref):

```@repl userguide
∇(p)
gradient( q )
r = p-q-2*p*q
hessian(ans)
jacobian([p,q], [2,1])
hessian(r, [1.0,1.0])
```

Other specific applications are described in the
[Examples](@ref).

## Mixtures

As mentioned above, `Taylor1{T}`, `HomogeneousPolynomial{T}` and `TaylorN{T}`
are parameterized structures such that `T<:AbstractSeries`, the latter
is a subtype of `Number`. Then, we may actually define Taylor expansions in
``N+1`` variables, where one of the variables (the `Taylor1` variable) is
somewhat special.

```@repl userguide
x, y = set_variables("x y", order=3)
t1N = Taylor1([zero(x), one(x)], 5)
```

The last line defines a `Taylor1{TaylorN{Float64}}` variable, which is of order
5 in `t` and order 3 in `x` and `y`. Then, we can evaluate functions involving
such polynomials:

```@repl userguide
cos(2.1+x+t1N)
```

This kind of expansions are of interest when studying the dependence of
parameters, for instance in the context of bifurcation theory or when considering
the dependence of the solution of a differential equation on the initial conditions,
around a given solution. In this case, `x` and `y` represent small variations
around a given value of the parameters, or around some specific initial condition.
Such constructions are exploited in the package [`TaylorIntegration.jl`](https://github.com/PerezHz/TaylorIntegration.jl).

```@meta
CurrentModule = nothing
```
