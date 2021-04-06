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
etc. Yet, it is possible to have the natural ordering with respect
to the degree; see below. This is a dense representation of the polynomial.
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
using a big-𝒪 notation. The convention followed when different orders are
combined is consistent with the mathematics and the big-𝒪 notation, that is,
to propagate the lowest order. In some cases, it is desirable to not display
the big-𝒪 notation. The function [`displayBigO`](@ref) allows to
control whether it is displayed or not.
```@repl userguide
displayBigO(false) # turn-off displaying big O notation
t
displayBigO(true) # turn it on
t
```

Similarly, it is possible to control if the format of the
displayed series through the function [`use_show_default`](@ref);
`use_show_default(true)` uses the `Base.show_default`, while
`use_show_default(false)` uses the custom display form (default).
```@repl userguide
use_show_default(true) # use Base.show method
t
use_show_default(false) # use custom `show`
t
```

The definition of `shift_taylor(a)` uses the method
[`Taylor1([::Type{Float64}], order::Int)`](@ref), which is a
shortcut to define the independent variable of a Taylor expansion,
of given type and order (the default is `Float64`).
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
Taylor1(3) + Taylor1(5) == 2Taylor1(3)  # big-𝒪 convention applies
t^6  # t is of order 5
```

If no valid Taylor expansion can be computed an error is thrown, for instance,
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
more functions will be added in the future. Note that this way of obtaining the
Taylor coefficients is not a *lazy* way, in particular for many independent
variables. Yet, it is quite efficient, especially for the integration of
ordinary differential equations, which is among the applications we have in
mind (see
[TaylorIntegration.jl](https://github.com/PerezHz/TaylorIntegration.jl)).

```@repl userguide
exp(t)
log(1-t)
sqrt(1 + t)
imag(exp(tI)')
real(exp(Taylor1([0.0,1im],17))) - cos(Taylor1([0.0,1.0],17)) == 0.0
convert(Taylor1{Rational{Int}}, exp(t))
```

Again, errors are thrown whenever it is necessary.

```@repl userguide
sqrt(t)
log(t)
```

To obtain a specific coefficient, [`getcoeff`](@ref) can be used. Another
alternative is to request the specific degree using the vector notation,
where the index corresponds to the degree of the term.

```@repl userguide
expon = exp(t)
getcoeff(expon, 0) == expon[0]
rationalize(expon[3])
```

Differentiating and integrating is straightforward for polynomial expansions in
one variable, using [`differentiate`](@ref) and [`integrate`](@ref). (The
function [`derivative`](@ref) is a synonym of `differentiate`.) These
functions return the corresponding [`Taylor1`](@ref) expansions. Note that
the order of the derivative of a `Taylor1` corresponds to the order of the
original polynomial *minus 1*. For the integral, an integration constant may be
set by the user (the default is zero); the order of the integrated polynomial
for the integral is *kept unchanged*. The *value* of the ``n``-th (``n \ge 0``)
derivative is obtained using `differentiate(n,a)`, where `a` is a Taylor series;
likewise, the `Taylor1` polynomial of the ``n``-th derivative is obtained as
`differentiate(a,n)`; the resulting polynomial is of order `get_order(a)-n`.

```@repl userguide
differentiate(exp(t)) # exp(t) is of order 5; the derivative is of order 4
integrate(exp(t))  # the resulting TaylorSeries is of order 5
integrate( exp(t), 1.0)
integrate( differentiate( exp(-t)), 1.0 ) == exp(-t)
differentiate(1, exp(shift_taylor(1.0))) == exp(1.0)
differentiate(5, exp(shift_taylor(1.0))) == exp(1.0)    # 5-th differentiate of `exp(1+t)`
derivative(exp(1+t), 3)    # Taylor1 polynomial of the 3-rd derivative of `exp(1+t)`
```

To evaluate a Taylor series at a given point, Horner's rule is used via the
function `evaluate(a, dt)`. Here, `dt` is the increment from
the point ``t_0`` around which the Taylor expansion of `a` is calculated,
i.e., the series
is evaluated at ``t = t_0 + dt``. Omitting `dt` corresponds to ``dt = 0``;
see [`evaluate`](@ref).

```@repl userguide
evaluate(exp(shift_taylor(1.0))) - ℯ    # exp(t) around t0=1 (order 5), evaluated there (dt=0)
evaluate(exp(t), 1) - ℯ                 # exp(t) around t0=0 (order 5), evaluated at t=1
evaluate(exp( Taylor1(17) ), 1) - ℯ     # exp(t) around t0=0, order 17
tBig = Taylor1(BigFloat, 50)            # Independent variable with BigFloats, order 50
eBig = evaluate( exp(tBig), one(BigFloat) )
ℯ - eBig
```

Another way to obtain the value of a `Taylor1` polynomial `p` at a given value `x`, is to call `p` as if it was a function, i.e., `p(x)`:

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
equivalent to `evaluate(p)`.

Useful shortcuts are [`taylor_expand`](@ref) and [`update!`](@ref).
The former returns
the expansion of a function around a given value `t0`, mimicking the use
of `shift_taylor` above. In turn, `update!`
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
used to specify the maximum order of the `TaylorN` polynomials. Note that
one can specify the variables using a vector of symbols.

```@repl userguide
set_variables([:x, :y], order=10)
```

Similarly, subindexed variables are also available by specifying a single
variable name and the optional keyword argument `numvars`:

```@repl userguide
set_variables("α", numvars=3)
```

Alternatively to `set_variables`, [`get_variables`](@ref) can be used if one
does not want to change internal dictionaries. `get_variables` returns a vector
of `TaylorN` independent variables of a desired `order`
(lesser than `get_order` so the
internals doesn't have to change) with the length and variable names defined
by `set_variables` initially.

```@repl userguide
get_variables(2) # vector of independent variables of order 2
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
set_variables(:x, numvars=2); # symbols can be used
HomogeneousPolynomial([1,-1])
TaylorN([HomogeneousPolynomial([1,0]), HomogeneousPolynomial([1,2,3])],4)
```

The Taylor expansions are implemented around 0 for all variables; if the
expansion
is needed around a different value, the trick is a simple translation of
the corresponding independent variable, i.e. ``x \to x+a``.

As before, the usual arithmetic operators (`+`, `-`, `*`, `/`, `^`, `==`)
have been extended to work with [`TaylorN`](@ref) objects, including the
appropriate promotions to deal with numbers.
Note that some of the arithmetic operations have been extended for
[`HomogeneousPolynomial`](@ref), whenever the result is a
[`HomogeneousPolynomial`](@ref); division, for instance, is not extended.
The same convention used for `Taylor1` objects is used when combining
`TaylorN` polynomials of different order.

The elementary functions have also been
implemented, again by computing their coefficients recursively:

```@repl userguide
x, y = set_variables("x y", order=10);
exy = exp(x+y)
```

The function [`getcoeff`](@ref)
gives the normalized coefficient of the polynomial that corresponds to the
monomial specified by the tuple or vector `v` containing the powers.
For instance, for
the polynomial `exy` above, the coefficient of the monomial ``x^3 y^5`` is
obtained using `getcoeff(exy, (3,5))` or `getcoeff(exy, [3,5])`.
```@repl userguide
getcoeff(exy, (3,5))
rationalize(ans)
```

Similar to `Taylor1`, vector notation can be used to request specific
coefficients of `HomogeneousPolynomial` or `TaylorN` objects. For `TaylorN`
objects, the index refers to the degree of the `HomogeneousPolynomial`.
In the case of `HomogeneousPolynomial` the index refers to the position
of the hash table. The function [`show_monomials`](@ref) can be used to
obtain the coefficient a specific monomial, given the degree of the
`HomogeneousPolynomial`.

```@repl userguide
exy[8] # get the 8th order term
show_monomials(8)
exy[8][6] # get the 6th coeff of the 8th order term
```

Partial differentiation is also implemented for [`TaylorN`](@ref) objects,
through the function [`differentiate`](@ref), specifying the number
of the variable, or its symbol, as the second argument.

```@repl userguide
p = x^3 + 2x^2 * y - 7x + 2
q = y - x^4
differentiate( p, 1 )   # partial derivative with respect to 1st variable
differentiate( q, :y )  # partial derivative with respect to :y
```

If we ask for the partial derivative with respect to a non-defined variable,
an error is thrown.

```@repl userguide
differentiate( q, 3 )   # error, since we are dealing with 2 variables
```

To obtain more specific partial derivatives we have two specialized methods
that involve a tuple, which represents the number of derivatives with
respect to each variable (so the tuple's length has to be the
same as the actual number of variables). These methods either return
the `TaylorN` object in question, or the coefficient corresponding to
the specified tuple, normalized by the factorials defined by the tuple.
The latter is in essence the 0-th order coefficient of the former.

```@repl userguide
differentiate(p, (2,1)) # two derivatives on :x and one on :y
differentiate((2,1), p) # 0-th order coefficient of the previous expression
differentiate(p, (1,1)) # one derivative on :x and one on :y
differentiate((1,1), p) # 0-th order coefficient of the previous expression
```

Integration with respect to the `r`-th variable for
`HomogeneousPolynomial`s and `TaylorN` objects is obtained
using [`integrate`](@ref). Note that `integrate` for `TaylorN`
objects allows to specify a constant of integration, which must
be independent from the integrated variable. Again, the integration
variable may be specified by its symbol.

```@repl userguide
integrate( differentiate( p, 1 ), 1) # integrate with respect to the first variable
integrate( differentiate( p, 1 ), :x, 2) # integration with respect to :x, constant of integration is 2
integrate( differentiate( q, 2 ), :y, -x^4) == q
integrate( differentiate( q, 2 ), 2, y)
```

[`evaluate`](@ref) can also be used for [`TaylorN`](@ref) objects, using
it on vectors of
numbers (`Real` or `Complex`); the length of the vector must coincide with the
number of independent variables. [`evaluate`](@ref) also allows to specify only
one variable and a value.

```@repl userguide
evaluate(exy, [.1,.02]) == exp(0.12)
evaluate(exy, :x, 0.0) == exp(y)  # evaluate `exy` for :x -> 0
```

Analogously to `Taylor1`, another way to obtain the value of a `TaylorN`
polynomial `p` at a given point `x`, is to call it as if it were a function:
the syntax `p(x)` for `p::TaylorN` is equivalent to `evaluate(p,x)`, and
`p()` is equivalent to `evaluate(p)`.

```@repl userguide
exy([.1,.02]) == exp(0.12)
exy(:x, 0.0)
```

Internally, `evaluate` for `TaylorN` considers separately
the contributions of all `HomogeneousPolynomial`s by `order`,
which are finally added up *after* sorting them in place (which is the default)
in increasing order by `abs2`. This is done in order to
use as many significant figures as possible of all terms
in the final sum, which then should yield a more
accurate result. This default can be changed to a non-sorting
sum thought, which may be more performant or useful for
certain subtypes of `Number` which, for instance, do not have `isless`
defined. See
[this issue](https://github.com/JuliaDiff/TaylorSeries.jl/issues/242)
for a motivating example. This can be done using the keyword
`sorting` in `evaluate`, which expects a `Bool`, or using a
that boolean as the *first* argument in the function-like evaluation.

```@repl userguide
exy([.1,.02]) # default is `sorting=true`
evaluate(exy, [.1,.02]; sorting=false)
exy(false, [.1,.02])
```

In the examples shown above, the first entry corresponds to the
default case (`sorting=true`), which yields the same result as
`exp(0.12)`, and the remaining two illustrate
turning off sorting the terms. Note that the results are not
identical, since [floating point addition is not
associative](https://en.wikipedia.org/wiki/Associative_property#Nonassociativity_of_floating_point_calculation),
which may introduce rounding errors.

The functions `taylor_expand` and `update!` work as well for `TaylorN`.

```@repl userguide
xysq = x^2 + y^2
update!(xysq, [1.0, -2.0]) # expand around (1,-2)
xysq
update!(xysq, [-1.0, 2.0]) # shift-back
xysq == x^2 + y^2
```

Functions to compute the gradient, Jacobian and
Hessian have also been implemented; note that these
functions *are not* exported, so its use require the
prefix `TaylorSeries`. Using the
polynomials ``p = x^3 + 2x^2 y - 7 x + 2`` and ``q = y-x^4`` defined above,
we may use [`TaylorSeries.gradient`](@ref) (or `∇`); the results are of
type `Array{TaylorN{T},1}`. To compute the Jacobian and Hessian of a vector field
evaluated at a point, we use respectively [`TaylorSeries.jacobian`](@ref) and
[`TaylorSeries.hessian`](@ref):

```@repl userguide
∇(p)
TaylorSeries.gradient( q )
r = p-q-2*p*q
TaylorSeries.hessian(ans)
TaylorSeries.jacobian([p,q], [2,1])
TaylorSeries.hessian(r, [1.0,1.0])
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
