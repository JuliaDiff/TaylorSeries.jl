
    {meta}
    CurrentModule = TaylorSeries

# User guide

---

[TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl)
is a basic polynomial algebraic manipulator in one or more
variables; these two cases are treated separately.  Three new types are defined,
[`Taylor1`]({ref}), [`HomogeneousPolynomial`]({ref}) and [`TaylorN`]({ref}),
which correspond to
expansions in one independent variable, homogeneous polynomials of various variables,
and the polynomial
series in many independent variables, respectively. These types are subtypes
of `Number` and are defined parametrically.

The package is loaded as usual:

    {repl userguide}
    using TaylorSeries


## One variable

Taylor expansions in one variable are represented by the [`Taylor1`]({ref}) type, which
consists of a vector of coefficients (field `coeffs`) and the maximum
order considered for the expansion (field `order`). The
coefficients are arranged in ascending order with respect to the power of the
independent variable, so that
`coeffs[1]` is the constant term, `coeffs[2]` gives the first order term,
etc. This is a dense representation of the polynomial.
The order of the polynomial can be
omitted in the constructor, which is then fixed from the length of the
vector of coefficients; otherwise, the maximum
of the length of the vector of coefficients and the given integer is taken.

    {repl userguide}
    Taylor1([1, 2, 3]) # Polynomial of order 2 with coefficients 1, 2, 3
    Taylor1([0.0, 1im]) # Also works with complex numbers
    affine(a) = a + taylor1_variable(typeof(a),5)  ## a + t of order 5
    t = affine(0.0) # Independent variable `t`


Note that the information about the maximum order considered is displayed
using a big-O notation.

The definition of `affine(a)` uses the function [`taylor1_variable()`]({ref}), which is a
shortcut to define the independent variable of a Taylor expansion,
with a given type and given order. As we show below, this is one of the
easiest ways to work with the package.

The usual arithmetic operators (`+`, `-`, `*`, `/`, `^`, `==`) have been
extended to work with the [`Taylor1`]({ref}) type, including promotions that involve
`Number`s. The operations return a valid Taylor expansion with the same
maximum order; compare the last example below, where this is not possible:

    {repl userguide}
    t*(3t+2.5)
    1/(1-t)
    t*(t^2-4)/(t+2)
    tI = im*t
    t^6  # order is 5
    (1-t)^3.2
    (1+t)^t

If no valid Taylor expansion can be computed, an error is thrown.

    {repl userguide}
    1/t
    t^3.2

Several elementary functions have been implemented; these compute their
coefficients recursively. So far, these functions are `exp`, `log`,
`sqrt`, `sin`, `cos` and `tan`;
more will be added in the future. Note that this way of obtaining the
Taylor coefficients is not the *laziest* way, in particular for many independent
variables. Yet, it is quite efficient, especially for the integration of
ordinary differential equations, which is among the applications we have in mind.

    {repl userguide}
    exp(t)
    log(1-t)
    sqrt(1 + t)
    imag(exp(tI)')
    real(exp(Taylor1([0.0,1im],17))) - cos(Taylor1([0.0,1.0],17)) == 0.0
    convert(Taylor1{Rational{Int64}}, exp(t))

Again, errors are thrown whenever it is necessary.

    {repl userguide}
    sqrt(t)
    log(t)


Differentiating and integrating is straightforward for polynomial expansions in
one variable, using [`diffTaylor()`]({ref}) and [`integTaylor()`]({ref}). These
functions return the corresponding [`Taylor1`]({ref}) expansions.
The last coefficient of a derivative is set to zero to keep the
same order as the original polynomial; for the integral, an
integration constant may be set (the default is zero). The
order of the resulting polynomial is not changed. The value of the
$n$-th ($n \ge 0$)
derivative is obtained using `deriv(a,n)`, where `a` is a Taylor series;
the default is $n=1$; see [`deriv()`]({ref}).

    {repl userguide}
    diffTaylor(exp(t))
    integTaylor(exp(t))
    integTaylor( exp(t), 1.0)
    integTaylor( diffTaylor( exp(-t)), 1.0 ) == exp(-t)
    deriv( exp(affine(1.0))) == exp(1.0)
    deriv( exp(affine(1.0)), 5) == exp(1.0) # Fifth derivative of `exp(1+t)`


To evaluate a Taylor series at a point, Horner's rule is used via the function
`evaluate(a::Taylor1, dt::Number)`. Here, `dt` is the increment from
the point $t_0$ where the Taylor expansion of `a` is calculated, i.e., the series
is evaluated at $t = t_0 + dt$. Omitting `dt` corresponds to $dt = 0$.
See [`evaluate()`]({ref}).

    {repl userguide}
    evaluate(exp(affine(1.0))) - e # exp(t) around t0=1 (order 5), evaluated there (dt=0)
    evaluate(exp(t), 1) - e # exp(t) around t0=0 (order 5), evaluated at t=1
    evaluate( exp( taylor1_variable(17) ), 1) - e # exp(t) around t0=0 (order 17),
    tBig = Taylor1([zero(BigFloat),one(BigFloat)],50) # With BigFloats
    eBig = evaluate( exp(tBig), one(BigFloat) )
    e - eBig


## Many variables

A polynomial in $N>1$ variables can be represented in (at least) two ways:
As a vector whose coefficients are homogeneous polynomials of fixed degree, or
as a vector whose coefficients are polynomials in $N-1$ variables. We have opted
to implement the first option, which seems to show better performance. An elegant
(lazy) implementation of the second representation was discussed on the
[julia-users](https://groups.google.com/forum/#!msg/julia-users/AkK_UdST3Ig/sNrtyRJHK0AJ) list.

[`TaylorN`]({ref}) is thus constructed as a vector of parameterized homogeneous polynomials
defined by the type [`HomogeneousPolynomial`]({ref}), which in turn is a vector of
coefficients of given order (degree). This implementation imposes that the user
has to specify the (maximum) order and the number of independent
variables, which is done using the [`set_variables()`]({ref}) function.
A vector of the resulting Taylor variables is returned:

    {repl userguide}
    x, y = set_variables("x y")
    typeof(x)
    x.order
    x.coeffs

As shown, the resulting objects are of [`TaylorN{Float64}`]({ref}) type.
There is an optional `order` keyword argument for [`set_variables()`]({ref}):

    {repl userguide}
    set_variables("x y", order=10)
    x # x variable

Numbered variables are also available by specifying a single
variable name and the optional keyword argument `numvars`:

    {repl userguide}
    set_variables("α", numvars=3)

The function [`show_params_TaylorN()`]({ref}) displays the current values of the
parameters, in an info block.

    julia> show_params_TaylorN()
    INFO: Parameters for `TaylorN` and `HomogeneousPolynomial`:
    Maximum order       = 6
    Number of variables = 3
    Variable names      = UTF8String["α₁","α₂","α₃"]

Internally, changing these parameters defines dictionaries that
translate the position of the coefficients of a `HomogeneousPolynomial`
into the corresponding
multi-variable monomials. Fixing these values from the start is imperative.

The easiest way to construct a [`TaylorN`]({ref}) object is by defining symbols for
the independent variables, as above. Again, the Taylor expansions are implemented
around 0 for all variables; if the expansion
is needed around a different value, the trick is a simple translation of
the corresponding
independent variable $x \to x+a$.

Other ways of constructing [`TaylorN`]({ref}) polynomials involve
using [`HomogeneousPolynomial`]({ref})
objects directly, which is uncomfortable:

    {repl userguide}
    set_variables("x", numvars=2);
    HomogeneousPolynomial([1,-1])
    TaylorN([HomogeneousPolynomial([1,0]), HomogeneousPolynomial([1,2,3])],4)

As before, the usual arithmetic operators (`+`, `-`, `*`, `/`, `^`, `==`)
have been extended to work with [`TaylorN`]({ref}) objects, including the appropriate
promotions to deal with numbers. (Some of the arithmetic operations have
also been extended for
[`HomogeneousPolynomial`]({ref}), whenever the result is a
[`HomogeneousPolynomial`]({ref}); division, for instance, is not extended.)
Also, the elementary functions have been
implemented, again by computing their coefficients recursively:


    {repl userguide}
    x, y = set_variables("x y", order=10);
    exy = exp(x+y)

The function `get_coeff(a,v)` (see [`get_coeff()`]({ref}))
gives the coefficient of `x` that corresponds to the monomial
specified by the vector of powers `v`:

    {repl userguide}
    get_coeff(exy, [3,5])
    rationalize(ans)

Partial differentiation is also implemented for [`TaylorN`]({ref}) objects,
through [`diffTaylor()`]({ref}); integration is yet to be implemented.

    {repl userguide}
    f(x,y) = x^3 + 2x^2 * y - 7x + 2
    g(x,y) = y - x^4
    diffTaylor( f(x,y), 1 )   # partial derivative with respect to 1st variable
    diffTaylor( g(x,y), 2 )

If we ask for the partial derivative with respect to a non-defined variable,
an error is thrown.

    {repl userguide}
    diffTaylor( g(x,y), 3 )   # error, since we are dealing with 2 variables

[`evaluate()`]({ref}) can also be used for [`TaylorN`]({ref}) objects, using
it on vectors of
numbers (`Real` or `Complex`); the length of the vector must coincide with the number
of independent variables.

    {repl userguide}
    evaluate(exy, [.1,.02]) == e^0.12

Functions to compute the gradient, Jacobian and
Hessian have also been implemented. Using the
functions $f(x,y) = x^3 + 2x^2 y - 7 x + 2$ and $g(x,y) = y-x^4$ defined above,
we may use [`gradient()`]({ref}) or `∇` (`\nabla+TAB`); the results are of
type `Array{TaylorN{T},1}`. To compute the Jacobian and Hessian of a vector field
evaluated at a point, we use respectively [`jacobian()`]({ref}) and
[`hessian()`]({ref}):

    {repl userguide}
    f1 = f(x,y)
    g1 = g(x,y)
    ∇(f1)
    gradient( g1 )
    fg = f1-g1-2*f1*g1
    hessian(ans)
    jacobian([f1,g1], [2,1])
    hessian(fg, [1.0,1.0])

Some specific applications are given in the next [section](Examples).

