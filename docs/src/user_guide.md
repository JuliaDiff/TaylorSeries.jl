# User guide

---

    {meta}
    CurrentModule = TaylorSeries
    DocTestSetup = quote
        using TaylorSeries
        affine(a) = a + taylor1_variable(typeof(a),5)
        t = taylor1_variable(5)
        tI = im * t
        x, y = set_variables("x y", order=10);
        exy = exp(x+y)
        f(x,y) = x^3 + 2x^2 * y - 7x + 2
        g(x,y) = y - x^4
    end

`TaylorSeries.jl` can be thought of as a polynomial algebraic manipulator in one or more
variables; these two cases are treated separately.  Three new types are defined,
[`Taylor1`]({ref}), [`HomogeneousPolynomial`]({ref}) and [`TaylorN`]({ref}),
which correspond to
expansions in one independent variable, homogeneous polynomials of various variables,
and the polynomial
series in many independent variables, respectively. These types are subtypes
of `Number` and are defined parametrically.

The package is loaded as usual:

```julia
julia> using TaylorSeries

```

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

```julia
julia> Taylor1([1, 2, 3]) # Polynomial of order 2 with coefficients 1, 2, 3
 1 + 2 t + 3 t² + 𝒪(t³)

julia> Taylor1([0.0, 1im]) # Also works with complex numbers
 ( 1.0 im ) t + 𝒪(t²)

julia> affine(a) = a + taylor1_variable(typeof(a),5)  ## a + t of order 5
affine (generic function with 1 method)

julia> t = affine(0.0) # Independent variable `t`
 1.0 t + 𝒪(t⁶)

```

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

```julia
julia> t*(3t+2.5)
 2.5 t + 3.0 t² + 𝒪(t⁶)

julia> 1/(1-t)
 1.0 + 1.0 t + 1.0 t² + 1.0 t³ + 1.0 t⁴ + 1.0 t⁵ + 𝒪(t⁶)

julia> t*(t^2-4)/(t+2)
 - 2.0 t + 1.0 t² + 𝒪(t⁶)

julia> tI = im*t
 ( 1.0 im ) t + 𝒪(t⁶)

julia> t^6  # order is 5
 0.0 + 𝒪(t⁶)

julia> (1-t)^3.2
 1.0 - 3.2 t + 3.5200000000000005 t² - 1.4080000000000004 t³ + 0.07040000000000009 t⁴ + 0.011264000000000012 t⁵ + 𝒪(t⁶)

julia> (1+t)^t
 1.0 + 1.0 t² - 0.5 t³ + 0.8333333333333333 t⁴ - 0.75 t⁵ + 𝒪(t⁶)

```

If no valid Taylor expansion can be computed, an error is thrown.

```julia
julia> 1/t
ERROR: ArgumentError: Division does not define a Taylor1 polynomial
or its first non-zero coefficient is Inf/NaN.
Order k=0 => coeff[1]=Inf.
 in divfactorization at /Users/benet/.julia/v0.4/TaylorSeries/src/Taylor1.jl:239
 in / at /Users/benet/.julia/v0.4/TaylorSeries/src/Taylor1.jl:215
 in / at ./promotion.jl:170

julia> t^3.2
ERROR: ArgumentError: The 0th order Taylor1 coefficient must be non-zero
to raise the Taylor1 polynomial to a non-integer exponent.
 in ^ at /Users/benet/.julia/v0.4/TaylorSeries/src/Taylor1.jl:355

```

Several elementary functions have been implemented; these compute their
coefficients recursively. So far, these functions are `exp`, `log`, `sqrt`, `sin`, `cos`
and `tan`;
more will be added in the future. Note that this way of obtaining the
Taylor coefficients is not the *laziest* way, in particular for many independent
variables. Yet, it is quite efficient, especially for the integration of
ordinary differential equations, which is among the applications we have in mind.

```julia
julia> exp(t)
 1.0 + 1.0 t + 0.5 t² + 0.16666666666666666 t³ + 0.041666666666666664 t⁴ + 0.008333333333333333 t⁵ + 𝒪(t⁶)

julia> log(1-t)
 - 1.0 t - 0.5 t² - 0.3333333333333333 t³ - 0.25 t⁴ - 0.2 t⁵ + 𝒪(t⁶)

julia> sqrt(t)
ERROR: ArgumentError: First non-vanishing Taylor1 coefficient must correspond
to an **even power** in order to expand `sqrt` around 0.
 in sqrt at /Users/benet/.julia/v0.4/TaylorSeries/src/Taylor1.jl:427

julia> sqrt(1 + t)
 1.0 + 0.5 t - 0.125 t² + 0.0625 t³ - 0.0390625 t⁴ + 0.02734375 t⁵ + 𝒪(t⁶)

julia> imag(exp(tI)')
 - 1.0 t + 0.16666666666666666 t³ - 0.008333333333333333 t⁵ + 𝒪(t⁶)

julia> real(exp(Taylor1([0.0,1im],17))) - cos(Taylor1([0.0,1.0],17)) == 0.0
true

julia> convert(Taylor1{Rational{Int64}}, exp(t))
 1//1 + 1//1 t + 1//2 t² + 1//6 t³ + 1//24 t⁴ + 1//120 t⁵ + 𝒪(t⁶)

```

Differentiating and integrating is straightforward for polynomial expansions in
one variable, using `diffTaylor()` and `integTaylor()`. These
functions return the corresponding [`Taylor1`]({ref}) expansions.
The last coefficient of a derivative is set to zero to keep the
same order as the original polynomial; for the integral, an
integration constant may be set (the default is zero). The
order of the resulting polynomial is not changed. The value of the
$n$-th ($n \ge 0$)
derivative is obtained using `deriv(a,n)`, where `a` is a Taylor series;
the default is $n=1$.

```julia
julia> diffTaylor(exp(t))
 1.0 + 1.0 t + 0.5 t² + 0.16666666666666666 t³ + 0.041666666666666664 t⁴ + 𝒪(t⁶)

julia> integTaylor(exp(t))
 1.0 t + 0.5 t² + 0.16666666666666666 t³ + 0.041666666666666664 t⁴ + 0.008333333333333333 t⁵ + 𝒪(t⁶)

julia> integTaylor( exp(t), 1.0)
 1.0 + 1.0 t + 0.5 t² + 0.16666666666666666 t³ + 0.041666666666666664 t⁴ + 0.008333333333333333 t⁵ + 𝒪(t⁶)

julia> integTaylor( diffTaylor( exp(-t)), 1.0 ) == exp(-t)
true

julia> deriv( exp(affine(1.0))) == exp(1.0)
true

julia> deriv( exp(affine(1.0)), 5) == exp(1.0) # Fifth derivative of `exp(1+t)`
true

```

To evaluate a Taylor series at a point, Horner's rule is used via the function
`evaluate(a::Taylor1, dt::Number)`. Here, `dt` is the increment from
the point $t_0$ where the Taylor expansion of `a` is calculated, i.e., the series
is evaluated at $t = t_0 + dt$. Omitting `dt` corresponds to $dt = 0$.
See [`evaluate()`]({ref}).

```julia
julia> evaluate(exp(affine(1.0))) - e # exp(t) around t0=1 (order 5), evaluated there (dt=0)

0.0

julia> evaluate(exp(t), 1) - e # exp(t) around t0=0 (order 5), evaluated at t=1
-0.0016151617923783057

julia> evaluate( exp( taylor1_variable(17) ), 1) - e # exp(t) around t0=0 (order 17),
0.0

julia> tBig = Taylor1([zero(BigFloat),one(BigFloat)],50) # With BigFloats
 1.000000000000000000000000000000000000000000000000000000000000000000000000000000 t + 𝒪(t⁵¹)

julia> eBig = evaluate( exp(tBig), one(BigFloat) )
2.718281828459045235360287471352662497757247093699959574966967627723419298053556

julia> e - eBig
6.573322999985292556154129119543257102601105719980995128942636339920549561322098e-67

```


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

```julia
julia> x, y = set_variables("x y")
2-element Array{TaylorSeries.TaylorN{Float64},1}:
  1.0 x + 𝒪(‖x‖⁷)
  1.0 y + 𝒪(‖x‖⁷)

julia> x
 1.0 x + 𝒪(‖x‖⁷)

julia> typeof(x)
TaylorSeries.TaylorN{Float64}

julia> x.order
6

julia> x.coeffs
7-element Array{TaylorSeries.HomogeneousPolynomial{Float64},1}:
    0.0
  1.0 x
    0.0
    0.0
    0.0
    0.0
    0.0

```

As shown, the resulting objects are of [`TaylorN{Float64}`]({ref}) type.
There is an optional `order` keyword argument for [`set_variables()`]({ref}):

```julia
julia> set_variables("x y", order=10)
2-element Array{TaylorSeries.TaylorN{Float64},1}:
  1.0 x + 𝒪(‖x‖¹¹)
  1.0 y + 𝒪(‖x‖¹¹)

```

Numbered variables are also available by specifying a single
variable name and the optional keyword argument `numvars`:

```julia
julia> set_variables("α", numvars=3)
3-element Array{TaylorSeries.TaylorN{Float64},1}:
  1.0 α₁ + 𝒪(‖x‖⁷)
  1.0 α₂ + 𝒪(‖x‖⁷)
  1.0 α₃ + 𝒪(‖x‖⁷)

```

The function [`show_params_TaylorN()`({ref})] displays the current values of the
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

```julia
julia> set_variables("x", numvars=2);

julia> HomogeneousPolynomial([1,-1])
 1 x₁ - 1 x₂

julia> TaylorN([HomogeneousPolynomial([1,0]), HomogeneousPolynomial([1,2,3])],4)
 1 x₁ + 1 x₁² + 2 x₁ x₂ + 3 x₂² + 𝒪(‖x‖⁵)

```

As before, the usual arithmetic operators (`+`, `-`, `*`, `/`, `^`, `==`)
have been extended to work with [`TaylorN`]({ref}) objects, including the appropriate
promotions to deal with numbers. (Some of the arithmetic operations have
also been extended for
[`HomogeneousPolynomial`]({ref}), whenever the result is a
[`HomogeneousPolynomial`]({ref}); division, for instance, is not extended.)
Also, the elementary functions have been
implemented, again by computing their coefficients recursively:

```julia
julia> x, y = set_variables("x y", order=10);

julia> exy = exp(x+y)
 1.0 + 1.0 x + 1.0 y + 0.5 x² + 1.0 x y + 0.5 y² + 0.16666666666666666 x³ + 0.5 x² y + 0.5 x y² + 0.16666666666666666 y³ + 0.041666666666666664 x⁴ + 0.16666666666666666 x³ y + 0.25 x² y² + 0.16666666666666666 x y³ + 0.041666666666666664 y⁴ + 0.008333333333333333 x⁵ + 0.041666666666666664 x⁴ y + 0.08333333333333333 x³ y² + 0.08333333333333333 x² y³ + 0.041666666666666664 x y⁴ + 0.008333333333333333 y⁵ + 0.0013888888888888887 x⁶ + 0.008333333333333331 x⁵ y + 0.020833333333333332 x⁴ y² + 0.027777777777777776 x³ y³ + 0.020833333333333332 x² y⁴ + 0.008333333333333331 x y⁵ + 0.0013888888888888887 y⁶ + 0.00019841269841269839 x⁷ + 0.0013888888888888885 x⁶ y + 0.004166666666666666 x⁵ y² + 0.006944444444444443 x⁴ y³ + 0.006944444444444443 x³ y⁴ + 0.004166666666666666 x² y⁵ + 0.0013888888888888885 x y⁶ + 0.00019841269841269839 y⁷ + 2.4801587301587298e-5 x⁸ + 0.00019841269841269836 x⁷ y + 0.0006944444444444443 x⁶ y² + 0.0013888888888888887 x⁵ y³ + 0.0017361111111111108 x⁴ y⁴ + 0.0013888888888888887 x³ y⁵ + 0.0006944444444444443 x² y⁶ + 0.00019841269841269836 x y⁷ + 2.4801587301587298e-5 y⁸ + 2.7557319223985884e-6 x⁹ + 2.4801587301587295e-5 x⁸ y + 9.920634920634918e-5 x⁷ y² + 0.0002314814814814814 x⁶ y³ + 0.0003472222222222221 x⁵ y⁴ + 0.0003472222222222221 x⁴ y⁵ + 0.0002314814814814814 x³ y⁶ + 9.920634920634918e-5 x² y⁷ + 2.4801587301587295e-5 x y⁸ + 2.7557319223985884e-6 y⁹ + 2.7557319223985883e-7 x¹⁰ + 2.7557319223985884e-6 x⁹ y + 1.2400793650793647e-5 x⁸ y² + 3.306878306878306e-5 x⁷ y³ + 5.787037037037036e-5 x⁶ y⁴ + 6.944444444444443e-5 x⁵ y⁵ + 5.787037037037036e-5 x⁴ y⁶ + 3.306878306878306e-5 x³ y⁷ + 1.2400793650793647e-5 x² y⁸ + 2.7557319223985884e-6 x y⁹ + 2.7557319223985883e-7 y¹⁰ + 𝒪(‖x‖¹¹)

```

The function [`get_coeff(a,v)`]({ref})
gives the coefficient of `x` that corresponds to the monomial
specified by the vector of powers `v`:

```julia
julia> get_coeff(exy, [3,5]) == 1/720
false

julia> rationalize(get_coeff(exy, [3,5]))
1//720

```

Partial differentiation is also implemented for [`TaylorN`]({ref}) objects,
through [`diffTaylor`]({ref}); integration is yet to be implemented.

```julia
julia> f(x,y) = x^3 + 2x^2 * y - 7x + 2
f (generic function with 1 method)

julia> g(x,y) = y - x^4
g (generic function with 1 method)

julia> diffTaylor( f(x,y), 1 )   # partial derivative with respect to 1st variable
 - 7.0 + 3.0 x² + 4.0 x y + 𝒪(‖x‖¹¹)

julia> diffTaylor( g(x,y), 2 )
 1.0 + 𝒪(‖x‖¹¹)

julia> diffTaylor( g(x,y), 3 )   # error, since we are dealing with 2 variables
ERROR: AssertionError: 1 <= r <= get_numvars()
 in diffTaylor at /Users/benet/.julia/v0.4/TaylorSeries/src/TaylorN.jl:774
 in diffTaylor at /Users/benet/.julia/v0.4/TaylorSeries/src/TaylorN.jl:804

```

[`evaluate`]({ref}) can also be used for [`TaylorN`]({ref}) objects, using
it on vectors of
numbers (`Real` or `Complex`); the length of the vector must coincide with the number
of independent variables.

```julia
julia> evaluate(exy, [.1,.02]) == e^0.12
true

```

Functions to compute the gradient, Jacobian and
Hessian have also been implemented. Using the
functions $f(x,y) = x^3 + 2x^2 y - 7 x + 2$ and $g(x,y) = y-x^4$ defined above,
we may use [`gradient()`]({ref}) or `∇` (`\nabla+TAB`); the results are of
type `Array{TaylorN{T},1}`. To compute the Jacobian and Hessian of a vector field
evaluated at a point, we use respectively [`jacobian()`]({ref}) and [`hessian()`]({ref}):

```julia
julia> f1 = f(x,y)
 2.0 - 7.0 x + 1.0 x³ + 2.0 x² y + 𝒪(‖x‖¹¹)

julia> g1 = g(x,y)
 1.0 y - 1.0 x⁴ + 𝒪(‖x‖¹¹)

julia> gradient( g1 )
2-element Array{TaylorSeries.TaylorN{Float64},1}:
  - 4.0 x³ + 𝒪(‖x‖¹¹)
       1.0 + 𝒪(‖x‖¹¹)

julia> ∇(f1)
2-element Array{TaylorSeries.TaylorN{Float64},1}:
  - 7.0 + 3.0 x² + 4.0 x y + 𝒪(‖x‖¹¹)
                    2.0 x² + 𝒪(‖x‖¹¹)

julia> fg = f1-g1-2*f1*g1
 2.0 - 7.0 x - 5.0 y + 14.0 x y + 1.0 x³ + 2.0 x² y + 5.0 x⁴ - 2.0 x³ y - 4.0 x² y² - 14.0 x⁵ + 2.0 x⁷ + 4.0 x⁶ y + 𝒪(‖x‖¹¹)

julia> jacobian([f1,g1], [2,1])
2x2 Array{Float64,2}:
  13.0  8.0
 -32.0  1.0

julia> hessian(fg, [1.0,1.0])
2x2 Array{Float64,2}:
 -26.0  20.0
  20.0  -8.0

```
