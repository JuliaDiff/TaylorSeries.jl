<script type="text/x-mathjax-config">
  MathJax.Hub.Config({
    TeX: { equationNumbers: { autoNumber: "AMS" } }
  });
  MathJax.Hub.Config({
    TeX: { extensions: ["AMSmath.js", "AMSsymbols.js", "autobold.js", "autoload-all.js"] }
  });
  MathJax.Hub.Config({
    tex2jax: {
      inlineMath: [['$','$'], ['\\(','\\)']],
      processEscapes: true
    }
  });
</script>
<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS_HTML">
</script>

# User guide

---

`TaylorSeries.jl` can be thought of as a polynomial algebraic manipulator in one or more
variables; these two cases are treated separately.  Three new types are defined,
`Taylor1`, `HomogeneousPolynomial` and `TaylorN`, which correspond to
expansions in one independent variable, homogeneous polynomials of various variables, and the polynomial
series in many independent variables, respectively. These types are subtypes
of `Number` and are defined parametrically.

The package is loaded as usual:

```julia
julia> using TaylorSeries
```

## One variable

Taylor expansions in one variable are represented by the `Taylor1` type, which
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
 1 + 2⋅t + 3⋅t² + 𝒪(t³)

julia> Taylor1([0.0, 1im]) # Also works with complex numbers
 ( 1.0 im )⋅t + 𝒪(t²)

julia> affine(a) = a + taylor1_variable(typeof(a),5)  ## a + t of order 5
affine (generic function with 1 method)

julia> t = affine(0.0) # Independent variable `t`
 1.0⋅t + 𝒪(t⁶)
```
Note that the information about the maximum order considered is displayed
using a big-O notation.

The definition of `affine(a)` uses the function `taylor1_variable`, which is a
shortcut to define the independent variable of a Taylor expansion,
with a given type and given order. As we show below, this is one of the
easiest ways to work with the package.

The usual arithmetic operators (`+`, `-`, `*`, `/`, `^`, `==`) have been
extended to work with the `Taylor1` type, including promotions that involve
`Number`s. The operations return a valid Taylor expansion with the same
maximum order; compare the last example below, where this is not possible:

```julia
julia> t*(3t+2.5)
 2.5⋅t + 3.0⋅t² + 𝒪(t⁶)

julia> 1/(1-t)
 1.0 + 1.0⋅t + 1.0⋅t² + 1.0⋅t³ + 1.0⋅t⁴ + 1.0⋅t⁵ + 𝒪(t⁶)

julia> t*(t^2-4)/(t+2)
 - 2.0⋅t + 1.0⋅t² + 𝒪(t⁶)

julia> tI = im*t
 ( 1.0 im )⋅t + 𝒪(t⁶)

julia> t^6  # order is 5
 0.0 + 𝒪(t⁶)

julia> (1-t)^3.2
 1.0 - 3.2⋅t + 3.5200000000000005⋅t² - 1.4080000000000004⋅t³ + 0.07040000000000009⋅t⁴ + 0.011264000000000012⋅t⁵ + 𝒪(t⁶)

julia> (1+t)^t
 1.0 + 1.0⋅t² - 0.5⋅t³ + 0.8333333333333333⋅t⁴ - 0.75⋅t⁵ + 𝒪(t⁶)

julia> t^3.2
ERROR: The 0th order Taylor1 coefficient must be non-zero
to raise the Taylor1 polynomial to a non-integer exponent
 in ^ at /Users/benet/Fisica/6-IntervalArithmetics/TaylorSeries.jl/src/utils_Taylor1.jl:280
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
 1.0 + 1.0⋅t + 0.5⋅t² + 0.16666666666666666⋅t³ + 0.041666666666666664⋅t⁴ + 0.008333333333333333⋅t⁵ + 𝒪(t⁶)

julia> log(1-t)
 - 1.0⋅t - 0.5⋅t² - 0.3333333333333333⋅t³ - 0.25⋅t⁴ - 0.2⋅t⁵ + 𝒪(t⁶)

julia> sqrt(t)
ERROR: First non-vanishing Taylor1 coefficient must correspond
to an **even power** in order to expand `sqrt` around 0
 in sqrt at /Users/benet/Fisica/6-IntervalArithmetics/TaylorSeries.jl/src/utils_Taylor1.jl:351

julia> sqrt(1 + t)
 1.0 + 0.5⋅t - 0.125⋅t² + 0.0625⋅t³ - 0.0390625⋅t⁴ + 0.02734375⋅t⁵ + 𝒪(t⁶)

julia> imag(exp(tI)')
 - 1.0⋅t + 0.16666666666666666⋅t³ - 0.008333333333333333⋅t⁵ + 𝒪(t⁶)

julia> real(exp(Taylor1([0.0,1im],17))) - cos(Taylor1([0.0,1.0],17)) == 0.0
true

julia> convert(Taylor1{Rational{Int64}}, exp(t))  # output differes in v0.4
 1//1 + 1//1⋅t + 1//2⋅t² + 1//6⋅t³ + 1//24⋅t⁴ + 1//120⋅t⁵ + 𝒪(t⁶)
```

Differentiating and integrating is straightforward for polynomial expansions in
one variable. The last coefficient of a derivative is set to zero to keep the
same order as the original polynomial; for the integral, an
integration constant may be set to a different value (the default is zero). The
order of the resulting polynomial is not changed. The $n$-th ($n \ge 0$)
derivative is obtained using `deriv(a,n)`, where `a` is a Taylor series;
the default is $n=1$.

```julia
julia> diffTaylor(exp(t))
 1.0 + 1.0⋅t + 0.5⋅t² + 0.16666666666666666⋅t³ + 0.041666666666666664⋅t⁴ + 𝒪(t⁶)

julia> integTaylor(exp(t))
 1.0⋅t + 0.5⋅t² + 0.16666666666666666⋅t³ + 0.041666666666666664⋅t⁴ + 0.008333333333333333⋅t⁵ + 𝒪(t⁶)

julia> integTaylor( ans, 1.0)
 1.0 + 0.5⋅t² + 0.16666666666666666⋅t³ + 0.041666666666666664⋅t⁴ + 0.008333333333333333⋅t⁵ + 𝒪(t⁶)

julia> integTaylor( diffTaylor( exp(-t)), 1.0 ) == exp(-t)
true

julia> deriv( exp(affine(1.0))) == exp(1.0)
true

julia> deriv( exp(affine(1.0)), 5) == exp(1.0) # Fifth derivative of `exp(1+t)`
true
```

To evaluate a Taylor series at a point, Horner's rule is used via the function
`evalTaylor(a::Taylor, dt::Number)`. Here, $dt$ is the increment from
the point $t_0$ where the Taylor expansion is calculated, i.e., the series
is evaluated at $t = t_0 + dt$. Omitting $dt$ corresponds to $dt = 0$.

```julia
julia> evalTaylor(exp(affine(1.0))) - e # exp(t) around t0=1 (order 5), evaluated there (dt=0)
0.0

julia> evalTaylor(exp(t), 1) - e # exp(t) around t0=0 (order 5), evaluated at t=1
-0.0016151617923783057

julia> evalTaylor( exp( taylor1_variable(17) ), 1) - e # exp(t) around t0=0 (order 17), evaluated at t=1
0.0

julia> tBig = Taylor1([zero(BigFloat),one(BigFloat)],50) # With BigFloats
 1e+00⋅t + 𝒪(t⁵¹)

julia> evalTaylor( exp(tBig), one(BigFloat) )
2.718281828459045235360287471352662497757247093699959574966967627723419298053556e+00 with 256 bits of precision

julia> e - ans
6.573322999985292556154129119543257102601105719980995128942636339920549561322098e-67 with 256 bits of precision
```


## Many variables

A polynomial in $N>1$ variables can be represented in two distinct ways:
As a vector whose coefficients are homogeneous polynomials of fixed degree, or
as a vector whose coefficients are polynomials in $N-1$ variables. We have opted
to implement the first option, which seems to show better performance. An elegant
(lazy) implementation of the second representation was discussed on the
[julia-users](https://groups.google.com/forum/#!msg/julia-users/AkK_UdST3Ig/sNrtyRJHK0AJ) list.

`TaylorN` is thus constructed as a vector of parameterized homogeneous polynomials
defined by the type `HomogeneousPolynomial`, which in turn is a vector of
coefficients of given order (degree). This implementation imposes that the user
has to specify the maximum order and the maximum number of independent
variables. This is done using `set_params_TaylorN(maxord, numVars)`, where
`maxord` is the maximum order (degree)
of the polynomial considered, and `numVars` is the number of variables;
`set_maxOrder()` and `set_numVars()` allow to change independently each
parameter. By default, these parameters are set to 6 and 2, respectively. The
function ``show_params_TaylorN()` displays the current values of these parameters.

```julia
julia> show_params_TaylorN()
INFO: `TaylorN` and `HomogeneousPolynomial` parameters:
    Maximum order       = 6
    Number of variables = 2

julia> set_params_TaylorN(8,2)
Warning: redefining constant _params_taylorN
(8,2)
```

Technically (internally), these values define dictionaries that
translate the position of the coefficients of a `HomogeneousPolynomial`
into the corresponding
multi-variable monomials. Fixing these values from the start is imperative.

The easiest way of constructing a `TaylorN` object is by defining a symbol for
the independent variables, using the function
`taylorN_variable(T::Type{T<:Top}, nv::Int64, order::Int64)`;
omitting the type yields a `TaylorN{Float64}` object, and omitting the order
defaults to the maximum order. Again, the Taylor expansions are implemented
around 0 for all variables; if the expansion
is needed around a different value, the trick is a simple translation of
the corresponding
independent variable $x \to x+a$.

```julia
julia> x = taylorN_variable(1)
 1.0⋅x₁ + 𝒪(‖x‖⁹)

julia> y = taylorN_variable(2,6)
 1.0⋅x₂ + 𝒪(‖x‖⁷)

julia> typeof(x)
TaylorN{Float64} (constructor with 1 method)

julia> x.order
8

julia> x.coeffs
9-element Array{HomogeneousPolynomial{Float64},1}:
     0.0
  1.0⋅x₁
     0.0
     0.0
     0.0
     0.0
     0.0
     0.0
     0.0

julia> get_maxOrder(y)
6
```

Other ways of constructing `TaylorN` polynomials involve using `HomogeneousPolynomial`
objects directly, which is uncomfortable:

```julia
julia> HomogeneousPolynomial([1,-1])
 1⋅x₁ - 1⋅x₂

julia> TaylorN( [HomogeneousPolynomial([1,0]), HomogeneousPolynomial([1,2,3])], 4)
 1⋅x₁ + 1⋅x₁² + 2⋅x₁⋅x₂ + 3⋅x₂² + 𝒪(‖x‖⁵)
```

As before, the usual arithmetic operators (`+`, `-`, `*`, `/`, `^`, `==`)
have been extended to work with `TaylorN` objects, including the appropriate
promotions to deal with numbers. (Some of the arithmetic operations have
also been extended for
`HomogeneousPolynomial`, whenever the result is a `HomogeneousPolynomial`;
division, for instance, is not extended.) Also, the elementary functions have been
implemented, again by computing their coefficients recursively:

```julia
julia> exy = exp(x+y)
 1.0 + 1.0⋅x₁ + 1.0⋅x₂ + 0.5⋅x₁² + 1.0⋅x₁⋅x₂ + 0.5⋅x₂² + 0.16666666666666666⋅x₁³ + 0.5⋅x₁²⋅x₂ + 0.5⋅x₁⋅x₂² + 0.16666666666666666⋅x₂³ + 0.041666666666666664⋅x₁⁴ + 0.16666666666666666⋅x₁³⋅x₂ + 0.25⋅x₁²⋅x₂² + 0.16666666666666666⋅x₁⋅x₂³ + 0.041666666666666664⋅x₂⁴ + 0.008333333333333333⋅x₁⁵ + 0.041666666666666664⋅x₁⁴⋅x₂ + 0.08333333333333333⋅x₁³⋅x₂² + 0.08333333333333333⋅x₁²⋅x₂³ + 0.041666666666666664⋅x₁⋅x₂⁴ + 0.008333333333333333⋅x₂⁵ + 0.0013888888888888887⋅x₁⁶ + 0.008333333333333331⋅x₁⁵⋅x₂ + 0.020833333333333332⋅x₁⁴⋅x₂² + 0.027777777777777776⋅x₁³⋅x₂³ + 0.020833333333333332⋅x₁²⋅x₂⁴ + 0.008333333333333331⋅x₁⋅x₂⁵ + 0.0013888888888888887⋅x₂⁶ + 0.00019841269841269839⋅x₁⁷ + 0.0013888888888888885⋅x₁⁶⋅x₂ + 0.004166666666666666⋅x₁⁵⋅x₂² + 0.006944444444444443⋅x₁⁴⋅x₂³ + 0.006944444444444443⋅x₁³⋅x₂⁴ + 0.004166666666666666⋅x₁²⋅x₂⁵ + 0.0013888888888888885⋅x₁⋅x₂⁶ + 0.00019841269841269839⋅x₂⁷ + 2.4801587301587298e-5⋅x₁⁸ + 0.00019841269841269836⋅x₁⁷⋅x₂ + 0.0006944444444444443⋅x₁⁶⋅x₂² + 0.0013888888888888887⋅x₁⁵⋅x₂³ + 0.0017361111111111108⋅x₁⁴⋅x₂⁴ + 0.0013888888888888887⋅x₁³⋅x₂⁵ + 0.0006944444444444443⋅x₁²⋅x₂⁶ + 0.00019841269841269836⋅x₁⋅x₂⁷ + 2.4801587301587298e-5⋅x₂⁸ + 𝒪(‖x‖⁹)
```

Note above that `y` has been promoted internally so the result corresponds
to the maximum of the orders of `x` and `y`. The function `get_coeff(a,v)`
gives the coefficient of `x` that corresponds to the monomial
specified by the vector of powers `v`:
```julia
julia> get_coeff(exy, [3,5]) == 1/720
false

julia> rationalize(get_coeff(exy, [3,5]))
1//720
```

Partial differentiation is also implemented for
`TaylorN` objects,
using `diffTaylor`; integration is yet to be implemented.

```julia
julia> f(x,y) = x^3 + 2x^2 * y - 7x + 2
f (generic function with 1 method)

julia> g(x,y) = y - x^4
g (generic function with 1 method)

julia> diffTaylor( f(x,y), 1 )   # partial derivative with respect to 1st variable
 - 7.0 + 3.0⋅x₁² + 4.0⋅x₁⋅x₂ + 𝒪(‖x‖⁹)

julia> diffTaylor( g(x,y), 2 )
 1.0 + 𝒪(‖x‖⁹)

julia> diffTaylor( g(x,y), 3 )   # error, since we are dealing with 2 variables
ERROR: assertion failed: 1 <= r <= _params_taylorN.numVars
 in diffTaylor at /Users/benet/Fisica/6-IntervalArithmetics/TaylorSeries.jl/src/utils_TaylorN.jl:679
 in diffTaylor at /Users/benet/Fisica/6-IntervalArithmetics/TaylorSeries.jl/src/utils_TaylorN.jl:709
```

`evalTaylor` can also be used for `TaylorN` objects, both for vectors of
numbers or vectors of `Taylor1` objects; the length of the vector must
coincide with the number
of independent variables. The implementation still needs some improvements.

```julia
julia> evalTaylor(x+y, [t, 2t])  # x+y, with x=t, y=2t
 3.0⋅t + 𝒪(t⁶)

julia> evalTaylor(exy, [t,2t])
 1.0 + 3.0⋅t + 4.5⋅t² + 4.5⋅t³ + 3.375⋅t⁴ + 2.025⋅t⁵ + 𝒪(t⁶)

julia> exp(3.0*taylor1_variable(8))
 1.0 + 3.0⋅t + 4.5⋅t² + 4.5⋅t³ + 3.375⋅t⁴ + 2.025⋅t⁵ + 1.0125⋅t⁶ + 0.4339285714285714⋅t⁷ + 0.16272321428571426⋅t⁸ + 𝒪(t⁹)
```

Functions to compute the gradient, Jacobian and
Hessian have also been implemented. Using the
functions $f(x,y) = x^3 + 2x^2 y - 7 x + 2$ and $g(x,y) = y-x^4$ defined above,
we may use `∇` (`\nabla+TAB`) or `TaylorSeries.gradient`; the results are of
type `Array{TaylorN{T},1}`. To compute the Jacobian or Hessian of a vector field
evaluated at a point, we use `jacobian` and `hessian`:

```julia
julia> f1 = f(x,y)
 2.0 - 7.0⋅x₁ + 1.0⋅x₁³ + 2.0⋅x₁²⋅x₂ + 𝒪(‖x‖⁹)

julia> g1 = g(x,y)
 1.0⋅x₂ - 1.0⋅x₁⁴ + 𝒪(‖x‖⁹)

julia> ∇(f1)
2-element Array{TaylorN{Float64},1}:
  - 7.0 + 3.0⋅x₁² + 4.0⋅x₁⋅x₂ + 𝒪(‖x‖⁹)
                      2.0⋅x₁² + 𝒪(‖x‖⁹)

julia> gradient( g1 )
2-element Array{TaylorN{Float64},1}:
  - 4.0⋅x₁³ + 𝒪(‖x‖⁹)
        1.0 + 𝒪(‖x‖⁹)

julia> jacobian([f1,g1], [2,1])
2x2 Array{Float64,2}:
  13.0  8.0
 -32.0  1.0

julia> fg = f1-g1-2*f1*g1
 2.0 - 7.0⋅x₁ - 5.0⋅x₂ + 14.0⋅x₁⋅x₂ + 1.0⋅x₁³ + 2.0⋅x₁²⋅x₂ + 5.0⋅x₁⁴ - 2.0⋅x₁³⋅x₂ - 4.0⋅x₁²⋅x₂² - 14.0⋅x₁⁵ + 2.0⋅x₁⁷ + 4.0⋅x₁⁶⋅x₂ + 𝒪(‖x‖⁹)

julia> hessian(ans)
2x2 Array{Float64,2}:
  0.0  14.0
 14.0   0.0

julia> evalTaylor(fg, [x+1.0, y+1.0])
 - 2.0 - 12.0⋅x₁ + 5.0⋅x₂ - 13.0⋅x₁² + 20.0⋅x₁⋅x₂ - 4.0⋅x₂² + 29.0⋅x₁³ + 48.0⋅x₁²⋅x₂ - 8.0⋅x₁⋅x₂² + 65.0⋅x₁⁴ + 78.0⋅x₁³⋅x₂ - 4.0⋅x₁²⋅x₂² + 52.0⋅x₁⁵ + 60.0⋅x₁⁴⋅x₂ + 18.0⋅x₁⁶ + 24.0⋅x₁⁵⋅x₂ + 2.0⋅x₁⁷ + 4.0⋅x₁⁶⋅x₂ + 𝒪(‖x‖⁹)

julia> hessian(fg, [1.0,1.0])
2x2 Array{Float64,2}:
 -26.0  20.0
  20.0  -8.0
```


## Examples

### 1. Four-square identity

The first example described below, shows that the four-square identity holds:
\\begin{eqnarray}
(a_1+a_2+a_3+a_4)\\cdot(b_1+b_2+b_3+b_4) & = &
     (a_1 b_1 - a_2 b_2 - a_3 b_3 -a_4 b_4)^2 + \\qquad \\nonumber \\\\
\\label{eq:Euler}
  & & (a_1 b_2 - a_2 b_1 - a_3 b_4 -a_4 b_3)^2 + \\\\
  & & (a_1 b_3 - a_2 b_4 - a_3 b_1 -a_4 b_2)^2 + \\nonumber \\\\
  & & (a_1 b_4 - a_2 b_3 - a_3 b_2 -a_4 b_1)^2, \\nonumber
\\end{eqnarray}
as proved by Euler. The code can we found in one of the tests of the package.

First, we reset the maximum degree of the polynomial to 4, since the RHS
of the equation
has *a priori* terms of fourth order, and the number of independent variables to
8.

```julia
julia> set_params_TaylorN(4,8)
Warning: redefining constant _params_taylorN
(4,8)
```

Next, we define the 8 independent variables

```julia
julia> for i=1:4
           ai = symbol(string("a",i))
           bi = symbol(string("b",i))
           @eval ($ai) = taylorN_variable(Int,$i,4)
           @eval ($bi) = taylorN_variable(Int,4+($i),4)
       end

julia> a1
 1⋅x₁ + 𝒪(‖x‖⁵)

julia> b1
 1⋅x₅ + 𝒪(‖x‖⁵)
```

followed by the distinct terms that appear in (\\ref{eq:Euler}):
```julia
julia> expr_lhs1 = a1^2 + a2^2 + a3^2 + a4^2 ;

julia> expr_lhs2 = b1^2 + b2^2 + b3^2 + b4^2 ;

julia> expr_rhs1 = (a1*b1 - a2*b2 - a3*b3 - a4*b4)^2 ;

julia> expr_rhs2 = (a1*b2 + a2*b1 + a3*b4 - a4*b3)^2 ;

julia> expr_rhs3 = (a1*b3 - a2*b4 + a3*b1 + a4*b2)^2 ;

julia> expr_rhs4 = (a1*b4 + a2*b3 - a3*b2 + a4*b1)^2 ;
```
and, finally, check that the LHS is equal to the RHS:
```julia
julia> lhs = expr_lhs1 * expr_lhs2
 1⋅x₁²⋅x₅² + 1⋅x₂²⋅x₅² + 1⋅x₃²⋅x₅² + 1⋅x₄²⋅x₅² + 1⋅x₁²⋅x₆² + 1⋅x₂²⋅x₆² + 1⋅x₃²⋅x₆² + 1⋅x₄²⋅x₆² + 1⋅x₁²⋅x₇² + 1⋅x₂²⋅x₇² + 1⋅x₃²⋅x₇² + 1⋅x₄²⋅x₇² + 1⋅x₁²⋅x₈² + 1⋅x₂²⋅x₈² + 1⋅x₃²⋅x₈² + 1⋅x₄²⋅x₈² + 𝒪(‖x‖⁵)

julia> rhs = expr_rhs1 + expr_rhs2 + expr_rhs3 + expr_rhs4
 1⋅x₁²⋅x₅² + 1⋅x₂²⋅x₅² + 1⋅x₃²⋅x₅² + 1⋅x₄²⋅x₅² + 1⋅x₁²⋅x₆² + 1⋅x₂²⋅x₆² + 1⋅x₃²⋅x₆² + 1⋅x₄²⋅x₆² + 1⋅x₁²⋅x₇² + 1⋅x₂²⋅x₇² + 1⋅x₃²⋅x₇² + 1⋅x₄²⋅x₇² + 1⋅x₁²⋅x₈² + 1⋅x₂²⋅x₈² + 1⋅x₃²⋅x₈² + 1⋅x₄²⋅x₈² + 𝒪(‖x‖⁵)

julia> lhs == rhs
true
```
The identity is satisfied $\\square$.

### 2. Fateman's test

Richard J. Fateman, from Berkley, proposed as a stringent test
of polynomial multiplication
the evaluation of $s*(s+1)$, where $s = (1+x+y+z+w)^{20}$. This is
implemented in
the function `fateman1`. We shall also evaluate the form $s^2+s$ in `fateman2`,
which involves fewer operations (and makes a fairer comparison to what
Mathematica does). Below we use Julia v0.4, since it is faster than v0.3.

```julia
julia> set_params_TaylorN(40,4)   # maxOrder = 40; numVars = 4
(40,4)

julia> function fateman1(ndeg::Int)
           T = Int128
           unoH = HomogeneousPolynomial(one(T), 0)
           s = TaylorN( [unoH, HomogeneousPolynomial([one(T),one(T),one(T),one(T)],1)], ndeg )
           s = s^ndeg
           # s is converted to order 2*ndeg
           s = TaylorN(s, 2ndeg)
           return s * (s+TaylorN(unoH, 2*ndeg))
       end
fateman1 (generic function with 1 method)

julia> @time f1 = fateman1(0);
   7.379 milliseconds (8346 allocations: 361 KB)

julia> @time f1 = fateman1(20);
   9.339 seconds      (6150 allocations: 58698 KB, 0.07% gc time)

julia> get_coeff(f1,[1,6,7,20])
128358585324486316800

julia> ans > typemax(Int)  # this is the reason for using Int128
true

julia> function fateman2(ndeg::Int)
           T = Int128
           unoH = HomogeneousPolynomial(one(T), 0)
           s = TaylorN( [unoH, HomogeneousPolynomial([one(T),one(T),one(T),one(T)],1)], ndeg )
           s = s^ndeg
           # s is converted to order 2*ndeg
           s = TaylorN(s, 2ndeg)
           return s^2 + s
       end
fateman2 (generic function with 1 method)

julia> @time f2 = fateman2(0);
  10.840 milliseconds (6805 allocations: 297 KB, 41.42% gc time)

julia> @time f2 = fateman2(20);
   4.737 seconds      (4870 allocations: 42049 KB, 0.09% gc time)

julia> get_coeff(f2,[1,6,7,20])
128358585324486316800

julia> sum(TaylorSeries.sizeTable) # number of distinct monomials
135751
```

The tests above show the necessity of using integers of type `Int128`, and that
`fateman2` is
about twice as fast as `fateman1`. We
mention that our implementation of `fateman2` is roughly
1.5 times slower than Mathematica, which takes 3.22 sec.
