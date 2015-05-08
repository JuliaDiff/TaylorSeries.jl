# TaylorSeries.jl

A Julia package for Taylor expansions in one or more independent variables.

[![Build Status](https://travis-ci.org/lbenet/TaylorSeries.jl.svg?branch=master)](https://travis-ci.org/lbenet/TaylorSeries.jl)
[![TaylorSeries](http://pkg.julialang.org/badges/TaylorSeries_nightly.svg)](http://pkg.julialang.org/?pkg=TaylorSeries&ver=nightly)

#### Authors
- Luis Benet, Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)
- David P. Sanders, Facultad de Ciencias, Universidad Nacional Autónoma de México (UNAM)

Comments on the code and suggestions are welcome and appreciated.

### Basic of usage

This package computes Taylor expansions using automatic differentiation techniques (see e.g. Alex Haro's [nice paper][1] for a description of the recurrence relations used). The module implements the basic operations (+,-,*,/,^), some elementary functions (exp, log, sin, cos, tan), integration, differentiation and evaluation (Horner's rule) of the series on a point. The expansions work for one or more independent variables.

Installing `TaylorSeries` is done as usual with `Pkg.add("TaylorSeries")`.

The package introduces three new types, `Taylor1`, `HomogeneousPolynomial` and `TaylorN`, for the expansions in one and more independent-variable expansions, respectively. These are subtypes of `Number`.

`Taylor1` is defined by one of the constructors

- `Taylor1{T<:Number}(a::Taylor1{T}, order::Integer)`
- `Taylor1{T<:Number}(a::Array{T,1}, order::Integer)`
- `Taylor1{T<:Number}(a::T, order::Integer)`

Thus, one-variable Taylor polynomials can be constructed from its coefficients (increasing with the order of the polynomial term), which may be of type `Taylor1`, a vector or simply a number (`<:Number`). The order of the polynomial can be omitted, which is then fixed from the length of the coefficients; otherwise it is the `max` among the length of the coefficients or the integer specified. In order to simplify defining an independent `Taylor1` variable, we have defined the shortcut `taylor1_variable(T::Type, order::Int=1 )`.


`TaylorN` is constructed as a vector of homogeneous polynomials defined by the type `HomogeneousPolynomial`, which is a vector of coefficients of given order (degree). (This new implementation is the **major difference** with respect to the previous released version (v0.0.3), with all adaptations that this implies; this change results in much better performance.) The functions `set_maxOrder()` and `set_numVars()`set the overall constants defining the maximum order (maximum degree) of the expansions and the number of variables, respectively; `set_params_TaylorN` is more efficient for setting both parameters at once. Setting this parameters is *essential* to properly construct the dictionaries that translate the position of the coefficients of `HomogeneousPolynomial` into the corresponding multi-variable monomials. It is a good idea to fix them before defining the `TaylorN` or `HomogeneousPolynomial` objects, to avoid unpleasant surprises. This implementation is
certainly not the laziest one, though it achieves very good performance.

`TaylorN` can be thus constructed by:

- `TaylorN{T<:Number}(x::TaylorN{T}, order::Int)`
- `TaylorN{T<:Number}(v::Array{HomogeneousPolynomial{T},1}, order::Int)`
- `TaylorN{T<:Number}(x::HomogeneousPolynomial{T}, order::Int)`
- `TaylorN{T<:Number}(x::T, order::Int)`

Again, the order can be omitted, in which case it is fixed by the `HomogeneousPolynomial` of maximum order. In order to simplify defining the independent variables, we have defined `taylorN_variable(T::Type, nv::Int, order::Int=1 )`.

In the examples shown below, the output corresponds to Julia version 0.3; `t` is used as the independent `Taylor1` variable whereas `x_{i}` is used for the i-th variable of `HomogeneousPolynomial` and `TaylorN`.

By default, the Taylor expansions are implemented around 0; if the expansion is needed around a different value $ x_0 $, the trick is a simple translation $ x \to x+a $; see the definition of `xT(a)` below.

```julia
julia> using TaylorSeries

help?> Taylor1
INFO: Loading help data...
DataType   : Taylor1{T<:Number} (constructor with 6 methods)
  supertype: Number
  fields   : (:coeffs,:order)

julia> Taylor1([0,1,2])
 1⋅t + 2⋅t² 

julia> dump(ans)
Taylor1{Int64}
  coeffs: Array(Int64,(3,)) [0,1,2]
  order: Int64 2

julia> Taylor1([0.0,1im])
 ( 1.0 im )⋅t

julia> xT(a) = Taylor1([a, one(a)], 5) # expansion (of order 5) around an arbitrary value `a`
xT (generic function with 1 method)

julia> xT(1.0)
 1.0 + 1.0⋅t

julia> xT0 = taylor1_variable(5) # abbreviation for the independent variable for expansions around 0
 1.0⋅t

help?> TaylorN
DataType   : TaylorN{T<:Number} (constructor with 8 methods)
  supertype: Number
  fields   : (:coeffs,:order)

julia> set_maxOrder(6)
(6,2)

julia> get_numVars()
2

julia> show_params_TaylorN()
INFO: `TaylorN` and `HomogeneousPolynomial` parameters:
    Maximum order       = 6
    Number of variables = 2

julia> HomogeneousPolynomial([1,-1])
 1⋅x₁ - 1⋅x₂

julia> TaylorN( [HomogeneousPolynomial([1,0]), HomogeneousPolynomial([1,2,3])], 4)
 1⋅x₁ + 1⋅x₁² + 2⋅x₁⋅x₂ + 3⋅x₂²

julia> xTN(a) = a + taylorN_variable(typeof(a),1,6)  # shifted variable 1 of TaylorN of order 6 (default is get_maxOrder())
xTN (generic function with 1 method)

julia> yTN(a) = a + taylorN_variable(typeof(a),2)
yTN (generic function with 1 method)

julia> yTN(1.0pi)
 3.141592653589793 + 1.0⋅x₂

```

The usual arithmetic operators (`+`, `-`, `*`, `/`, `^`, `==`) have been extended to work with `Taylor1` and `TaylorN` type, including combinations of polynomials and numbers. (Some of the arithmetics operations have also been extended for `HomogeneousPolynomial`, whenever the result is a `HomogeneousPolynomial`.)

```julia
julia> xT0*(3xT0+2.5)
 2.5⋅t + 3.0⋅t² 

julia> 1/(1-xT0)
 1.0 + 1.0⋅t + 1.0⋅t² + 1.0⋅t³ + 1.0⋅t⁴ + 1.0⋅t⁵ 

julia> xT0*(xT0^2-4)/(xT0+2)
 - 2.0⋅t + 1.0⋅t² 

julia> Taylor1([1],5) == 1.0
true

julia> xTI = 1im*xT(0)
 ( 1 im )⋅t

julia> xTI'
 - ( 1 im )⋅t

julia> xTI^4 / xT0^4
  ( 1.0 )

julia> (1-xT0)^3.2
 1.0 - 3.2⋅t + 3.5200000000000005⋅t² - 1.4080000000000004⋅t³ + 0.07040000000000009⋅t⁴ + 0.011264000000000012⋅t⁵ 

julia> (1+xT0)^xT0
 1.0 + 1.0⋅t² - 0.5⋅t³ + 0.8333333333333333⋅t⁴ - 0.75⋅t⁵ 

julia> xTN(1)-yTN(1) # 1+x-(1+y)
 1⋅x₁ - 1⋅x₂

julia> xTN(1)*yTN(1)  # (1+x)*(1+y)
 1 + 1⋅x₁ + 1⋅x₂ + 1⋅x₁⋅x₂

julia> 1/xTN(1)
 1.0 - 1.0⋅x₁ + 1.0⋅x₁² - 1.0⋅x₁³ + 1.0⋅x₁⁴ - 1.0⋅x₁⁵ + 1.0⋅x₁⁶

```

Some elemental functions have been implemented using automatic differentiation techniques, i.e., by computing recursively their coefficients. Some examples include `exp`, `log`, `sqrt`, `sin`, `cos` and `tan`; more functions will be added in the future.

```julia
julia> exp(xT0)
 1.0 + 1.0⋅t + 0.5⋅t² + 0.16666666666666666⋅t³ + 0.041666666666666664⋅t⁴ + 0.008333333333333333⋅t⁵ 

julia> convert(Taylor1{Rational{Int64}}, ans)
 1//1 + 1//1⋅t + 1//2⋅t² + 1//6⋅t³ + 1//24⋅t⁴ + 1//120⋅t⁵ 

julia> log(1-xT0)
 - 1.0⋅t - 0.5⋅t² - 0.3333333333333333⋅t³ - 0.25⋅t⁴ - 0.2⋅t⁵ 

julia> sqrt(xT0)
ERROR: First non-vanishing Taylor1 coefficient must correspond to an **even power**
to expand `sqrt` around 0.
 in sqrt at /Users/benet/Fisica/6-IntervalArithmetics/TaylorSeries/src/utils_Taylor1.jl:339

julia> sqrt( xT(1.0) )   # identic to sqrt( 1+x0 )
 1.0 + 0.5⋅t - 0.125⋅t² + 0.0625⋅t³ - 0.0390625⋅t⁴ + 0.02734375⋅t⁵ 

julia> imag(exp(xTI)')
 - 1.0⋅t + 0.16666666666666666⋅t³ - 0.008333333333333333⋅t⁵ 

julia> real(exp(Taylor1([0.0,1im],17))) - cos(Taylor1([0.0,1.0],17)) == 0.0
true

julia> exp(xTN(0.0)+yTN(0.0))
 1.0 + 1.0⋅x₁ + 1.0⋅x₂ + 0.5⋅x₁² + 1.0⋅x₁⋅x₂ + 0.5⋅x₂² + 0.16666666666666666⋅x₁³ + 0.5⋅x₁²⋅x₂ + 0.5⋅x₁⋅x₂² + 0.16666666666666666⋅x₂³ + 0.041666666666666664⋅x₁⁴ + 0.16666666666666666⋅x₁³⋅x₂ + 0.25⋅x₁²⋅x₂² + 0.16666666666666666⋅x₁⋅x₂³ + 0.041666666666666664⋅x₂⁴ + 0.008333333333333333⋅x₁⁵ + 0.041666666666666664⋅x₁⁴⋅x₂ + 0.08333333333333333⋅x₁³⋅x₂² + 0.08333333333333333⋅x₁²⋅x₂³ + 0.041666666666666664⋅x₁⋅x₂⁴ + 0.008333333333333333⋅x₂⁵ + 0.0013888888888888887⋅x₁⁶ + 0.008333333333333331⋅x₁⁵⋅x₂ + 0.020833333333333332⋅x₁⁴⋅x₂² + 0.027777777777777776⋅x₁³⋅x₂³ + 0.020833333333333332⋅x₁²⋅x₂⁴ + 0.008333333333333331⋅x₁⋅x₂⁵ + 0.0013888888888888887⋅x₂⁶

julia> exp(xTN(0.0)+yTN(0.0)) - exp(xTN(0.0))*exp(yTN(0.0))  # should yield 0
 - 1.734723475976807e-18⋅x₁⁵⋅x₂ - 1.734723475976807e-18⋅x₁⋅x₂⁵

julia> ans == 0.0
true
```

Differentiating and integrating is trivial for polynomial expansions in one variable. The last coefficient of the differential is simply set to zero; for the integral, an integration constant may be set to a different value (the default is zero). The order of the resulting polynomial is not changed. The n-th (n > 0) derivative is obtained with `deriv(a,n)`, where `a` is a Taylor series; default is n=1. Partial differentiation is also implemented for expansion in more than one variable. Note that integration is not yet implemented for `TaylorN`.

```julia
julia> diffTaylor(exp(xT0))
 1.0 + 1.0⋅t + 0.5⋅t² + 0.16666666666666666⋅t³ + 0.041666666666666664⋅t⁴ 

julia> integTaylor(exp(xT0))
 1.0⋅t + 0.5⋅t² + 0.16666666666666666⋅t³ + 0.041666666666666664⋅t⁴ + 0.008333333333333333⋅t⁵ 

julia> integTaylor( ans, 1.0)
 1.0 + 0.5⋅t² + 0.16666666666666666⋅t³ + 0.041666666666666664⋅t⁴ + 0.008333333333333333⋅t⁵

julia> integTaylor( diffTaylor( exp(-xT0)), 1.0 ) == exp(-xT0)
true

julia> deriv( exp(xT(1.0))) == e
true

julia> deriv( exp(xT(1.0)), 5) == e
true

julia> diffTaylor( 3+xTN(0.0)+2yTN(0.0),1 )   # partial derivative with respect to 1st variable
 1.0

julia> diffTaylor( 3+xTN(0.0)+2yTN(0.0),2 )   # partial derivative with respect to 2nd variable
 2.0
```

For the evaluation of a Taylor series, we use Horner's rule in `evalTaylor(a::Taylor1, dx::Number)`. Here, $ dx $ is the difference with respect to the point $ x_0 $ where the Taylor expansion is calculated, i.e., the series is evaluated at $ x = x_0 + dx $. Omitting $ dx $ corresponds to $ dx = 0 $. This function is generalized to admit `TaylorN` objects and vectors of numbers or `Taylor1` objects; the length of these vectors must coincide with the number of independent variables.

```julia
julia> evalTaylor( exp( xT(1.0) )) - e ## exp(x) around x0=1 (order 5), evaluated there (dx=0)
0.0

julia> evalTaylor( exp(xT0), 1.0) - e  ## exp(x) around x0=0 (order 5), evaluated at x=1
-0.0016151617923783057

julia> evalTaylor( exp( taylor1_variable(17)), 1.0) - e  ## exp(x) around x0=0 (order 17), evaluated at x=1
0.0

julia> evalTaylor( exp( taylor1_variable(BigFloat,50) ), one(BigFloat) )
2.718281828459045235360287471352662497757247093699959574966967627723419298053556e+00 with 256 bits of precision

julia> convert(BigFloat,e) - ans
6.573322999985292556154129119543257102601105719980995128942636339920549561322098e-67 with 256 bits of precision

julia> evalTaylor(xTN(0.0)+yTN(0.0), [xT0, 2xT0])   # x+y, with x=x0, y=2x0
 3.0⋅t
```

We have also implemented functions to compute gradients, Jacobians and Hessians. The following computes the gradients of $ f(x,y) = x^3 + 2x^2 y - 7 x + 2 $ and $ g(x,y) = y-x^4 $, using either `∇` (\nabla+TAB) or `TaylorSeries.gradient`; the results are of type `Array{TaylorN{T},1}`. In order to compute the Jacobian of a vector field evaluated on a point we use `jacobian`, and for the Hessian of a function we employ `hessian`.

```julia
julia> f(x,y) = x^3 + 2x^2 * y - 7x + 2
f (generic function with 1 method)

julia> g(x,y) = y - x^4
f (generic function with 1 method)

julia> f1 = f(xTN(0),yTN(0))
 2 - 7⋅x₁ + 1⋅x₁³ + 2⋅x₁²⋅x₂

julia> g1 = g(xTN(0),yTN(0))
 1⋅x₂ - 1⋅x₁⁴
 
julia> ∇(f1)
2-element Array{TaylorN{Int64},1}:
2-element Array{TaylorN{Int64},1}:
  - 7 + 3⋅x₁² + 4⋅x₁⋅x₂
                  2⋅x₁²

julia> jacobian([f1,g1], [2,1])
2x2 Array{Int64,2}:
  13  8
 -32  1

julia> hessian(f1-g1-2*f1*g1)
2x2 Array{Int64,2}:
  0  14
 14   0

```

Some concrete applications of the package will be found in the directory [examples][4]; for the time being we have only included [Kepler's problem integrated using Taylor's method][5]).

Finally, it is worth pointing out the existing Julia packages [Polynomial][2] and [PowerSeries][3] have similar functionality as `TaylorSeries` for one-variable expansions using somewhat different approaches.


#### Acknowledgments
This project began (using `python`) during a Masters' course in the postgraduate programs in Physics and in Mathematics at UNAM, during the second half of 2013. We thank the participants of the course for putting up with the half-baked material and contributing energy and ideas.

We acknowledge financial support from DGAPA-UNAM PAPIME grants PE-105911 and PE-107114, and PAPIIT grant IG-101113. LB acknowledges support through a *Cátedra Moshinsky* (2013).

[1]: http://www.maia.ub.es/~alex/admcds/admcds.pdf
[2]: https://github.com/vtjnash/Polynomial.jl
[3]: https://github.com/jwmerrill/PowerSeries.jl
[4]: https://github.com/lbenet/TaylorSeries.jl/tree/master/examples/
[5]: http://nbviewer.ipython.org/github/lbenet/TaylorSeries.jl/blob/master/examples/1-KeplerProblem.ipynb
