# TaylorSeries.jl

A julia package for Taylor expansions in one or more independent variables.

[![Build Status](https://travis-ci.org/lbenet/TaylorSeries.jl.svg?branch=master)](https://travis-ci.org/lbenet/TaylorSeries.jl)
[![TaylorSeries](http://pkg.julialang.org/badges/TaylorSeries_nightly.svg)](http://pkg.julialang.org/?pkg=TaylorSeries&ver=nightly)

#### Authors
- Luis Benet, Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)
- David P. Sanders, Facultad de Ciencias, Universidad Nacional Autónoma de México (UNAM)

Comments on the code and suggestions are welcome and appreciated.

### Basic of usage

This package computes Taylor expansions using automatic differentiation techniques (see e.g. Alex Haro's [nice paper][1] for a description of the recurrence relations used). The module implements the basic operations (+,-,*,/,^), some elementary functions (exp, log, sin, cos, tan), integration, differentiation and evaluation (Horner's rule) of the series on a point. The expansions work for one or more independent variables.

Installing `TaylorSeries` is done as usuas with `Pkg.add("TaylorSeries")`.

The package introduces three new types, `Taylor`, `HomogPol` and `TaylorN`, for the expansions in one and more independent-variable expansions, respectively. These are subtypes of `AbstractSeries{T,N} <: Number`, with `T` the type of the coefficients and `N` the number of independent variables.

`Taylor` is defined by one of the constructors

- `Taylor{T<:Number}(a::Taylor{T}, order::Integer)`
- `Taylor{T<:Number}(a::Array{T,1}, order::Integer)`
- `Taylor{T<:Number}(a::T, order::Integer)`

Thus, one-variable Taylor polynomials can be constructed from its coefficients (increasing with the order of the polynomial term), which may be of type `Taylor`, a vector or simply a number (`<:Number`). The order of the polynomial can be omitted, which is then fixed from the length of the coefficients; otherwise it is the `max` among the length of the coefficients or the integer specified.

`TaylorN` is constructed as a vector of homogeneous polynomials defined by the type `HomogPol`, which is a vector of coefficients of given order (degree). (This new implementation is the **major difference** with respect to the previous version, with all adaptations that this implies; this change results in much better performance.) Two basic functions, `set_maxOrder()` and `set_numVars()`, fix the overall constants defining the maximum order (maximum degree) of the expansions and the number of variables, respectively. These are essential to construct the dictionaries that translate the position of the coefficients of `HomogPol` into the corresponding multi-variable monomials. It is a good idea to fix these values before defining the `TaylorN` or `HomogPol` objects; if the number of variables is changed, all `TaylorN` objects must be redefined. As mentioned above, the coefficients of the homogeneous polynomials `HomogPol` are stored as a vector; this will be improved in the future.

`TaylorN` can be thus constructed by:

- `TaylorN{T<:Number}(x::TaylorN{T}, order::Int)`
- `TaylorN{T<:Number}(v::Array{HomogPol{T},1}, order::Int)`
- `TaylorN{T<:Number}(x::HomogPol{T}, order::Int)`
- `TaylorN{T<:Number}(x::T, order::Int)`

Again, the order can be omitted, in which case it is fixed by the `HomogPol` of maximum order. In order to simplify defining independent variables, we have defined `taylorvar(T::Type, nv::Int, order::Int=1 )`.

By default, the Taylor expansions are implemented around 0; if the expansion is needed around a different value $ x_0 $, the trick is a simple translation $ x \to x+a $; see the definition of `xT(a)` below.

```julia
julia> using TaylorSeries

help?> Taylor   # ?Taylor
INFO: Loading help data...
DataType   : Taylor{T<:Number} (constructor with 6 methods)
  supertype: AbstractSeries{T<:Number,1}
  fields   : (:coeffs,:order)

julia> Taylor([0,1,2])
Taylor{Int64}([0,1,2],2)

julia> Taylor([0.0,1im])
Taylor{Complex{Float64}}(Complex{Float64}[0.0+0.0im,0.0+1.0im],1)

julia> xT(a) = Taylor([a, one(a)], 5) # expansion (of order 5) around an arbitrary value `a`
xT (generic function with 1 method)

julia> xT(1.0)
Taylor{Float64}([1.0,1.0,0.0,0.0,0.0,0.0],5)

julia> xT0 = xT(0.0) # abbreviation of the independent variable for expansions around 0
Taylor{Float64}([0.0,1.0,0.0,0.0,0.0,0.0],5)

help?> TaylorN   # ?TaylorN
DataType   : TaylorN{T<:Number} (constructor with 8 methods)
  supertype: AbstractSeries{T<:Number,2}
  fields   : (:coeffs,:order)

help?> TaylorSeries.AbstractSeries  # ?TaylorSeries.AbstractSeries
DataType   : AbstractSeries{T<:Number,N}
  supertype: Number
  subtypes : {HomogPol{T<:Number},TaylorN{T<:Number},Taylor{T<:Number}}

julia> set_maxOrder(6)
INFO: MAXORDER is now 6; hash tables regenerated.
6

julia> get_numVars()
2

julia> HomogPol([1,-1])
HomogPol{Int64}([1,-1],1)

julia> TaylorN( [HomogPol([1,0]), HomogPol([1,2,3])], 4)
TaylorN{Int64}( [HomogPol{Int64}([1,0],1), HomogPol{Int64}([1,2,3],2)], 4)

julia> xTN(a) = a + taylorvar(typeof(a),1,5)  # shifted variable 1 of TaylorN of order 6 (default is 1)
xTN (generic function with 1 method)

julia> yTN(a) = a + taylorvar(typeof(a),2,5)  # shifted variable 2 of TaylorN of order 6 (default is 1)
yTN (generic function with 1 method)

julia> yTN(1.0pi)
TaylorN{Float64}( [HomogPol{Float64}([3.14159],0), HomogPol{Float64}([0.0,1.0],1)], 5)
```

The usual arithmetic operators (`+`, `-`, `*`, `/`, `^`, `==`) have been extended to work with `Taylor` and `TaylorN` type, including combinations of polynomials and numbers. (Some of the arithmetics operations have also been extended for `HomogPol`, whenever the result is a `HomogPol`.)

```julia
julia> xT0*(3xT0+2.5)
Taylor{Float64}([0.0,2.5,3.0,0.0,0.0,0.0],5)

julia> 1/(1-xT0)
Taylor{Float64}([1.0,1.0,1.0,1.0,1.0,1.0],5)

julia> xT0*(xT0^2-4)/(xT0+2)
Taylor{Float64}([-0.0,-2.0,1.0,0.0,0.0,0.0],5)

julia> Taylor([1],5) == 1.0
true

julia> xTI = Taylor([0,im],5)
Taylor{Complex{Int64}}(Complex{Int64}[0+0im,0+1im,0+0im,0+0im,0+0im,0+0im],5)

julia> xTI'
Taylor{Complex{Int64}}(Complex{Int64}[0+0im,0-1im,0+0im,0+0im,0+0im,0+0im],5)

julia> xTI^4 / xT0^4
Taylor{Complex{Float64}}(Complex{Float64}[1.0+0.0im,0.0+0.0im,0.0+0.0im,0.0+0.0im,0.0+0.0im,0.0+0.0im],5)

julia> (1-xT0)^3.2
Taylor{Float64}([1.0,-3.2,3.52,-1.408,0.0704,0.011264],5)

julia> print(ans)
Taylor{Float64}([1.0,-3.2,3.5200000000000005,-1.4080000000000004,0.07040000000000009,0.011264000000000012],5)

julia> (1+xT0)^xT0
Taylor{Float64}([1.0,0.0,1.0,-0.5,0.833333,-0.75],5)

julia> xTN(1)-yTN(1) # 1+x-(1+y)
TaylorN{Int64}( [HomogPol{Int64}([1,-1],1)], 5)

julia> xTN(1)*yTN(1)  # (1+x)*(1+y)
TaylorN{Int64}( [HomogPol{Int64}([1],0), HomogPol{Int64}([1,1],1), HomogPol{Int64}([0,1,0],2)], 5)

julia> 1/xTN(1)
TaylorN{Float64}( [HomogPol{Float64}([1.0],0), HomogPol{Float64}([-1.0,0.0],1), HomogPol{Float64}([1.0,0.0,0.0],2), HomogPol{Float64}([-1.0,0.0,0.0,0.0],3), HomogPol{Float64}([1.0,0.0,0.0,0.0,0.0],4), HomogPol{Float64}([-1.0,0.0,0.0,0.0,0.0,0.0],5)], 5)
```

Some elemental functions have been implemented using automatic differentiation techniques, i.e., by computing recursively their coefficients. Some examples include `exp`, `log`, `sqrt`, `sin`, `cos` and `tan`; more functions will be added in the future.

```julia
julia> exp(xT0)
Taylor{Float64}([1.0,1.0,0.5,0.166667,0.0416667,0.00833333],5)

julia> log(1-xT0)
Taylor{Float64}([0.0,-1.0,-0.5,-0.333333,-0.25,-0.2],5)

julia> sqrt(xT0)
ERROR: First non-vanishing Taylor coefficient must be an EVEN POWER
to expand SQRT around 0.

 in sqrt at /Users/benet/.julia/v0.3/TaylorSeries/src/utils_Taylor1.jl:265

julia> sqrt( xT(1.0) )   # identic to sqrt( 1+x0 )
Taylor{Float64}([1.0,0.5,-0.125,0.0625,-0.0390625,0.0273438],5)

julia> imag(exp(xTI)')
Taylor{Float64}([-0.0,-1.0,-0.0,0.166667,-0.0,-0.00833333],5)

julia> real(exp(Taylor([0.0,1im],17))) - cos(Taylor([0.0,1.0],17)) == 0.0
true

julia> exp(xTN(0.0)+yTN(0.0))
TaylorN{Float64}( [HomogPol{Float64}([1.0],0), HomogPol{Float64}([1.0,1.0],1), HomogPol{Float64}([0.5,1.0,0.5],2), HomogPol{Float64}([0.166667,0.5,0.5,0.166667],3), HomogPol{Float64}([0.0416667,0.166667,0.25,0.166667,0.0416667],4), HomogPol{Float64}([0.00833333,0.0416667,0.0833333,0.0833333,0.0416667,0.00833333],5)], 5)

julia> exp(xTN(0.0)+yTN(0.0)) - exp(xTN(0.0))*exp(yTN(0.0))  # should yield 0
TaylorN{Float64}( [HomogPol{Float64}([0.0],0)], 5)

julia> ans == 0.0
true
```

Differentiating and integrating is trivial for polynomial expansions in one variable. The last coefficient of the differential is simply set to zero; for the integral, an integration constant may be set to a different value (the default is zero). The order of the resulting polynomial is not changed. The n-th ($ n \ge 0 $) derivative is obtained with `deriv(a,n)`, where `a` is a Taylor series; default is n=1. Partial differentiation is also implemented for expansion in more than one variable. Note that integration is not yet implemented for `TaylorN`.

```julia
julia> diffTaylor(exp(xT0))
Taylor{Float64}([1.0,1.0,0.5,0.166667,0.0416667,0.0],5)

julia> integTaylor(exp(xT0))
Taylor{Float64}([0.0,1.0,0.5,0.166667,0.0416667,0.00833333],5)

julia> integTaylor( ans, 1.0)
Taylor{Float64}([1.0,1.0,0.5,0.166667,0.0416667,0.00833333],5)

julia> integTaylor( diffTaylor( exp(-xT0)), 1.0 ) == exp(-xT0)
true

julia> deriv( exp(xT(1.0))) == e
true

julia> deriv( exp(xT(1.0)), 5) == e
true

julia> diffTaylor( 3+xTN(0.0)+2yTN(0.0),1 )   # partial derivative with respect to 1st variable
TaylorN{Float64}( [HomogPol{Float64}([1.0],0)], 5)

julia> diffTaylor( 3+xTN(0.0)+2yTN(0.0),2 )   # partial derivative with respect to 2nd variable
TaylorN{Float64}( [HomogPol{Float64}([2.0],0)], 5)
```

For the evaluation of a Taylor series, we use Horner's rule in `evalTaylor(a::Taylor, dx::Number)`. Here, $ dx $ is the difference with respect to the point $ x_0 $ where the Taylor expansion is calculated, i.e., the series is evaluated at $ x = x_0 + dx $. Omitting $ dx $ corresponds to $ dx = 0 $. This function is generalized to admit `TaylorN` objects and vectors of numbers or `Taylor` objects; the length of these vectors must coincide with the number of independent variables.

```julia
julia> evalTaylor( exp( xT(1.0) )) - e ## exp(x) around x0=1 (order 5), evaluated there (dx=0)
0.0

julia> evalTaylor( exp(xT0), 1) - e    ## exp(x) around x0=0 (order 5), evaluated at x=1
-0.0016151617923783057

julia> evalTaylor( exp( Taylor(xT0, 17) ), 1) - e  ## exp(x) around x0=0 (order 17), evaluated at x=1
0.0

julia> evalTaylor( exp( Taylor([zero(BigFloat),one(BigFloat)],50) ), one(BigFloat) )
2.718281828459045235360287471352662497757247093699959574966967627723419298053556e+00 with 256 bits of precision

julia> convert(BigFloat,e) - ans
6.573322999985292556154129119543257102601105719980995128942636339920549561322098e-67 with 256 bits of precision

julia> evalTaylor(xTN(0.0)+yTN(0.0), [xT0, 2xT0])   # x+y, with x=x0, y=2x0
Taylor{Float64}([0.0,3.0,0.0,0.0,0.0,0.0],5)
```

In order to get a *nicer* displayed series, we have introduced `pretty_print`; we distinguish one-variable expansions with the variable $ x_{0} $, and more than one through the variables $ x_{1}, x_{2}, \dots $.

```julia
julia> pretty_print( xT0 )
5-order Taylor{Float64}:
  1.0 * x_{0}

julia> exp(xTI)
Taylor{Complex{Float64}}(Complex{Float64}[0.540302+0.841471im,0.540302+0.841471im,0.270151+0.420735im,0.0900504+0.140245im,0.0225126+0.0350613im,0.00450252+0.00701226im],5)

julia> pretty_print( ans )
5-order Taylor{Complex{Float64}}:
  ( 1.0 )  + ( 1.0 im ) * x_{0}  - ( 0.5 ) * x_{0}^2  - ( 0.16666666666666666 im ) * x_{0}^3  + ( 0.041666666666666664 ) * x_{0}^4  + ( 0.008333333333333333 im ) * x_{0}^5

julia> convert(Taylor{Rational{Int64}}, exp(xT0))
Taylor{Rational{Int64}}(Rational{Int64}[1//1,1//1,1//2,1//6,1//24,1//120],5)

julia> pretty_print(ans)
5-order Taylor{Rational{Int64}}:
 1//1 + 1//1 * x_{0} + 1//2 * x_{0}^2 + 1//6 * x_{0}^3 + 1//24 * x_{0}^4 + 1//120 * x_{0}^5

julia> pretty_print( sin(xTN(0.0)+yTN(0.0)) )
5-order TaylorN{Float64} in 2 variables:
 1.0 * x_{1} + 1.0 * x_{2} - 0.16666666666666666 * x_{1}^3 - 0.5 * x_{1}^2 * x_{2} - 0.5 * x_{1} * x_{2}^2 - 0.16666666666666666 * x_{2}^3 + 0.008333333333333333 * x_{1}^5 + 0.041666666666666664 * x_{1}^4 * x_{2} + 0.08333333333333333 * x_{1}^3 * x_{2}^2 + 0.08333333333333333 * x_{1}^2 * x_{2}^3 + 0.041666666666666664 * x_{1} * x_{2}^4 + 0.008333333333333333 * x_{2}^5
```

We have also implemented functions to compute gradients, Jacobians and Hessians. The following computes the gradients of $ f(x,y) = x^3 + 2x^2 y - 7 x + 2 $ and $ g(x,y) = y-x^4 $, using either `∇` or `TaylorSeries.gradient`; the results are of type `Array{TaylorN{T},1}`. In order to compute the Jacobian of a vector field evaluated on a point we use `jacobian`, and for the Hessian of a function we employ `hessian`.

```julia
julia> f(x,y) = x^3 + 2x^2 * y - 7x + 2
f (generic function with 1 method)

julia> g(x,y) = y - x^4
f (generic function with 1 method)

julia> f1 = f(xTN(0),yTN(0)); g1 = g(xTN(0),yTN(0))
TaylorN{Int64}( [HomogPol{Int64}([0,1],1), HomogPol{Int64}([-1,0,0,0,0],4)], 6)

julia> ∇(f1)
2-element Array{TaylorN{Int64},1}:
 TaylorN{Int64}( [HomogPol{Int64}([-7],0), HomogPol{Int64}([3,4,0],2)], 6)
                          TaylorN{Int64}( [HomogPol{Int64}([2,0,0],2)], 6)

julia> pretty_print(ans)
6-order TaylorN{Int64} in 2 variables:
 - 7 + 3 * x_{1}^2 + 4 * x_{1} * x_{2}

6-order TaylorN{Int64} in 2 variables:
 2 * x_{1}^2


julia> pretty_print( gradient( g1 ) )
6-order TaylorN{Int64} in 2 variables:
 - 4 * x_{1}^3

6-order TaylorN{Int64} in 2 variables:
 1


julia> jacobian([f1,g1], [2,1])
2x2 Array{Int64,2}:
  13  8
 -32  1

julia> hessian(f1-g1-2*f1*g1)
2x2 Array{Int64,2}:
  0  14
 14   0

```

Some concrete applications of the package will be found in the directory [examples][4]; for the time being we included the integration of [Kepler's problem using Taylor's method][5]).

Finally, it is worth pointing out the existing julia packages [Polynomial][2] and [PowerSeries][3] have similar functionality as `TaylorSeries` for one-variable expansions using somewhat different approaches.


#### Acknowledgments
This project began (using `python`) during a Masters' course in the postgraduate programs in Physics and in Mathematics at UNAM, during the second half of 2013. We thank the participants of the course for putting up with the half-baked material and contributing energy and ideas.

We acknowledge financial support from DGAPA-UNAM PAPIME grants PE-105911 and PE-107114. LB acknowledges support through a *Cátedra Moshinsky* (2013).

[1]: http://www.maia.ub.es/~alex/admcds/admcds.pdf
[2]: https://github.com/vtjnash/Polynomial.jl
[3]: https://github.com/jwmerrill/PowerSeries.jl
[4]: https://github.com/lbenet/TaylorSeries.jl/tree/master/examples/
[5]: http://nbviewer.ipython.org/github/lbenet/TaylorSeries.jl/blob/master/examples/1-KeplerProblem.ipynb
