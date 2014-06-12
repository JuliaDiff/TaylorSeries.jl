# TaylorSeries.jl

A julia package for Taylor expansions in one or more independent variables.

[![Build Status](https://travis-ci.org/lbenet/TaylorSeries.jl.svg?branch=master)](https://travis-ci.org/lbenet/TaylorSeries.jl)

#### Authors
- Luis Benet, Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)
- David P. Sanders, Facultad de Ciencias, Universidad Nacional Autónoma de México (UNAM)

Comments on the code and suggestions for improvements are welcome and appreciated.

### Basic of usage

This package computes Taylor expansions using automatic differentiation techniques (see e.g. Alex Haro's [nice paper][1] for a description of the techniques). The module implements the basic operations (+,-,*,/,^), some elementary functions (exp, log, sin, cos, tan), integration, differentiation and evaluation (Horner's rule) of the series on a point. The expansions work for one or more independent variables.

The package introduces two types, `Taylor` and `TaylorN`, for the expansions in one or more independent-variable expansions, respectively. Both are subtypes of `AbstractSeries{T,N}`, with `T` the type of the coefficients and `N` the number of independent variables, which is introduced as a subtype of `Number`.

`Taylor` is defined by one of the constructors

- `Taylor{T<:Number}(a::Taylor{T}, order::Integer)`
- `Taylor{T<:Number}(a::Array{T,1}, order::Integer)`
- `Taylor{T<:Number}(a::T, order::Integer)`

Thus, Taylor polynomials can be constructed from its coefficients (increasing with the order of the polynomial term), which may be of type `Taylor`, a vector or even a number (`<:Number`). The order of the polynomial can be omitted, which is then fixed from the length of the coefficients (otherwise it is the `max` among the length of the coefficients or the integer specified). 

For `TaylorN` there are two basic functions, `set_maxOrder()` and `set_numVars()`, which fix the overall constants defining the maximum order of the expansions and the number of variables, respectively. These are essential to construct the dictionaries that translate the position of the coefficients vector into the indexes of the multi-variable monomials. It is a good idea to fix this values before defining the `TaylorN` objects; if the number of variables is changed, all `TaylorN` objects must be redefined. (This could be improved in the future.)

`TaylorN` can be constructed by:

- `TaylorN{T<:Number}(x::TaylorN{T}, order::Int)`
- `TaylorN{T<:Number}(coeffs::Array{T,1}, order::Int)`

Again, the order can be omitted, in which case it is fixed by the value set by `set_maxOrder()`, which by default is 4. Similarly, the number of independent variables is by default 2. 

By default, the Taylor expansions are implemented around 0; if the expansion around a different value $x_0$ is needed, the trick is simply $x\to x+a$; see the definition of `xT(a)` below.

```julia
julia> using TaylorSeries

 help> Taylor   #?Taylor
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

 help> TaylorN   #?TaylorN
DataType   : TaylorN{T<:Number} (constructor with 6 methods)
  supertype: AbstractSeries{T<:Number,2}
  fields   : (:coeffs,:order,:numVars)

 help> TaylorSeries.AbstractSeries
DataType   : AbstractSeries{T<:Number,N}
  supertype: Number
  subtypes : {TaylorN{T<:Number},Taylor{T<:Number}}

julia> set_maxOrder(6)
INFO: MAXORDER is now 6; hash tables regenerated.
6

julia> get_numVars()
2

julia> TaylorN([1,2,3,4])
TaylorN{Int64}([1,2,3,4,0,0,0,0,0,0,0,0,0,0,0],4,2)

julia> xTN(a) = TaylorN([a,one(a)])
xTN (generic function with 1 method)

julia> yTN(a) = TaylorN([a,zero(a),one(a)])
yTN (generic function with 1 method)

julia> yTN(1.0pi)
TaylorN{Float64}([3.14159,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0  …  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],6,2)
```

The usual arithmetic operators (`+`, `-`, `*`, `/`, `^`, `==`) have been extended to work with `Taylor` and `TaylorN` type, including combinations of polynomials and numbers. 

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

julia> xTI^4/xT0^4
WARNING: Factorizing the polynomial.
The last k=4 Taylor coefficients ARE SET to 0.
Taylor{Complex{Float64}}(Complex{Float64}[1.0+0.0im,0.0+0.0im,0.0+0.0im,0.0+0.0im,0.0+0.0im,0.0+0.0im],5)

julia> (1-xT0)^3.2
Taylor{Float64}([1.0,-3.2,3.52,-1.408,0.0704,0.011264],5)

julia> (1+xT0)^xT0
Taylor{Float64}([1.0,0.0,1.0,-0.5,0.833333,-0.75],5)

julia> xTN(1)-yTN(1)
TaylorN{Int64}([0,1,-1,0,0,0,0,0,0,0  …  0,0,0,0,0,0,0,0,0,0],6,2)

julia> xTN(1)*yTN(1)   # (1+x)*(1+y)
TaylorN{Int64}([1,1,1,0,1,0,0,0,0,0  …  0,0,0,0,0,0,0,0,0,0],6,2)

julia> 1/xTN(1)
TaylorN{Float64}([1.0,-1.0,0.0,1.0,0.0,0.0,-1.0,0.0,0.0,0.0  …  0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0],6,2)

julia> println(ans)
TaylorN{Float64}([1.0,-1.0,0.0,1.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0],6,2)
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

 in sqrt at /Users/benet/.julia/v0.3/TaylorSeries/src/utils_Taylor1.jl:278

julia> sqrt( xT(1.0) )   # identic to sqrt( 1+x0 )
Taylor{Float64}([1.0,0.5,-0.125,0.0625,-0.0390625,0.0273438],5)

julia> imag(exp(xTI)')
Taylor{Float64}([-0.0,-1.0,-0.0,0.166667,-0.0,-0.00833333],5)

julia> real(exp(Taylor([0.0,1im],17))) - cos(Taylor([0.0,1.0],17)) == 0.0
true

julia> exp(xTN(0)+yTN(0))
TaylorN{Float64}([1.0,1.0,1.0,0.5,1.0,0.5,0.166667,0.5,0.5,0.166667  …  0.0833333,0.0416667,0.00833333,0.00138889,0.00833333,0.0208333,0.0277778,0.0208333,0.00833333,0.00138889],6,2)

julia> ans - exp(xTN(0))*exp(yTN(0)) == 0
true
```

Differentiating and integrating is trivial for polynomial expansions in one variable. The last coefficient of the differential is set to zero; for the integral, the integration constant may be set to a different value (the default is zero). The order of the resulting polynomial is not changed. The n-th ($n\ge 0$) derivative is obtained with `deriv(a,n)`, where `a` is a Taylor series; default is n=1. Partial differentiation is also implemented for expansion in more than one variable.

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

julia> diffTaylor( 3+xTN(0.0)+2yTN(0.0),1 )   # ∂(3+x+2y)/∂x
TaylorN{Float64}([1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0  …  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],6,2)

julia> diffTaylor( 3+xTN(0.0)+2yTN(0.0),2 )   # ∂(3+x+2y)/∂x
TaylorN{Float64}([2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0  …  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],6,2)
```

For the evaluation of a Taylor series, we use Horner's rule and `evalTaylor(a::Taylor, dx::Number)`. Here, $dx$ is the difference with respect to the point $x_0$ where the Taylor expansion is calculated, i.e., the series is evaluated at $x = x_0 + dx$. Omitting $dx$ corresponds to $dx=0$. This function is generalized to admit `TaylorN` objects and vectors (whose length is the number of independent variables) of numbers or `Taylor` objects.

```julia
julia> evalTaylor( exp( xT(1.0) )) - e ## exp(x) around x0=1 (order 5), evaluated there (dx=0)
0.0

julia> evalTaylor( exp(xT0), 1) - e    ## exp(x) around x0=0 (order 5), evaluated at x=1
-0.0016151617923783057

julia> evalTaylor( exp( Taylor(xT0, 17) ), 1) - e  ## exp(x) around x0=0 (order 17), evaluated at x=1
0.0

julia> evalTaylor( exp( Taylor([0,BigFloat("1.0")],50) ), BigFloat("1.0") )
2.718281828459045235360287471352662497757247093699959574966967627723419298053556e+00 with 256 bits of precision

julia> convert(BigFloat,e) - ans
6.573322999985292556154129119543257102601105719980995128942636339920549561322098e-67 with 256 bits of precision

julia> evalTaylor(xTN(0.0)+yTN(0.0), [xT0, 2xT0])   # x+y, with x=x0, y=2x0
Taylor{Float64}([0.0,3.0,0.0,0.0,0.0,0.0],5)
```

Finally, in order to get a *compact* and more readable display, we have introduced `pretty_print`; we distinguish one-variable expansions with the variable $x0$, and more than one through the variables $x1, x2, \dots$.

```julia
julia> pretty_print( xT0 )
5-order Taylor{Float64}:
  1.0 * x0 

julia> exp(xTI)
Taylor{Complex{Float64}}(Complex{Float64}[0.540302+0.841471im,0.540302+0.841471im,0.270151+0.420735im,0.0900504+0.140245im,0.0225126+0.0350613im,0.00450252+0.00701226im],5)

julia> pretty_print( ans )
5-order Taylor{Complex{Float64}}:
    ( 0.5403023058681398 + 0.8414709848078965 im )  + ( 0.5403023058681398 + 0.8414709848078965 im ) * x0  + ( 0.2701511529340699 + 0.42073549240394825 im ) * x0^2  + ( 0.09005038431135663 + 0.1402451641346494 im ) * x0^3  + ( 0.022512596077839158 + 0.03506129103366235 im ) * x0^4  + ( 0.004502519215567832 + 0.00701225820673247 im ) * x0^5 

julia> convert(Taylor{Rational{Int64}}, exp(xT0))
Taylor{Rational{Int64}}(Rational{Int64}[1//1,1//1,1//2,1//6,1//24,1//120],5)

julia> pretty_print(ans)
5-order Taylor{Rational{Int64}}:
  1//1 + 1//1 * x0 + 1//2 * x0^2 + 1//6 * x0^3 + 1//24 * x0^4 + 1//120 * x0^5 

julia> pretty_print( sin(xTN(0.0)+yTN(0.0)) )
6-order TaylorN{Float64} in 2 variables:
  1.0 * x1 + 1.0 * x2 - 0.16666666666666666 * x1^3 - 0.5 * x1^2 * x2 - 0.5 * x1 * x2^2 - 0.16666666666666666 * x2^3 + 0.008333333333333333 * x1^5 + 0.041666666666666664 * x1^4 * x2 + 0.08333333333333333 * x1^3 * x2^2 + 0.08333333333333333 * x1^2 * x2^3 + 0.041666666666666664 * x1 * x2^4 + 0.008333333333333333 * x2^5
```

Some concrete applications of the package can be found in the directory [examples/][4].

Finally, it is worth pointing out the existing julia packages [Polynomial][2] and [PowerSeries][3], which have similar functionality as `TaylorSeries` for one-variable expansions, but using somewhat different approaches.


#### Acknowledgments 
This project began (using `python`) during a Masters' course in the postgraduate programs in Physics and in Mathematics at UNAM, during the second half of 2013. We thank the participants of the course for putting up with the half-baked material and contributing energy and ideas.

We acknowledge financial support from DGAPA-UNAM PAPIME grants PE-105911 and PE-107114. LB acknowledges support through a *Cátedra Moshinsky* (2013).

[1]: http://www.maia.ub.es/~alex/admcds/admcds.pdf
[2]: https://github.com/vtjnash/Polynomial.jl
[3]: https://github.com/jwmerrill/PowerSeries.jl
[4]: https://github.com/lbenet/TaylorSeries.jl/tree/master/examples/

