# TaylorSeries.jl

A julia package for Taylor expansions in one independent variable.


### An overlook

This package computes Taylor expansions (around 0) using automatic differentiation techniques (see e.g. Alex Haro's [nice paper][1]). The module implements the basic operations (+,-,*,/,^), some elementary functions (exp, log, sin, cos, tan), integration, differentiation and (a Horner's rule) evaluation on a point. 

The basic constructors for the Taylor type are:

- `Taylor{T<:Number}(a::Taylor{T}, order::Integer)`
- `Taylor{T<:Number}(a::Array{T,1}, order::Integer)` 
- `Taylor{T<:Number}(a::T, order::Integer)`

Taylor polynomials can be constructed from its coefficients (increasing with the order of the polynomial), which may be of type Taylor, a vector or even a number (`<:Number`). The order of the polynomial can be omitted, which is then guessed from the length of the coefficients (otherwise it is the `max` among the length of the coefficients or the integer specified). 

By default, the Taylor expansions are implemented around 0; if the expansion around a different value $x_0$ is needed, the trick is $x\to x+a$; see the definition of `x(x0)` below. 

```julia
julia> using TaylorSeries

julia> ?Taylor
Loading help data...
DataType   : Taylor{T<:Number} (constructor with 6 methods)
  supertype: Any
  fields   : (:coeffs,:order)

julia> Taylor([0,1,2])
Taylor{Int64}([0,1,2],2)

julia> Taylor(0.0)
Taylor{Float64}([0.0],0)

julia> Taylor([0.0,1im])
Taylor{Complex{Float64}}(Complex{Float64}[0.0+0.0im,0.0+1.0im],1)

julia> x(a) = Taylor([a, 1.0], 5) # expansion (order 5) around an arbitrary value a
x (generic function with 1 method)

julia> x(1.0)
Taylor{Float64}([1.0,1.0,0.0,0.0,0.0,0.0],5)

julia> x0 = x(0.0) # abbreviation for expansions around 0
Taylor{Float64}([0.0,1.0,0.0,0.0,0.0,0.0],5)

```

The usual arithmetic operators (`+`, `-`, `*`, `/`, `^`, `==`) have been extended to work with the Taylor type, including combinations of polynomials and scalars. 

```julia
julia> x0*(3x0+2.5)
Taylor{Float64}([0.0,2.5,3.0,0.0,0.0,0.0],5)

julia> 1/(1-x0)
Taylor{Float64}([1.0,1.0,1.0,1.0,1.0,1.0],5)

julia> x0*(x0^2-4)/(x0+2)
Taylor{Float64}([-0.0,-2.0,1.0,0.0,0.0,0.0],5)

julia> Taylor([1],5) == 1.0
true

julia> xI = Taylor([0,im],5)
Taylor{Complex{Int64}}(Complex{Int64}[0+0im,0+1im,0+0im,0+0im,0+0im,0+0im],5)

julia> xI'
Taylor{Complex{Int64}}(Complex{Int64}[0+0im,0-1im,0+0im,0+0im,0+0im,0+0im],5)

julia> xI^4/x0^4
WARNING: Applying L'Hôpital. The last k=4 Taylor coefficients ARE SET to 0.
Taylor{Complex{Float64}}(Complex{Float64}[1.0+0.0im,0.0+0.0im,0.0+0.0im,0.0+0.0im,0.0+0.0im,0.0+0.0im],5)

julia> (1-x0)^3.2
Taylor{Float64}([1.0,-3.2,3.52,-1.408,0.0704,0.011264],5)

julia> (1+x0)^x0
Taylor{Float64}([1.0,0.0,1.0,-0.5,0.833333,-0.75],5)

```

Some elemental functions have been implemented using automatic differentiation techniques, i.e., by computing recursively their coefficients. These are: `exp`, `log`, `square`, `sqrt`, `sin`, `cos` and `tan`; more functions will be added in the future.
```julia
julia> exp( -x0 )
Taylor{Float64}([1.0,-1.0,0.5,-0.166667,0.0416667,-0.00833333],5)

julia> exp( -x(2.0) )
Taylor{Float64}([0.135335,-0.135335,0.0676676,-0.0225559,0.00563897,-0.00112779],5)

julia> log(1-x0)
Taylor{Float64}([0.0,-1.0,-0.5,-0.333333,-0.25,-0.2],5)

julia> sqrt(x0)
ERROR: First non-vanishing Taylor coefficient must be an EVEN POWER
to expand SQRT around 0.

 in sqrt at /Users/benet/.julia/v0.3/TaylorSeries/src/TaylorSeries.jl:281

julia> sqrt( x(1.0) )
Taylor{Float64}([1.0,0.5,-0.125,0.0625,-0.0390625,0.0273438],5)

julia> sqrt( 1+x0 )
Taylor{Float64}([1.0,0.5,-0.125,0.0625,-0.0390625,0.0273438],5)

julia> imag(exp(xI)')
Taylor{Float64}([-0.0,-1.0,-0.0,0.166667,-0.0,-0.00833333],5)

julia> real(exp(xI)) - cos(x0) == 0
true

```

Differentiating and integrating is trivial for polynomial expansions. The last coefficient of the differential is set to zero; for the integral, the integration constant may be set to a different value (the default is zero). The order of the resulting polynomial is not changed. The n-th ($n\ge 0$) derivative is obtained with `deriv(a,n)`, where `a` is a Taylor series; default is n=1.

```julia
julia> diffTaylor(exp(x0))
Taylor{Float64}([1.0,1.0,0.5,0.166667,0.0416667,0.0],5)

julia> integTaylor(exp(x0))
Taylor{Float64}([0.0,1.0,0.5,0.166667,0.0416667,0.00833333],5)

julia> integTaylor( exp(x0), 1.0)
Taylor{Float64}([1.0,1.0,0.5,0.166667,0.0416667,0.00833333],5)

julia> integTaylor( diffTaylor( exp(-x0)), 1.0 ) == exp(-x0)
true

julia> deriv( exp(x(1.0))) == e
true

julia> deriv( exp(x(1.0)), 5) == e
true
```

For the evaluation of a Taylor series, we use Horner's rule and `evalTaylor(a::Taylor, dx::Number)`. Here, $dx$ is the difference with respect to the point $x_0$ where the Taylor expansion was calculated, i.e., the series is evaluated at $x = x_0 + dx$. Omitting $dx$ corresponds to $dx=0$.

```julia
julia> evalTaylor( exp( x(1.0) )) - e ## exp(x) around x0=1 (order 5), evaluated there (dx=0)
0.0

julia> evalTaylor( exp(x0), 1) - e    ## exp(x) around x0=0 (order 5), evaluated at x=1
-0.0016151617923783057

julia> evalTaylor( exp( Taylor(x0, 17) ), 1) - e  ## exp(x) around x0=0 (order 17), evaluated at x=1
0.0

julia> evalTaylor( exp( Taylor([0,BigFloat("1.0")],50) ), BigFloat("1.0") )
2.718281828459045235360287471352662497757247093699959574966967627723419298053556e+00 with 256 bits of precision

julia> convert(BigFloat,e) - ans
6.573322999985292556154129119543257102601105719980995128942636339920549561322098e-67 with 256 bits of precision

```

In order to get a *compact* and more readable display, we have modified `showcompact`.
```julia
julia> showcompact( x0 )
5-order Taylor{Float64}:
  1.0 * x 

julia> exp(xI)
Taylor{Complex{Float64}}(Complex{Float64}[1.0+0.0im,0.0+1.0im,-0.5+0.0im,0.0-0.166667im,0.0416667+0.0im,0.0+0.00833333im],5)

julia> showcompact( exp(xI) )
5-order Taylor{Complex{Float64}}
 ( 1.0 )  + ( 1.0 im ) * x  - ( 0.5 ) * x^2  - ( 0.16666666666666666 im ) * x^3  + ( 0.041666666666666664 ) * x^4  + ( 0.008333333333333333 im ) * x^5 ) 

julia> print( exp(xI) )
Taylor{Complex{Float64}}(Complex{Float64}[1.0 + 0.0im,0.0 + 1.0im,-0.5 + 0.0im,0.0 - 0.16666666666666666im,0.041666666666666664 + 0.0im,0.0 + 0.008333333333333333im],5)

julia> convert(Taylor{Rational{Int64}}, exp(x0))
Taylor{Rational{Int64}}(Rational{Int64}[1//1,1//1,1//2,1//6,1//24,1//120],5)

julia> showcompact(ans)
5-order Taylor{Rational{Int64}}:
  1//1  + 1//1 * x + 1//2 * x^2 + 1//6 * x^3 + 1//24 * x^4 + 1//120 * x^5 

```

It is worth pointing out the existing julia packages [Polynomial][2] and [PowerSeries][3], which have similar functionality as this one, but somewhat different approaches.

Comments on the code are appreciated.

[1]: http://www.maia.ub.es/~alex/admcds/admcds.pdf
[2]: https://github.com/vtjnash/Polynomial.jl
[3]: https://github.com/jwmerrill/PowerSeries.jl

