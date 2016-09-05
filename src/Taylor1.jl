# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#


## Constructors ##
doc"""
    Taylor1{T<:Number} <: Number

DataType for polynomial expansions in one independent variable.

**Fields:**

- `coeffs :: Array{T,1}` Expansion coefficients; the $i$-th
component is the coefficient of degree $i-1$ of the expansion.
- `order  :: Int64` Maximum order (degree) of the polynomial.
"""
immutable Taylor1{T<:Number} <: Number
    coeffs :: Array{T,1}
    order :: Int

    ## Inner constructor ##
    function Taylor1(coeffs::Array{T,1}, order::Int)
        lencoef = length(coeffs)
        order = max(order, lencoef-1)
        if order == lencoef-1
            return new(coeffs, order)
        else
            resize!(coeffs, order+1)
            for i = lencoef+1:order+1
                coeffs[i] = zero(T)
            end
            return new(coeffs, order)
        end
    end
end

## Outer constructors ##
Taylor1{T<:Number}(x::Taylor1{T}, order::Int) = Taylor1{T}(x.coeffs, order)
Taylor1{T<:Number}(x::Taylor1{T}) = x
Taylor1{T<:Number}(coeffs::Array{T,1}, order::Int) = Taylor1{T}(coeffs, order)
Taylor1{T<:Number}(coeffs::Array{T,1}) = Taylor1{T}(coeffs, length(coeffs)-1)
Taylor1{T<:Number}(x::T, order::Int) = Taylor1{T}([x], order)
# Taylor1{T<:Number}(x::T) = Taylor1{T}([x], 0)

# Shortcut to define Taylor1 independent variables
"""
    Taylor1(T, [order=1])
    Taylor1([order=1])

Shortcut to define the independent variable of a `Taylor1{T}` polynomial of
given `order`. If the type `T` is ommitted, `Float64` is assumend.
"""
Taylor1{T<:Number}(::Type{T}, order::Int=1) = Taylor1{T}( [zero(T), one(T)], order)
Taylor1(order::Int=1) = Taylor1(Float64, order)

## get_coeff ##
"""
    get_coeff(a, n)

Return the coefficient of order `n::Int` of a `a::Taylor1` polynomial.
"""
get_coeff(a::Taylor1, n::Int) = (@assert 0 <= n <= a.order+1;
    return a.coeffs[n+1])

## Type, length ##
eltype{T<:Number}(::Taylor1{T}) = T
length{T<:Number}(a::Taylor1{T}) = a.order

## Conversion and promotion rules ##
convert{T<:Number}(::Type{Taylor1{T}}, a::Taylor1) =
    Taylor1(convert(Array{T,1}, a.coeffs), a.order)
function convert{T<:Integer, S<:AbstractFloat}(::Type{Taylor1{Rational{T}}},
    a::Taylor1{S})
    v = Array(Rational{T}, length(a.coeffs))
    for i in eachindex(v)
        # v[i] = convert(Rational{T}, a.coeffs[i])
        v[i] = rationalize(a.coeffs[i], tol=eps(one(S)))
    end
    Taylor1(v)
end
convert{T<:Number}(::Type{Taylor1{T}}, b::Array{T,1}) = Taylor1(b, length(b)-1)
convert{T<:Number, S<:Number}(::Type{Taylor1{T}}, b::Array{S,1}) =
    Taylor1(convert(Array{T,1},b))
convert{T<:Number, S<:Number}(::Type{Taylor1{T}}, b::S) = Taylor1([convert(T,b)], 0)
convert{T<:Number}(::Type{Taylor1{T}}, b::T) = Taylor1([b], 0)

promote_rule{T<:Number}(::Type{Taylor1{T}}, ::Type{Taylor1{T}}) = Taylor1{T}
promote_rule{T<:Number, S<:Number}(::Type{Taylor1{T}}, ::Type{Taylor1{S}}) =
    Taylor1{promote_type(T, S)}
promote_rule{T<:Number}(::Type{Taylor1{T}}, ::Type{Array{T,1}}) = Taylor1{T}
promote_rule{T<:Number, S<:Number}(::Type{Taylor1{T}}, ::Type{Array{S,1}}) =
    Taylor1{promote_type(T, S)}
promote_rule{T<:Number}(::Type{Taylor1{T}}, ::Type{T}) = Taylor1{T}
promote_rule{T<:Number, S<:Number}(::Type{Taylor1{T}}, ::Type{S}) =
    Taylor1{promote_type(T, S)}

## Auxiliary function ##
function firstnonzero{T<:Number}(ac::Vector{T})
    nonzero::Int = length(ac)
    for i in eachindex(ac)
        if ac[i] != zero(T)
            nonzero = i-1
            break
        end
    end
    nonzero
end
firstnonzero{T<:Number}(a::Taylor1{T}) = firstnonzero(a.coeffs)

function fixshape{T<:Number}(a::Taylor1{T}, b::Taylor1{T})
    if a.order == b.order
        return a, b
    elseif a.order < b.order
        return Taylor1(a, b.order), b
    end
    return a, Taylor1(b, a.order)
end

function fixshape{T<:Number, S<:Number}(a::Taylor1{T}, b::Taylor1{S})
    fixshape(promote(a,b)...)
end

## real, imag, conj and ctranspose ##
for f in (:real, :imag, :conj)
    @eval ($f){T<:Number}(a::Taylor1{T}) = Taylor1(($f)(a.coeffs), a.order)
end
ctranspose{T<:Number}(a::Taylor1{T}) = conj(a)

## zero and one ##
zero{T<:Number}(a::Taylor1{T}) = Taylor1(zero(T), a.order)
one{T<:Number}(a::Taylor1{T}) = Taylor1(one(T), a.order)

## Equality ##
function ==(a::Taylor1, b::Taylor1)
    a, b = fixshape(a, b)
    return a.coeffs == b.coeffs
end

# Tests `isinf` and `isnan` for all polynomial coefficients
for f in (:isinf, :isnan)
    @eval begin
        function ($f)(a::Taylor1)
            test = false
            for i in eachindex(a.coeffs)
                @inbounds test = ($f)(a.coeffs[i])
                test && break
            end
            test
        end
    end
end

## Addition and substraction ##
for f in (:+, :-)
    @eval begin
        function ($f)(a::Taylor1, b::Taylor1)
            a, b = fixshape(a, b)
            v = similar(a.coeffs)
            @simd for i in eachindex(a.coeffs)
                @inbounds v[i] = ($f)(a.coeffs[i], b.coeffs[i])
            end
            return Taylor1(v, a.order)
        end
        function ($f)(a::Taylor1)
            v = similar(a.coeffs)
            @simd for i in eachindex(v)
                @inbounds v[i] = ($f)(a.coeffs[i])
            end
            return Taylor1(v, a.order)
        end
        function ($f)(a::Taylor1, b::Union{Real,Complex})
            @inbounds aux = ($f)(a.coeffs[1], b)
            v = Array(typeof(aux), length(a.coeffs))
            @simd for i in eachindex(v)
                @inbounds v[i] = a.coeffs[i]
            end
            @inbounds v[1] = aux
            Taylor1(v, a.order)
        end
        function ($f)(a::Union{Real,Complex}, b::Taylor1)
            @inbounds aux = ($f)(a, b.coeffs[1])
            v = Array(typeof(aux), length(b.coeffs))
            @simd for i in eachindex(v)
                @inbounds v[i] = ($f)(b.coeffs[i])
            end
            @inbounds v[1] = aux
            Taylor1(v, b.order)
        end
    end
end

## Multiplication ##
*(a::Bool, b::Taylor1) = *(promote(a,b)...)
*(a::Taylor1, b::Bool) = b*a
function *(a::Union{Real,Complex}, b::Taylor1)
    @inbounds aux = a * b.coeffs[1]
    v = Array(typeof(aux), length(b.coeffs))
    @simd for i in eachindex(v)
        @inbounds v[i] = a * b.coeffs[i]
    end
    Taylor1(v, b.order)
end
*(a::Taylor1, b::Union{Real,Complex}) = b * a
doc"""
```
*(a, b)
```

Return the Taylor expansion of $a \cdot b$, of order `max(a.order,b.order)`, for
`a::Taylor1`, `b::Taylor1` polynomials.

For details on making the Taylor expansion, see [`TaylorSeries.mulHomogCoef`](@ref).
"""
function *(a::Taylor1, b::Taylor1)
    a, b = fixshape(a, b)
    coeffs = similar(a.coeffs)
    @inbounds coeffs[1] = a.coeffs[1] * b.coeffs[1]
    @inbounds for k = 1:a.order
        coeffs[k+1] = mulHomogCoef(k, a.coeffs, b.coeffs)
    end
    Taylor1(coeffs, a.order)
end

# Homogeneous coefficient for the multiplication
doc"""
    mulHomogCoef(kcoef, ac, bc)

Compute the `k`-th expansion coefficient of $c = a\cdot b$ given by

\begin{equation*}
c_k = \sum_{j=0}^k a_j b_{k-j},
\end{equation*}

with $a$ and $b$ `Taylor1` polynomials.

Inputs are the `kcoef`-th coefficient, and the vectors of the expansion coefficients
`ac` and `bc`, corresponding respectively to `a` and `b`.
"""
function mulHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, bc::Array{T,1})
    kcoef == 0 && return ac[1] * bc[1]
    coefhomog = zero(T)
    @inbounds for i = 0:kcoef
        coefhomog += ac[i+1] * bc[kcoef-i+1]
    end
    coefhomog
end

## Division ##
/{T<:Real}(a::Taylor1, b::T) = a * inv(b)
/{T<:Complex}(a::Taylor1, b::T) = a * inv(b)
doc"""
```
/(a, b)
```

Return the Taylor expansion of $a/b$, of order `max(a.order,b.order)`, for
`a::Taylor1`, `b::Taylor1` polynomials.

For details on making the Taylor expansion, see [`TaylorSeries.divHomogCoef`](@ref).
"""
function /(a::Taylor1, b::Taylor1)
    a, b = fixshape(a, b)
    # order and coefficient of first factorized term
    orddivfact, cdivfact = divfactorization(a, b)
    T = typeof(cdivfact)
    v1 = convert(Array{T,1}, a.coeffs)
    v2 = convert(Array{T,1}, b.coeffs)
    coeffs = zeros(T, a.order+1)
    @inbounds coeffs[1] = cdivfact
    @inbounds for k = orddivfact+1:a.order
        coeffs[k-orddivfact+1] = divHomogCoef(k, v1, v2, coeffs, orddivfact)
    end
    Taylor1(coeffs, a.order)
end

function divfactorization(a1::Taylor1, b1::Taylor1)
    # order of first factorized term; a1 and b1 assumed to be of the same order
    a1nz = firstnonzero(a1)
    b1nz = firstnonzero(b1)
    orddivfact = min(a1nz, b1nz)
    if orddivfact > a1.order
        orddivfact = a1.order
    end
    cdivfact = a1.coeffs[orddivfact+1] / b1.coeffs[orddivfact+1]

    # Is the polynomial factorizable?
    if isinf(cdivfact) || isnan(cdivfact)
        throw(ArgumentError(
        """Division does not define a Taylor1 polynomial
        or its first non-zero coefficient is Inf/NaN.
        Order k=$(orddivfact) => coeff[$(orddivfact+1)]=$(cdivfact)."""))
    end

    return orddivfact, cdivfact
end

# Homogeneous coefficient for the division
doc"""
    divHomogCoef(kcoef, ac, bc, coeffs, ordfact)

Compute the `k-th` expansion coefficient of $c = a / b$ given by

\begin{equation*}
c_k =  \frac{1}{b_0} (a_k - \sum_{j=0}^{k-1} c_j b_{k-j}),
\end{equation*}

with $a$ and $b$ `Taylor1` polynomials.

Inputs are the `kcoef`-th coefficient, the vectors of the expansion coefficients
`ac` and `bc`, corresponding respectively to `a` and `b`, the
already calculated expansion coefficients `coeffs` of `c`, and `ordfact`
which is the order of the factorized term of the denominator,
whenever `b_0` is zero.
"""
function divHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, bc::Array{T,1},
    coeffs::Array{T,1}, ordfact::Int)
    #
    @inbounds kcoef == ordfact && return ac[ordfact+1] / bc[ordfact+1]
    coefhomog = mulHomogCoef(kcoef, coeffs, bc)
    @inbounds coefhomog = (ac[kcoef+1]-coefhomog) / bc[ordfact+1]
    coefhomog
end

## Division functions: rem and mod
for op in (:mod, :rem)
    @eval begin
        function ($op){T<:Real, S<:Real}(a::Taylor1{T}, x::S)
            coeffs = copy(a.coeffs)
            @inbounds coeffs[1] = ($op)(a.coeffs[1], x)
            return Taylor1(coeffs, a.order)
        end
    end
end

function mod2pi{T<:Real}(a::Taylor1{T})
    coeffs = copy(a.coeffs)
    @inbounds coeffs[1] = mod2pi( a.coeffs[1] )
    return Taylor1( coeffs, a.order)
end

## abs function ##
"""
    abs(a)

Return `a` or `-a` depending on the 0-th order coefficient of the
`Taylor1` polynomial `a`.
If `a.coeffs[1]` is zero, an `ArgumentError` is thrown.
"""
function abs{T<:Real}(a::Taylor1{T})
    if a.coeffs[1] > zero(T)
        return a
    elseif a.coeffs[1] < zero(T)
        return -a
    else
        throw(ArgumentError(
        """The 0th order Taylor1 coefficient must be non-zero
        (abs(x) is not differentiable at x=0)."""))
    end
end

## Int power ##
doc"""
```
^(a, x)
```

Return the Taylor expansion of $a^x$ for `a::Taylor1` polynomial and `x::Number`.
If `x` is non integer and the 0-th order coefficient is zero, an
`ArgumentError` is thrown.
"""
function ^{T<:Number}(a::Taylor1{T}, n::Integer)
    n == 0 && return one(a)
    n == 1 && return a
    n == 2 && return square(a)
    n < 0 && return inv( a^(-n) )
    return power_by_squaring(a, n)
end

function ^{T<:Integer}(a::Taylor1{T}, n::Integer)
    n == 0 && return one(a)
    n == 1 && return a
    n == 2 && return square(a)
    n < 0 && throw(DomainError())
    return power_by_squaring(a, n)
end

## power_by_squaring; slightly modified from base/intfuncs.jl
## Licensed under MIT "Expat"
function power_by_squaring(x::Taylor1, p::Integer)
    p == 1 && return copy(x)
    p == 0 && return one(x)
    p == 2 && return square(x)
    t = trailing_zeros(p) + 1
    p >>= t

    while (t -= 1) > 0
        x = square(x)
    end

    y = x
    while p > 0
        t = trailing_zeros(p) + 1
        p >>= t
        while (t -= 1) >= 0
            x = square(x)
        end
        y *= x
    end

    return y
end

## Rational power ##
^{T<:Integer}(a::Taylor1,x::Rational{T}) = a^(x.num/x.den)

^(a::Taylor1, b::Taylor1) = exp( b*log(a) )

## Real power ##
function ^{S<:Real}(a::Taylor1, x::S)
    x == zero(x) && return one(a)
    x == one(x)/2 && return sqrt(a)
    isinteger(x) && return a^round(Int,x)

    # First non-zero coefficient
    l0nz = firstnonzero(a)
    l0nz > a.order && return zero(a)

    # The first non-zero coefficient of the result; must be integer
    lnull = x*l0nz
    !isinteger(lnull) &&
        throw(ArgumentError(
        """The 0th order Taylor1 coefficient must be non-zero
        to raise the Taylor1 polynomial to a non-integer exponent."""))

    # Reaching this point, it is possible to implement the power of the Taylor1
    # polynomial. The last l0nz coefficients are set to zero.
    lnull = trunc(Int,lnull)
    #l0nz > 0 && warn("The last k=$(l0nz) Taylor1 coefficients ARE SET to 0.")
    @inbounds aux = (a.coeffs[l0nz+1])^x
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    coeffs = zeros(T, a.order+1)
    @inbounds coeffs[lnull+1] = aux
    k0 = lnull+l0nz
    @inbounds for k = k0+1:a.order
        coeffs[k-l0nz+1] = powHomogCoef(k, v, x, coeffs, l0nz)
    end

    Taylor1(coeffs,a.order)
end
^{T<:Complex}(a::Taylor1, x::T) = exp( x*log(a) )

# Homogeneous coefficients for real power
doc"""
    powHomogCoef(kcoef, ac, x, coeffs, knull)

Compute the `k-th` expansion coefficient of $c = a^x$, given by

\begin{equation*}
c_k = \frac{1}{k a_0} \sum_{j=0}^{k-1} ((k-j)x -j) - j)a_{k-j} c_j,
\end{equation*}

with $a$ a `Taylor1` polynomial, and `x` a number.

Inputs are the `kcoef`-th coefficient, the vector of the expansion coefficients
`ac` of `a`, the exponent `x`, the already calculated expansion
coefficients `coeffs` of `c`, and `knull`
which is the order of the first non-zero coefficient of `a`.
"""
function powHomogCoef{T<:Number, S<:Real}(kcoef::Int, ac::Array{T,1}, x::S,
    coeffs::Array{T,1}, knull::Int)

    kcoef == knull && return (ac[knull+1])^x
    coefhomog = zero(T)
    for i = 0:kcoef-knull-1
        aux = x*(kcoef-i)-i
        @inbounds coefhomog += aux*ac[kcoef-i+1]*coeffs[i+1]
    end
    aux = kcoef - knull*(x+1)
    @inbounds coefhomog = coefhomog / (aux*ac[knull+1])

    coefhomog
end

## Square ##
function square(a::Taylor1)
    coeffs = Array(eltype(a), a.order+1)
    coeffs[1] = a.coeffs[1]^2
    @inbounds for k = 1:a.order
        coeffs[k+1] = squareHomogCoef(k, a.coeffs)
    end
    Taylor1(coeffs,a.order)
end

# Homogeneous coefficients for square
doc"""
    squareHomogCoef(kcoef, ac)

Compute the `k-th` expansion coefficient of $c = a^2$, given by

\begin{eqnarray*}
c_k &=& 2 \sum_{j=0}^{(k-1)/2} a_{k-j} a_j, \text{ if $k$ is odd, or }  \\\\
c_k &=& 2 \sum_{j=0}^{(k-2)/2} a_{k-j} a_j + (a_{k/2})^2, \text{ if $k$ is even, }
\end{eqnarray*}

with $a$ a `Taylor1` polynomial.

Inputs are the `kcoef`-th coefficient and the vector of the expansion coefficients
`ac` of `a`.
"""
function squareHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1})
    kcoef == 0 && return ac[1]^2
    coefhomog = zero(T)
    kodd = kcoef%2
    kend = div(kcoef - 2 + kodd, 2)
    @inbounds for i = 0:kend
        coefhomog += ac[i+1]*ac[kcoef-i+1]
    end
    coefhomog = convert(T,2) * coefhomog
    if kodd == 0
        @inbounds coefhomog += ac[div(kcoef,2)+1]^2
    end
    coefhomog
end

## Square root ##
doc"""
    sqrt(a)

Return  the Taylor expansion of $\sqrt(a)$, of order `a.order`,
for `a::Taylor1` polynomial.
If the first non-vanishing coefficient of `a` is an odd power, and
`ArgumentError` is thrown.

For details on making the Taylor expansion, see [`TaylorSeries.sqrtHomogCoef`](@ref).
"""
function sqrt(a::Taylor1)
    # First non-zero coefficient
    l0nz = firstnonzero(a)
    if l0nz > a.order
        return zero(a)
    elseif l0nz%2 == 1 # l0nz must be pair
        throw(ArgumentError(
        """First non-vanishing Taylor1 coefficient must correspond
        to an **even power** in order to expand `sqrt` around 0."""))
    end

    # Reaching this point, it is possible to implement the sqrt of the Taylor1 polynomial.
    # The last l0nz coefficients are set to zero.
    ##l0nz > 0 && warn("The last k=$(l0nz) Taylor1 coefficients ARE SET to 0.")
    lnull = div(l0nz, 2)
    @inbounds aux = sqrt(a.coeffs[l0nz+1])
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    coeffs = zeros(T, a.order+1)
    @inbounds coeffs[lnull+1] = aux
    @inbounds for k = lnull+1:a.order-l0nz
        coeffs[k+1] = sqrtHomogCoef(k, v, coeffs, lnull)
    end
    Taylor1(coeffs, a.order)
end

# Homogeneous coefficients for the square-root
doc"""
    sqrtHomogCoef(kcoef, ac, coeffs, knull)

Compute the `k-th` expansion coefficient of $c = \sqrt(a)$, given by

\begin{eqnarray*}
c_k &=& \frac{1}{2 c_0} ( a_k - 2 \sum_{j=0}^{(k-1)/2} c_{k-j}c_j),
\text{ if $k$ is odd, or } \\\\
c_k &=& \frac{1}{2 c_0} ( a_k - 2 \sum_{j=0}^{(k-2)/2} c_{k-j}c_j) - (c_{k/2})^2,
\text{ if $k$ is even,}
\end{eqnarray*}

with $a$ a `Taylor1` polynomial.

Inputs are the `kcoef`-th coefficient, the vector of the expansion coefficients
`ac` of `a`, the already calculated expansion coefficients `coeffs` of `c`,
and `knull`, which is half of the order of the first non-zero coefficient of `a`.
"""
function sqrtHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, coeffs::Array{T,1}, knull::Int)
    kcoef == knull && return sqrt(ac[2*knull+1])
    coefhomog = zero(T)
    kodd = (kcoef - knull)%2
    kend = div(kcoef - knull - 2 + kodd, 2)
    two = convert(T,2)
    @inbounds for i = knull+1:knull+kend
        coefhomog += coeffs[i+1]*coeffs[kcoef+knull-i+1]
    end
    @inbounds aux = ac[kcoef+knull+1] - two*coefhomog

    if kodd == 0
        @inbounds aux = aux - (coeffs[kend+knull+2])^2
    end
    @inbounds coefhomog = aux / (two*coeffs[knull+1])

    coefhomog
end

## Exp ## 
doc"""
    exp(a)

Return the Taylor expansion of $e^a$, of order `a.order`, for
`a::Taylor1` polynomial.

For details on making the Taylor expansion, see [`TaylorSeries.expHomogCoef`](@ref).
"""
function exp(a::Taylor1)
    @inbounds aux = exp( a.coeffs[1] )
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    coeffs = similar(v)
    @inbounds coeffs[1] = aux
    @inbounds for k = 1:a.order
        coeffs[k+1] = expHomogCoef(k, v, coeffs)
    end
    Taylor1( coeffs, a.order )
end

# Homogeneous coefficients for exp
doc"""
    expHomogCoef(kcoef, ac, coeffs)

Compute the `k-th` expansion coefficient of $c = \exp(a)$ given by

\begin{equation*}
c_k = \frac{1}{k} \sum_{j=0}^{k-1} (k-j) a_{k-j} c_j,
\end{equation*}

with $a$ a `Taylor1` polynomial.

Inputs are the `kcoef`-th coefficient, the vector of the expansion coefficients
`ac` of $a$, and the already calculated expansion coefficients `coeffs` of `c`.
"""
function expHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, coeffs::Array{T,1})
    kcoef == 0 && return exp(ac[1])
    coefhomog = zero(T)
    @inbounds for i = 0:kcoef-1
        coefhomog += (kcoef-i) * ac[kcoef-i+1] * coeffs[i+1]
    end
    coefhomog = coefhomog/kcoef
    coefhomog
end

## Log ## 
doc"""
    log(a)

Return the Taylor expansion of $\log(a)$, of order `a.order`, for `a::Taylor1` polynomial.

For details on making the Taylor expansion, see [`TaylorSeries.logHomogCoef`](@ref).
"""
function log(a::Taylor1)
    ( firstnonzero(a)>0 ) && throw(
        ArgumentError("Impossible to expand `log` around 0."))
    @inbounds aux = log( a.coeffs[1] )
    T = typeof(aux)
    ac = convert(Array{T,1}, a.coeffs)
    coeffs = similar(ac)
    @inbounds coeffs[1] = aux
    @inbounds for k = 1:a.order
        coeffs[k+1] = logHomogCoef(k, ac, coeffs)
    end
    Taylor1( coeffs, a.order )
end

# Homogeneous coefficients for log
doc"""
    logHomogCoef(kcoef, ac, coeffs)

Compute the `k-th` expansion coefficient of $c = \log(a)$, given by

\begin{equation*}
c_k = \frac{1}{a_0} (a_k - \frac{1}{k} \sum_{j=0}^{k-1} j a_{k-j} c_j ),
\end{equation*}

with $a$ a `Taylor1` polynomial.

Inputs are the `kcoef`-th coefficient, the vector of the expansion coefficients
`ac` of `a`, and the already calculated expansion coefficients `coeffs` of `c`.
"""
function logHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, coeffs::Array{T,1})
    kcoef == 0 && return log( ac[1] )
    coefhomog = zero(T)
    @inbounds for i = 1:kcoef-1
        coefhomog += (kcoef-i) * ac[i+1] * coeffs[kcoef-i+1]
    end
    @inbounds coefhomog = (ac[kcoef+1] -coefhomog/kcoef) / ac[1]
    coefhomog
end


### TRIGONOMETRIC FUNCTIONS ###

## Sin ## 
doc"""
    sin(a)

Return the Taylor expansion of $\sin(a)$, of order `a.order`, for
`a::Taylor1` polynomial.

For details on making the Taylor expansion, see [`TaylorSeries.sincosHomogCoef`](@ref).
"""
sin(a::Taylor1) = sincos(a)[1]

## Cos ## 
doc"""
    cos(a)

Return the Taylor expansion of $\cos(a)$, of order `a.order`,
for `a::Taylor1` polynomial

For details on making the Taylor expansion, see [`TaylorSeries.sincosHomogCoef`](@ref).
"""
cos(a::Taylor1) = sincos(a)[2]

## Sin and Cos ## 
function sincos(a::Taylor1)
    @inbounds aux = sin( a.coeffs[1] )
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    sincoeffs = similar(v)
    coscoeffs = similar(v)
    @inbounds sincoeffs[1] = aux
    @inbounds coscoeffs[1] = cos( a.coeffs[1] )
    @inbounds for k = 1:a.order
        sincoeffs[k+1], coscoeffs[k+1] = sincosHomogCoef(k, v, sincoeffs, coscoeffs)
    end
    return Taylor1(sincoeffs, a.order), Taylor1(coscoeffs, a.order)
end

# Homogeneous coefficients for sincos
doc"""
    sincosHomogCoef(kcoef, ac, scoeffs, ccoeffs)

Compute the `k-th` expansion coefficient of $s = \sin(a)$ and $c=\cos(a)$
simultaneously given by

\begin{eqnarray*}
s_k &=& \frac{1}{k} \sum_{j=0}^{k-1} (k-j) a_{k-j} c_j \\\ \\

c_k &=& -\frac{1}{k}\sum_{j=0}^{k-1} (k-j) a_{k-j} s_j
\end{eqnarray*}

with $a$ a `Taylor1` polynomial.

Inputs are the `kcoef`-th coefficient, the vector of the expansion coefficients
`ac` of `a`, and the already calculated expansion coefficients `scoeffs`
and `ccoeffs` of `sin(a)` and `cos(a)`, respectvely.
"""
function sincosHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1},
    scoeffs::Array{T,1}, ccoeffs::Array{T,1})

    kcoef == 0 && return sin( ac[1] ), cos( ac[1] )
    sincoefhom = zero(T)
    coscoefhom = zero(T)

    @inbounds for i = 1:kcoef
        x = i * ac[i+1]
        sincoefhom += x * ccoeffs[kcoef-i+1]
        coscoefhom -= x * scoeffs[kcoef-i+1]
    end

    sincoefhom = sincoefhom/kcoef
    coscoefhom = coscoefhom/kcoef
    return sincoefhom, coscoefhom
end

## tan ##
doc"""
    tan(a)

Return the Taylor expansion of $\tan(a)$, of order `a.order`, for
`a::Taylor1` polynomial.

For details on making the Taylor expansion, see [`TaylorSeries.tanHomogCoef`](@ref).
"""
function tan(a::Taylor1)
    aux = tan( a.coeffs[1] )
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    coeffs = similar(v)
    coeffst2 = similar(v)
    @inbounds coeffs[1] = aux
    @inbounds coeffst2[1] = aux^2
    @inbounds for k = 1:a.order
        coeffs[k+1] = tanHomogCoef(k, v, coeffst2)
        coeffst2[k+1] = squareHomogCoef(k, coeffs)
    end
    Taylor1( coeffs, a.order )
end

# Homogeneous coefficients for tan
doc"""
    tanHomogCoef(kcoef, ac, coeffst2)

Compute the `k-th` expansion coefficient of $c = \tan(a)$ given by

\begin{equation*}
c_k = a_k + \frac{1}{k} \sum_{j=0}^{k-1} (k-j) a_{k-j} p_j,
\end{equation*}

with $a$ a `Taylor1` polynomial and $p = c^2$.

Inputs are the `kcoef`-th coefficient, the vector of the expansion coefficients
`ac` of `a`, and the already calculated expansion coefficients `coeffst2`
of `c^2`.
"""
function tanHomogCoef{T<:Number}(kcoef::Int,ac::Array{T,1},coeffst2::Array{T,1})
    kcoef == 0 && return tan( ac[1] )
    coefhomog = zero(T)
    @inbounds for i = 0:kcoef-1
        coefhomog += (kcoef-i)*ac[kcoef-i+1]*coeffst2[i+1]
    end
    @inbounds coefhomog = ac[kcoef+1] + coefhomog/kcoef
    coefhomog
end

### INVERSE TRIGONOMETRIC FUNCTIONS ### 

## Arcsin ##
doc"""
    asin(a)

Return the Taylor expansion of $\arcsin(a)$, of order `a.order`, for
`a::Taylor1` polynomial.

For details on making the Taylor expansion, see [`TaylorSeries.asinHomogCoef`](@ref).
"""
function asin(a::Taylor1)
    @inbounds aux = asin( a.coeffs[1] )
    T = typeof(aux)
    ac = convert(Array{T,1}, a.coeffs)
    ac[1]^2 == one(T) && throw(ArgumentError("""
        Recursion formula diverges due to vanishing `sqrt`."""))
    rc = sqrt(one(T) - a^2).coeffs
    coeffs = zeros(ac)
    @inbounds coeffs[1] = aux
    @inbounds for k in 1:a.order
        coeffs[k+1] = asinHomogCoef(k, ac, rc, coeffs)
    end
    Taylor1( coeffs, a.order )
end

# Homogeneous coefficients for arcsin
doc"""
    asinHomogCoef(kcoef, ac, rc, coeffs)

Compute the `k-th` expansion coefficient of $s = \arcsin(a)$ given by

\begin{equation*}
s_k = \frac{1}{ \sqrt{r_0} } \big( a_k - \frac{1}{k} \sum_{j=1}^{k-1} j r_{k-j} s_j \big),
\end{equation*}

with $a$ a `Taylor1` polynomial and $r = \sqrt{1 - a^2}$.

Inputs are the `kcoef`-th coefficient, the vector of the expansion coefficients
`ac` of `a`, the already calculated expansion coefficients `rc`
of $r$ (see above), and the already calculated expansion coefficients
`coeffs` of `asin(a)`.
"""
function asinHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, rc::Array{T,1},
    coeffs::Array{T,1})
    kcoef == 0 && return asin( ac[1] )
    coefhomog = zero(T)
    @inbounds for i in 1:kcoef-1
        coefhomog += (kcoef-i) * rc[i+1] * coeffs[kcoef-i+1]
    end
    @inbounds coefhomog = (ac[kcoef+1] - coefhomog/kcoef) / rc[1]
    coefhomog
end

## Arccos ## 
doc"""
    acos(a)

Return the Taylor expansion of $\arccos(a)$, of order `a.order`, for
`a::Taylor1` polynomial.

For details on making the Taylor expansion, see [`TaylorSeries.acosHomogCoef`](@ref).
"""
function acos(a::Taylor1)
    @inbounds aux = asin( a.coeffs[1] )
    T = typeof(aux)
    ac = convert(Array{T,1}, a.coeffs)
    ac[1]^2 == one(T) && throw(ArgumentError("""
        Recursion formula diverges due to vanishing `sqrt`."""))
    rc = sqrt(one(T) - a^2).coeffs
    coeffs = zeros(ac)
    @inbounds coeffs[1] = aux
    @inbounds for k in 1:a.order
        coeffs[k+1] = asinHomogCoef(k , ac, rc, coeffs)
    end
    @inbounds coeffs[1] = -acos( a.coeffs[1] )
    Taylor1( -coeffs, a.order )
end

# Homogeneous coefficients for arccos
doc"""
    acosHomogCoef(kcoef, ac, rc, coeffs)

Compute the `k-th` expansion coefficient of $c = \arccos(a)$ given by

\begin{equation*}
c_k = - \frac{1}{ r_0 } \big( a_k - \frac{1}{k} \sum_{j=1}^{k-1} j r_{k-j} c_j \big),
\end{equation*}

with $a$ a `Taylor1` polynomial and $r = \sqrt{1 - a^2}$.

Inputs are the `kcoef`-th coefficient, the vector of the expansion coefficients
`ac` of `a`, the already calculated expansion coefficients `rc`
of $r$ (see above), and the already calculated expansion coefficients
`coeffs` of `acos(a)`.
"""
function acosHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, rc::Array{T,1}, coeffs::Array{T,1})
    kcoef == 0 && return asin( ac[1] )
    coefhomog = zero(T)
    @inbounds for i in 1:kcoef-1
        coefhomog += (kcoef-i) * rc[i+1] * coeffs[kcoef-i+1]
    end
    @inbounds coefhomog = -(ac[kcoef+1] - coefhomog/kcoef) / rc[1]
    coefhomog
end


## Arctan
doc"""
    atan(a)

Return the Taylor expansion of $\arctan(a)$, of order `a.order`, for
`a::Taylor1` polynomial.

For details on making the Taylor expansion, see [`TaylorSeries.atanHomogCoef`](@ref).
"""
function atan(a::Taylor1)
    @inbounds aux = atan( a.coeffs[1] )
    T = typeof(aux)
    rc = (one(T) + a^2).coeffs
    rc[1] == zero(T) && throw(ArgumentError("""
        Recursion formula has a pole."""))
    ac = convert(Array{T,1}, a.coeffs)
    coeffs = zeros(ac)
    @inbounds coeffs[1] = aux
    @inbounds for k in 1:a.order
        coeffs[k+1] = atanHomogCoef(k , ac, rc, coeffs)
    end
    Taylor1( coeffs, a.order )
end

# Homogeneous coefficients for arctan
doc"""
    atanHomogCoef(kcoef, ac, rc, coeffs)

Compute the `k-th` expansion coefficient of $c = \arctan(a)$ given by

\begin{equation*}
t_k = \frac{1}{r_0}(a_k - \frac{1}{k} \sum_{j=1}^{k-1} j r_{k-j} t_j) ,
\end{equation*}

with $a$ a `Taylor1` polynomial and $r = 1 + a^2$.

Inputs are the `kcoef`-th coefficient, the vector of the expansion coefficients
`ac` of `a`, the already calculated expansion coefficients `rc`
of $r$ (see above), and the already calculated expansion coefficients
`coeffs` of `asin(a)`.
"""
function atanHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, rc::Array{T,1}, coeffs::Array{T,1})
    kcoef == 0 && return atan( ac[1] )
    coefhomog = zero(T)
    @inbounds for i in 1:kcoef-1
        coefhomog += (kcoef-i) * rc[i+1] * coeffs[kcoef-i+1]
    end
    @inbounds coefhomog = (ac[kcoef+1] - coefhomog/kcoef) / rc[1]
    coefhomog
end

## Differentiating ##
"""
    derivative(a)

Return the `Taylor1` polynomial of the differential of `a::Taylor1`.
The last coefficient is set to zero.
"""
function derivative(a::Taylor1)
    coeffs = zero(a.coeffs)
    @inbounds coeffs[1] = a.coeffs[2]
    @inbounds for i = 1:a.order
        coeffs[i] = i*a.coeffs[i+1]
    end
    return Taylor1(coeffs, a.order)
end

"""
    derivative(n, a)

Return the value of the `n`-th derivative of the polynomial `a`.
"""
function derivative{T<:Number}(n::Int, a::Taylor1{T})
    @assert a.order >= n >= 0
    factorial( widen(n) ) * a.coeffs[n+1] :: T
end

## Integrating ##
"""
    integrate(a, [x])

Return the integral of `a::Taylor1`. The constant of integration
(0-th order coefficient) is set to `x`, which is zero if ommitted.
"""
function integrate{T<:Number, S<:Number}(a::Taylor1{T}, x::S)
    R = promote_type(T, typeof(a.coeffs[1] / 1), S)
    coeffs = zeros(R, a.order+1)
    @inbounds for i = 1:a.order
        coeffs[i+1] = a.coeffs[i] / i
    end
    @inbounds coeffs[1] = convert(R, x)
    return Taylor1(coeffs, a.order)
end
integrate{T<:Number}(a::Taylor1{T}) = integrate(a, zero(T))

## Evaluating ##
"""
    evaluate(a, [dx])

Evaluate a `Taylor1` polynomial using Horner's rule (hand coded). If `dx` is
ommitted, its value is considered as zero.
"""
function evaluate{T<:Number}(a::Taylor1{T}, dx::T)
    @inbounds suma = a.coeffs[end]
    @inbounds for k = a.order:-1:1
        suma = suma*dx + a.coeffs[k]
    end
    suma
end
function evaluate{T<:Number,S<:Number}(a::Taylor1{T}, dx::S)
    R = promote_type(T,S)
    @inbounds suma = convert(R, a.coeffs[end])
    @inbounds for k = a.order:-1:1
        suma = suma*dx + a.coeffs[k]
    end
    suma
end
evaluate{T<:Number}(a::Taylor1{T}) = a.coeffs[1]

doc"""
    evaluate(x, δt)

Evaluates each element of `x::Array{Taylor1{T},1}`, representing
the dependent variables of an ODE, at *time* δt.
"""
function evaluate{T<:Number, S<:Number}(x::Array{Taylor1{T},1}, δt::S)
    R = promote_type(T,S)
    return evaluate(convert(Array{Taylor1{R},1},x), convert(R,δt))
end
function evaluate{T<:Number}(x::Array{Taylor1{T},1}, δt::T)
    xnew = Array{T}( length(x) )
    evaluate!(x, δt, xnew)
    return xnew
end
evaluate{T<:Number}(a::Array{Taylor1{T},1}) = evaluate(a, zero(T))

doc"""
    evaluate!(x, δt, x0)

Evaluates each element of `x::Array{Taylor1{T},1}`,
representing the Taylor expansion for the dependent variables
of an ODE at *time* δt. It updates the vector `x0` with the
computed values.
"""
function evaluate!{T<:Number}(x::Array{Taylor1{T},1}, δt::T, x0::Array{T,1})
    @assert length(x) == length(x0)
    @inbounds for i in eachindex(x)
        x0[i] = evaluate( x[i], δt )
    end
    nothing
end


"""
    evaluate(a, x)

Substitute `x::Taylor1` as independent variable in a `a::Taylor1` polynomial.
"""
function evaluate{T<:Number,S<:Number}(a::Taylor1{T}, x::Taylor1{S})
    a, x = fixshape(a, x)
    @inbounds suma = a.coeffs[end]
    @inbounds for k = a.order:-1:1
        suma = suma*x + a.coeffs[k]
    end
    suma
end


"""
    A_mul_B!(Y, A, B)

Multiply A*B and save the result in Y.
"""
function A_mul_B!{T<:Number}(y::Vector{Taylor1{T}}, a::Union{Matrix{T},SparseMatrixCSC{T}},
    b::Vector{Taylor1{T}})

    n, k = size(a)
    @assert (length(y)== n && length(b)== k)

    # determine the maximal order of b
    order = maximum([b1.order for b1 in b])

    # Use matrices of coefficients (of proper size) and A_mul_B!
    B = zeros(T, k, order+1)
    for i = 1:k
        @inbounds ord = b[i].order
        @inbounds for j = 1:ord+1
            B[i,j] = b[i].coeffs[j]
        end
    end
    Y = Array(T, n, order+1)
    A_mul_B!(Y, a, B)
    @inbounds for i = 1:n
        y[i] = Taylor1( collect(Y[i,:]), order)
    end

    return y
end
