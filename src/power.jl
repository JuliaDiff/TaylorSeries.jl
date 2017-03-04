# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#


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
        while (t -= 1) ≥ 0
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
    l0nz = findfirst(a)
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
    coeffs = Array{eltype(a)}(a.order+1)
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

For details on making the Taylor expansion, see
[`TaylorSeries.sqrtHomogCoef`](@ref).
"""
function sqrt(a::Taylor1)
    # First non-zero coefficient
    l0nz = findfirst(a)
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
function sqrtHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, coeffs::Array{T,1},
        knull::Int)
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



## Int power ##
function ^(a::HomogeneousPolynomial, n::Integer)
    n == 0 && return one(a)
    n == 1 && return a
    n == 2 && return square(a)
    n < 0 && throw(DomainError())
    return power_by_squaring(a, n)
end

function ^{T<:Number}(a::TaylorN{T}, n::Integer)
    n == 0 && return one(a)
    n == 1 && return a
    n == 2 && return square(a)
    n < 0 && return inv( a^(-n) )
    return power_by_squaring(a, n)
end

function ^{T<:Integer}(a::TaylorN{T}, n::Integer)
    n == 0 && return one(a)
    n == 1 && return a
    n == 2 && return square(a)
    n < 0 && throw(DomainError())
    return power_by_squaring(a, n)
end

## power_by_squaring; slightly modified from base/intfuncs.jl
## Licensed under MIT "Expat"
for T in (:HomogeneousPolynomial, :TaylorN)
    @eval begin
        function power_by_squaring(x::($T), p::Integer)
            p == 1 && return x# copy(x)
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
                while (t -= 1) ≥ 0
                    x = square(x)
                end
                y *= x
            end

            return y
        end
    end
end

## Rational power ##
^(a::TaylorN, x::Rational) = a^(x.num/x.den)

^(a::TaylorN, b::TaylorN) = exp( b*log(a) )

## Real power ##
function ^{S<:Real}(a::TaylorN, x::S)
    x == zero(x) && return TaylorN( one(eltype(a)), 0 )
    x == one(x)/2 && return sqrt(a)
    isinteger(x) && return a^round(Int,x)
    @inbounds a0 = a.coeffs[1].coeffs[1]
    @assert a0 != zero(a0)
    aux = ( a0 )^x
    T = typeof(aux)
    coeffs = zeros(HomogeneousPolynomial{T}, a.order)
    @inbounds coeffs[1] = HomogeneousPolynomial( [aux], 0 )

    for ord in eachindex(coeffs)
        ord == a.order+1 && continue
        @inbounds for i = 0:ord-1
            tt = x*(ord-i)-i
            cpol = coeffs[i+1]
            apol = a.coeffs[ord-i+1]
            (iszero(cpol) || iszero(apol)) && continue
            coeffs[ord+1] += tt * cpol * apol
        end
        coeffs[ord+1] = coeffs[ord+1] / (ord*a0)
    end

    return TaylorN{T}(coeffs, a.order)
end
^{T<:Complex}(a::TaylorN, x::T) = exp( x*log(a) )

## Square ##
function square(a::HomogeneousPolynomial)
    T = eltype(a)

    order = 2*a.order
    order > get_order() && return HomogeneousPolynomial([zero(T)], get_order())

    res = HomogeneousPolynomial([zero(T)], order)
    mul!(res, a)
    return res
end


function square(a::TaylorN)
    T = eltype(a)
    coeffs = zeros(HomogeneousPolynomial{T}, a.order)
    @inbounds mul!(coeffs[1], a.coeffs[1])

    for ord in eachindex(coeffs)
        ord == a.order+1 && continue
        kodd = ord%2
        kord = div(ord-2+kodd, 2)
        for i = 0 : kord
            @inbounds mul!(coeffs[ord+1], a.coeffs[i+1], a.coeffs[ord-i+1])
        end
        @inbounds coeffs[ord+1] = 2 * coeffs[ord+1]
        kodd == 1 && continue
        kodd = div(ord,2)
        @inbounds mul!(coeffs[ord+1], a.coeffs[kodd+1] )
    end

    return TaylorN{T}(coeffs, a.order)
end


## sqrt ##
function sqrt(a::TaylorN)
    @inbounds p0 = sqrt( a.coeffs[1].coeffs[1] )
    if p0 == zero(p0)
        throw(ArgumentError(
        """The 0-th order TaylorN coefficient must be non-zero
        in order to expand `sqrt` around 0."""))
    end

    T = typeof(p0)
    coeffs = zeros(HomogeneousPolynomial{T}, a.order)
    @inbounds coeffs[1] = HomogeneousPolynomial( [p0], 0 )

    for ord in eachindex(coeffs)
        ord == a.order+1 && continue
        kodd = ord%2
        kord = div(ord-2+kodd, 2)
        @inbounds for i = 1:kord
            coeffs[ord+1] += coeffs[i+1] * coeffs[ord-i+1]
        end
        @inbounds coeffs[ord+1] = a.coeffs[ord+1] - 2 * coeffs[ord+1]
        if iseven(ord)
            @inbounds coeffs[ord+1]=coeffs[ord+1]-square( coeffs[div(ord,2)+1])
        end
        @inbounds coeffs[ord+1] = coeffs[ord+1] / (2 * p0)
    end

    return TaylorN{T}(coeffs, a.order)
end
