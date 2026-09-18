# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#


# Apparently necessary for v1.12
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval Base.literal_pow(::typeof(^), a::$T, ::Val{N}) where {N} = ^(a, N)
end

function ^(a::HomogeneousPolynomial, n::Integer)
    n == 0 && return one(a)
    n == 1 && return HomogeneousPolynomial(a.space, a.coeffs[:], order(a))
    n == 2 && return square(a)
    n < 0 && throw(DomainError())
    return power_by_squaring(a, n)
end


for T in (:Taylor1, :TaylorN)
    @eval function ^(a::$T, n::Integer)
        n == 0 && return one(a)
        n == 1 && return $T(a.coeffs, order(a))
        n == 2 && return square(a)
        return _pow(a, n)
    end

    @eval ^(a::$T, r::S) where {S<:Rational} = a^float(r)

    @eval ^(a::$T, b::$T) = exp( b*log(a) )

    @eval ^(a::$T, z::T) where {T<:Complex} = exp( z*log(a) )
end


## Real power ##
function ^(a::Taylor1{T}, r::S) where {T<:Number, S<:Real}
    a0 = constant_term(a)
    aux = a0^zero(r)
    iszero(r) && return Taylor1(aux, order(a))
    aa = aux*a
    r == 1 && return aa
    r == 2 && return square(aa)
    r == 0.5 && return sqrt(aa)
    return _pow(aa, r)
end

function ^(a::TaylorN{T}, r::S) where {T<:Number, S<:Real}
    a0 = constant_term(a)
    aux = a0^zero(r)
    iszero(r) && return TaylorN(a.space, aux, order(a))
    aa = aux*a
    r == 1 && return aa
    r == 2 && return square(aa)
    isinteger(r) && r >= 0 && return power_by_squaring(a, Integer(r))
    r == 0.5 && return sqrt(aa)
    if _isthinzero(a0)
        throw(DomainError(a,
        """The 0-th order TaylorN coefficient must be non-zero
        in order to expand `^` around 0."""))
    end
    return _pow(aa, r)
end


# _pow
_pow(a::Taylor1, n::Integer) = a^float(n)
_pow(a::Taylor1{T}, n::Integer) where {T<:NumberNotSeries} = a^float(n)

_pow(a::TaylorN, n::Integer) = power_by_squaring(a, n)

for T in (:Taylor1, :TaylorN)
    @eval _pow(a::$T{T}, n::Integer) where {T<:Integer} = power_by_squaring(a, n)

    @eval function _pow(a::$T{Rational{T}}, n::Integer) where {T<:Integer}
        n < 0 && return inv( a^(-n) )
        return power_by_squaring(a, n)
    end
end

function _pow(a::Taylor1{T}, r::S) where {T<:NumberNotSeries, S<:Real}
    aux = one(constant_term(a))^r
    iszero(r) && return Taylor1(aux, order(a))
    l0 = findfirst(a)
    lnull = trunc(Int, r*l0 )
    (lnull > order(a)) && return Taylor1( zero(aux), order(a))
    c_order = l0 == 0 ? order(a) : min(order(a), trunc(Int, r*order(a)))
    c = Taylor1(zero(aux), c_order)
    order_a = order(a)
    lastnz = isinteger(r) && r > 0 ? findlast(a) : order_a
    for k in eachindex(c)
        _pow_cached!(c, a, r, k, l0, lnull, lastnz, order_a)
    end
    return c
end

function _pow(a::Taylor1{T}, r::S) where {T<:Number, S<:Real}
    aux = one(constant_term(a))^r
    iszero(r) && return Taylor1(aux, order(a))
    l0 = findfirst(a)
    lnull = trunc(Int, r*l0 )
    (lnull > order(a)) && return Taylor1( zero(aux), order(a))
    c_order = l0 == 0 ? order(a) : min(order(a), trunc(Int, r*order(a)))
    c = Taylor1(zero(aux), c_order)
    aux0 = zero(c)
    for k in eachindex(c)
        pow!(c, a, aux0, r, k)
    end
    return c
end

function _pow(a::TaylorN{T}, r::S) where {T<:Number, S<:Real}
    isinteger(r) && r ≥ 0 && return power_by_squaring(a, Integer(r))
    aux = one(constant_term(a))^r
    c = TaylorN(a.space, zero(aux), order(a))
    aux0 = zero(c)
    for ord in eachindex(a)
        pow!(c, a, aux0, r, ord)
    end
    return c
end


# in-place form of power_by_squaring
# this method assumes `y`, `x` and `aux` are of same order
# TODO: add power_by_squaring! method for HomogeneousPolynomial and mixtures
for T in (:Taylor1, :TaylorN)
    @eval function power_by_squaring!(y::$T, x::$T, aux::$T, p::Integer)
        if p == 0
            for k in eachindex(y)
                one!(y, x, k)
            end
            return nothing
        end
        t = trailing_zeros(p) + 1
        p >>= t
        # aux = x
        for k in eachindex(aux)
            identity!(aux, x, k)
        end
        while (t -= 1) > 0
            # aux = square(aux)
            for k in reverse(eachindex(aux))
                sqr!(aux, k)
            end
        end
        # y = aux
        for k in eachindex(y)
            identity!(y, aux, k)
        end
        while p > 0
            t = trailing_zeros(p) + 1
            p >>= t
            while (t -= 1) ≥ 0
                # aux = square(aux)
                for k in reverse(eachindex(aux))
                    sqr!(aux, k)
                end
            end
            # y = y * aux
            mul!(y, aux)
        end
        return nothing
    end
end


# power_by_squaring; slightly modified from base/intfuncs.jl
# Licensed under MIT "Expat"
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval function Base.power_by_squaring(x::$T, p::Integer)
        @assert p ≥ 0
        (p == 0) && return one(x)
        (p == 1) && return $T(x.coeffs[:], order(x))
        (p == 2) && return square(x)
        (p == 3) && return x*square(x)
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

# power_by_squaring specializations for non-mixtures of Taylor1 and TaylorN;
# uses internally mutating method `power_by_squaring!`
for T in (:Taylor1, :TaylorN)
    @eval function Base.power_by_squaring(x::$T{T}, p::Integer) where {T<:NumberNotSeries}
        @assert p ≥ 0
        (p == 0) && return one(x)
        (p == 1) && return $T(x.coeffs[:], order(x))
        (p == 2) && return square(x)
        (p == 3) && return x*square(x)
        y = zero(x)
        aux = zero(x)
        power_by_squaring!(y, x, aux, p)
        return y
    end
end


# Homogeneous coefficients for real power
@doc doc"""
    pow!(c, a, aux, r::Real, k::Int)

Update the `k`-th expansion coefficient `c[k]` of `c = a^r`, for
both `c`, `a` and `aux` either `Taylor1` or `TaylorN`.

The coefficients are given by

```math
c_k = \frac{1}{k a_0} \sum_{j=0}^{k-1} \big(r(k-j) -j\big)a_{k-j} c_j.
```

For `Taylor1` polynomials, a similar formula is implemented which
exploits `k_0`, the order of the first non-zero coefficient of `a`.

""" pow!

@inline function _pow_cached!(c::Taylor1{T}, a::Taylor1{T},
        r::S, k::Int, l0::Int, lnull::Int, lastnz::Int,
        order_a::Int) where {T<:NumberNotSeries, S<:Real}
    zero!(c, k)
    l0 < 0 && return nothing
    !isinteger(r*l0) && throw(DomainError(a,
        """The 0-th order Taylor1 coefficient must be non-zero
        to raise the Taylor1 polynomial to a non-integer exponent."""))
    kprime = k-lnull
    (kprime < 0 || lnull > order_a) && return nothing
    isinteger(r) && r > 0 && (k > r*lastnz) && return nothing
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    kk = k+1
    @inbounds a_l0 = a_coeffs[l0+1]
    if k == lnull
        @inbounds c_coeffs[kk] = a_l0^float(r)
        return nothing
    end

    @inbounds acc = zero(c_coeffs[kk])
    if l0+kprime ≤ order_a
        @inbounds acc = r * kprime * c_coeffs[lnull+1] * a_coeffs[l0+kprime+1]
    end
    ilo = max(1, l0+kprime-order_a)
    ihi = min(k-lnull-1, order_a-lnull)
    @inbounds for i = ilo:ihi
        aaux = r*(kprime-i) - i
        acc += aaux * c_coeffs[i+lnull+1] * a_coeffs[l0+kprime-i+1]
    end
    @inbounds c_coeffs[kk] = acc / (kprime * a_l0)
    return nothing
end

function pow!(c::Taylor1{T}, a::Taylor1{T}, aux::Taylor1{T},
                r::S, k::Int) where {T<:NumberNotSeries, S<:Real}
    (r == 0) && return one!(c, a, k)
    (r == 1) && return identity!(c, a, k)
    (r == 2) && return sqr!(c, a, constant_term(aux), k)
    (r == 0.5) && return sqrt!(c, a, aux, k)
    l0 = findfirst(a)
    if l0 < 0
        zero!(c, k)
        return nothing
    end
    lnull = trunc(Int, r*l0)
    order_a = order(a)
    lastnz = isinteger(r) && r > 0 ? findlast(a) : order_a
    return _pow_cached!(c, a, r, k, l0, lnull, lastnz, order_a)
end

function pow!(c::TaylorN{T}, a::TaylorN{T}, aux::TaylorN{T},
                r::S, k::Int) where {T<:NumberNotSeriesN, S<:Real}
    isinteger(r) && r > 0 && return pow!(c, a, aux, Integer(r), k)
    (r == 0.5) && return sqrt!(c, a, aux, k)
    # 0-th order coeff
    if k == 0
        @inbounds c[0][1] = ( constant_term(a) )^r
        return nothing
    end
    # Sanity
    zero!(c, k)
    # The recursion formula
    @inbounds for i = 0:k-1
        aaux = r*(k-i) - i
        # c[k] += a[k-i]*c[i]*aaux
        mul_scalar!(c[k], aaux, a[k-i], c[i])
    end
    # c[k] <- c[k]/(k * constant_term(a))
    @inbounds div!(c[k], c[k], k * constant_term(a))
    return nothing
end

# Uses power_by_squaring!
function pow!(res::TaylorN{T}, a::TaylorN{T}, aux::TaylorN{T},
        r::S, k::Int) where {T<:NumberNotSeriesN, S<:Integer}
    (r == 0) && return one!(res, a, k)
    (r == 1) && return identity!(res, a, k)
    (r == 2) && return sqr!(res, a, constant_term(aux), k)
    power_by_squaring!(res, a, aux, r)
    return nothing
end

function pow!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        aux::Taylor1{TaylorN{T}}, r::S, ordT::Int) where
        {T<:NumberNotSeries, S<:Real}
    (r == 0) && return one!(res, a, ordT)
    (r == 1) && return identity!(res, a, ordT)
    (r == 2) && return sqr!(res, a, constant_term(aux), ordT)
    (r == 0.5) && return sqrt!(res, a, aux, ordT)
    # Sanity
    zero!(res, ordT)
    # First non-zero coefficient
    l0 = findfirst(a)
    l0 < 0 && return nothing
    # The first non-zero coefficient of the result; must be integer
    !isinteger(r*l0) && throw(DomainError(a,
        """The 0-th order Taylor1 coefficient must be non-zero
        to raise the Taylor1 polynomial to a non-integer exponent."""))
    lnull = trunc(Int, r*l0 )
    kprime = ordT-lnull
    (kprime < 0 || lnull > order(a)) && return nothing
    # Relevant for positive integer r, to avoid round-off errors
    isinteger(r) && r > 0 && (ordT > r*findlast(a)) && return nothing
    if ordT == lnull
        a0 = constant_term(a[l0])
        if isinteger(r) && r > 0
            # pow!(res[ordT], a[l0], aux[0], round(Integer, r), 1)
            power_by_squaring!(res[ordT], a[l0], aux[0], round(Integer, r))
            return nothing
        end
        _isthinzero(a0) && throw(DomainError(a[l0],
            """The 0-th order TaylorN coefficient must be non-zero
            in order to expand `^` around 0."""))
        # Recursion formula
        for ordQ in eachindex(a[l0])
            pow!(res[ordT], a[l0], aux[0], r, ordQ)
        end
        return nothing
    end
    # The recursion formula
    for i = 0:ordT-lnull-1
        ((i+lnull) > order(a) || (l0+kprime-i > order(a))) && continue
        aaux = r*(kprime-i) - i
        @inbounds mul_scalar!(res[ordT], aaux, res[i+lnull], a[l0+kprime-i])
    end
    # res[ordT] /= a[l0]*kprime
    @inbounds div_scalar!(res[ordT], 1/kprime, a[l0])
    return nothing
end

function pow!(c::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}},
        aux::Taylor1{Taylor1{T}}, r::S, k::Int) where
        {T<:NumberNotSeries, S<:Real}
    (r == 0) && return one!(c, a, k)
    (r == 1) && return identity!(c, a, k)
    (r == 2) && return sqr!(c, a, constant_term(aux), k)
    (r == 0.5) && return sqrt!(c, a, aux, k)
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    aux_coeffs = aux.coeffs
    kk = k+1
    @inbounds c_k = c_coeffs[kk]
    zero!(c_k)
    # First non-zero coefficient
    l0 = findfirst(a)
    l0 < 0 && return nothing
    # Index of first non-zero coefficient of the result; must be integer
    !isinteger(r*l0) && throw(DomainError(a,
        """The 0-th order Taylor1 coefficient must be non-zero
        to raise the Taylor1 polynomial to a non-integer exponent."""))
    lnull = trunc(Int, r*l0 )
    kprime = k-lnull
    order_a = order(a)
    (kprime < 0 || lnull > order_a) && return nothing
    # Relevant for positive integer r, to avoid round-off errors
    lastnz = isinteger(r) && r > 0 ? findlast(a) : order_a
    isinteger(r) && r > 0 && (k > r*lastnz) && return nothing
    @inbounds a_l0 = a_coeffs[l0+1]
    if k == lnull
        @inbounds aux0 = aux_coeffs[1]
        @inbounds for j in eachindex(a_l0)
            pow!(c_k, a_l0, aux0, float(r), j)
        end
        return nothing
    end
    # The recursion formula
    @inbounds aux_k = aux_coeffs[kk]
    ilo = max(0, l0+kprime-order_a)
    ihi = min(k-lnull-1, order_a-lnull)
    @inbounds for i = ilo:ihi
        rr = r*(kprime-i) - i
        mul_scalar!(aux_k, rr, c_coeffs[i+lnull+1], a_coeffs[l0+kprime-i+1])
        add!(c_k, c_k, aux_k)
    end
    # c[k] = c[k] / (kprime * a[l0])
    identity!(aux_k, c_k)
    @inbounds for j in eachindex(a_l0)
        div!(c_k, aux_k, a_l0, j)
    end
    @inbounds for j in eachindex(a_l0)
        div!(c_k, c_k, kprime, j)
    end
    return nothing
end

function pow!(c::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}},
        aux::Taylor1{Taylor1{T}}, r::S, k::Int) where
        {T<:NumberNotSeriesN, S<:Real}
    (r == 0) && return one!(c, a, k)
    (r == 1) && return identity!(c, a, k)
    (r == 2) && return sqr!(c, a, constant_term(aux), k)
    (r == 0.5) && return sqrt!(c, a, aux, k)
    # Sanity
    zero!(aux)
    zero!(c, k)
    # First non-zero coefficient
    l0 = findfirst(a)
    l0 < 0 && return nothing
    # Index of first non-zero coefficient of the result; must be integer
    !isinteger(r*l0) && throw(DomainError(a,
        """The 0-th order Taylor1 coefficient must be non-zero
        to raise the Taylor1 polynomial to a non-integer exponent."""))
    lnull = trunc(Int, r*l0 )
    kprime = k-lnull
    (kprime < 0 || lnull > order(a)) && return nothing
    # Relevant for positive integer r, to avoid round-off errors
    isinteger(r) && r > 0 && (k > r*findlast(a)) && return nothing
    # First non-zero coeff
    if k == lnull
        # @inbounds c[k] = (a[l0])^float(r)
        for j in eachindex(a[l0])
            pow!(c[k], a[l0], aux[0], float(r), j)
        end
        return nothing
    end
    # The recursion formula
    for i = 0:k-lnull-1
        ((i+lnull) > order(a) || (l0+kprime-i > order(a))) && continue
        rr = r*(kprime-i) - i
        # @inbounds c[k] += rr * c[i+lnull] * a[l0+kprime-i]
        @inbounds for j in eachindex(a[l0])
            mul_scalar!(aux[k], rr, c[i+lnull], a[l0+kprime-i], j)
            add!(c[k], c[k], aux[k], j)
        end
    end
    # @inbounds c[k] = c[k] / (kprime * a[l0])
    @inbounds for j in eachindex(c[k])
        identity!(aux[k], c[k], j)
    end
    @inbounds for j in eachindex(a[l0])
        div!(c[k], aux[k], a[l0], j)
    end
    @inbounds for j in eachindex(a[l0])
        div!(c[k], c[k], kprime, j)
    end
    return nothing
end
