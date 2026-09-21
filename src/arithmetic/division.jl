# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

## Division ##
function /(a::Taylor1{Rational{T}}, b::S) where {T<:Integer, S<:NumberNotSeries}
    R = typeof( a.coeffs[1] // b)
    v = FixedSizeVectorDefault{R}(undef, order(a)+1)
    v .= a.coeffs .// b
    return Taylor1(v)
end

function /(a::Taylor1{T}, b::S) where {T<:Number, S<:NumberNotSeries}
    R = typeof( a.coeffs[1] / b)
    v = FixedSizeVectorDefault{R}(undef, order(a)+1)
    v .= a.coeffs ./ b
    return Taylor1(v)
end

for T in (:HomogeneousPolynomial, :TaylorN)
    @eval function /(a::$T{T}, b::S) where {T<:NumberNotSeries, S<:NumberNotSeries}
        @inbounds aux = a.coeffs[1] / b
        v = FixedSizeVectorDefault{typeof(aux)}(undef, length(a.coeffs))
        v .= a.coeffs ./ b
        return $T(a.space, v, order(a))
    end

    @eval function /(b::$T{Taylor1{S}}, a::Taylor1{T}) where
            {T<:NumberNotSeries, S<:NumberNotSeries}
        @inbounds aux = b.coeffs[1] / a
        R = typeof(aux)
        coeffs = FixedSizeVectorDefault{R}(undef, length(b.coeffs))
        coeffs .= b.coeffs ./ a
        return $T(b.space, coeffs, order(b))
    end

    @eval function /(b::$T{Taylor1{T}}, a::S) where {T<:NumberNotSeries, S<:NumberNotSeries}
        @inbounds aux = b.coeffs[1] / a
        R = typeof(aux)
        coeffs = FixedSizeVectorDefault{R}(undef, length(b.coeffs))
        coeffs .= b.coeffs ./ a
        return $T(b.space, coeffs, order(b))
    end

    @eval function /(b::Taylor1{$T{S}}, a::$T{T}) where
            {T<:NumberNotSeries, S<:NumberNotSeries}
        @inbounds aux = b[0] / a
        v = Taylor1(zero(aux), order(b))
        @inbounds for k in eachindex(b)
            v[k] = b[k] / a
        end
        return v
    end
end

/(a::Taylor1{T}, b::Taylor1{S}) where {T<:Number, S<:Number} = /(promote(a,b)...)

function /(a::Taylor1{T}, b::Taylor1{T}) where {T<:Number}
    iszero(a) && !iszero(b) && return zero(a)
    if order(a) != order(b)
        a, b = fixorder(a, b)
    end
    # order and coefficient of first factorized term
    ordfact, cdivfact = divfactorization(a, b)
    R = typeof(cdivfact)
    if R == T
        aa = a
        bb = b
    else
        aa = convert(Taylor1{R}, a)
        bb = convert(Taylor1{R}, b)
    end
    c = Taylor1(cdivfact, order(a)-ordfact)
    for ord in eachindex(c)
        div!(c, aa, bb, ord) # updates c[ord]
    end
    return c
end

function /(a::Taylor1{Taylor1{T}}, b::Taylor1{S}) where
        {T<:NumberNotSeries, S<:NumberNotSeries}
    iszero(a) && !iszero(b) && return zero(a)
    cdivfact = constant_term(a) / b
    R = typeof(cdivfact)
    aa = R == T ? a : convert(Taylor1{R}, a)
    bb = R == S ? b : convert(R, b)
    c = Taylor1(cdivfact, order(a))
    for ord in eachindex(c)
        div!(c, aa, bb, ord) # updates c[ord]
    end
    return c
end

function /(a::Taylor1{T}, b::Taylor1{Taylor1{S}}) where
        {T<:NumberNotSeries, S<:NumberNotSeries}
    iszero(a) && !iszero(b) && return zero(a)
    cdivfact = a / constant_term(b)
    R = typeof(cdivfact)
    aa = R == T ? a : convert(R, a)
    bb = R == S ? b : convert(Taylor1{R}, b)
    c = Taylor1(cdivfact, order(b))
    for ord in eachindex(c)
        div!(c, aa, bb, ord) # updates c[ord]
    end
    return c
end


/(a::TaylorN{T}, b::TaylorN{S}) where
    {T<:NumberNotSeriesN, S<:NumberNotSeriesN} = /(promote(a,b)...)

function /(a::TaylorN{T}, b::TaylorN{T}) where {T<:NumberNotSeriesN}
    _check_same_space(a, b)
    @assert !_isthinzero(constant_term(b))
    if order(a) != order(b)
        a, b = fixorder(a, b)
    end
    # first coefficient
    @inbounds cdivfact = a[0] / constant_term(b)
    c = TaylorN(space(a), cdivfact, order(a))
    for ord in eachindex(c)
        div!(c, a, b, ord) # updates c[ord]
    end
    return c
end

function /(a::S, b::TaylorN{T}) where {S<:NumberNotSeriesN, T<:NumberNotSeriesN}
    @assert !_isthinzero(constant_term(b))
    R = typeof(a / constant_term(b))
    bb = convert(TaylorN{R}, b)
    res = TaylorN(space(b), zero(R), order(b))
    iszero(a) && !iszero(b) && return res
    aa = convert(R, a)
    for ord in eachindex(res)
        div!(res, aa, bb, ord)
    end
    return res
end

function /(a::Taylor1{TaylorN{T}}, b::Taylor1{TaylorN{T}}) where {T<:NumberNotSeries}
    _check_same_space(a[0], b[0])
    iszero(a) && !iszero(b) && return zero(a)
    if (order(a) != order(b)) || any(order.(a.coeffs) .!= order.(b.coeffs))
        a, b = fixorder(a, b)
    end
    # order and coefficient of first factorized term
    ordfact, cdivfact = divfactorization(a, b)
    R = numtype(cdivfact)
    if R == T
        aa = a
        bb = b
    else
        aa = convert(Taylor1{TaylorN{R}}, a)
        bb = convert(Taylor1{TaylorN{R}}, b)
    end
    res = Taylor1(cdivfact, order(a)-ordfact)
    for ordT in eachindex(res)
        div!(res, aa, bb, ordT)
    end
    return res
end

function /(a::S, b::Taylor1{TaylorN{T}}) where {S<:NumberNotSeries, T<:NumberNotSeries}
    R = promote_type(TaylorN{S}, TaylorN{T})
    res = convert(Taylor1{R}, zero(b))
    iszero(a) && !iszero(b) && return res
    for ordT in eachindex(res)
        div!(res, a, b, ordT)
    end
    return res
end

function /(a::TaylorN{T}, b::Taylor1{TaylorN{T}}) where {T<:NumberNotSeries}
    res = zero(b)
    iszero(a) && !iszero(b) && return res
    aa = Taylor1(a, order(b))
    for ordT in eachindex(res)
        div!(res, aa, b, ordT)
    end
    return res
end

## divfactorization ##

# Get order of first factorized term; a1 and b1 assumed to be of the same order
function _orderfactorizedterm(a1::Taylor1{T}, b1::Taylor1{T}) where {T}
    a1nz = findfirst(a1)
    b1nz = findfirst(b1)
    a1nz = a1nz ≥ 0 ? a1nz : order(a1)
    b1nz = b1nz ≥ 0 ? b1nz : order(a1)
    return min(a1nz, b1nz)
end

@inline function divfactorization(a1::Taylor1{T}, b1::Taylor1{T}) where {T}
    # order of first factorized term; a1 and b1 assumed to be of the same order
    ordfact = _orderfactorizedterm(a1, b1)
    cdivfact = a1.coeffs[ordfact+1] / b1.coeffs[ordfact+1]
    # Is the polynomial factorizable?
    TS._isthinzero(b1[ordfact]) && throw( ArgumentError(
        """Division does not define a Taylor1 polynomial;
        order k=$(ordfact) => coeff[$(ordfact)]=$(cdivfact).""") )

    return ordfact, cdivfact
end

# Similar to divfactorization, but writes the first order coefficient into `aux`
@inline function divfactorization!(aux::Taylor1{T}, a1::Taylor1{Taylor1{T}},
        b1::Taylor1{Taylor1{T}}, ordfact::Int) where {T<:NumberNotSeriesN}
    # Is the polynomial factorizable?
    TS._isthinzero(b1.coeffs[ordfact+1]) && throw( ArgumentError(
        """Division does not define a Taylor1 polynomial;
        order k=$(ordfact) => leading coefficient of the
        denominator is zero.""") )
    for k in eachindex(aux)
        zero!(aux, k)
        div!(aux, a1.coeffs[ordfact+1], b1.coeffs[ordfact+1], k)
    end
    return nothing
end



## TODO: Implement factorization (divfactorization) for TaylorN polynomials

@doc doc"""
    div!(c, a, b, k::Int)

Compute the `k-th` expansion coefficient `c[k]` of `c = a / b`,
where all `c`, `a` and `b` are either `Taylor1` or `TaylorN`.

The coefficients are given by

```math
c_k =  \frac{1}{b_0} \big(a_k - \sum_{j=0}^{k-1} c_j b_{k-j}\big).
```

For `Taylor1` polynomials, a similar formula is implemented which
exploits `k_0`, the order of the first non-zero coefficient of `a`.
""" div!

# @inline
function div!(c::Taylor1{T}, a::Taylor1{T}, b::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    kk = k+1
    @inbounds c_coeffs[kk] = zero(c_coeffs[kk])
    iszero(a) && !iszero(b) && return nothing
    # order and coefficient of first factorized term
    ordfact = _orderfactorizedterm(a, b)
    if k == 0
        @inbounds c_coeffs[1] = a_coeffs[ordfact+1] / b_coeffs[ordfact+1]
        return nothing
    end
    b_order = order(b)
    imin = max(0, k+ordfact-b_order)
    @inbounds acc = c_coeffs[imin+1] * b_coeffs[k+ordfact-imin+1]
    @inbounds for i = imin+1:k-1
        acc += c_coeffs[i+1] * b_coeffs[k+ordfact-i+1]
    end
    if k+ordfact ≤ b_order
        @inbounds acc = a_coeffs[k+ordfact+1] - acc
    else
        acc = -acc
    end
    @inbounds c_coeffs[kk] = acc / b_coeffs[ordfact+1]
    return nothing
end

# @inline
function div!(v::Taylor1{T}, a::Taylor1{T}, b::NumberNotSeries,
        k::Int) where {T<:Number}
    @inbounds v.coeffs[k+1] = a.coeffs[k+1] / b
    return nothing
end

function div!(v::Taylor1{T}, a::Taylor1{S}, b::NumberNotSeries) where
        {T<:NumberNotSeries, S<:NumberNotSeries}
    v_coeffs = v.coeffs
    a_coeffs = a.coeffs
    @inbounds for i in eachindex(v_coeffs)
        v_coeffs[i] = a_coeffs[i] / b
    end
    return nothing
end

function div!(v::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}},
        b::NumberNotSeries) where {T<:NumberNotSeries}
    v_coeffs = v.coeffs
    a_coeffs = a.coeffs
    @inbounds for i in eachindex(v_coeffs)
        div!(v_coeffs[i], a_coeffs[i], b)
    end
    return nothing
end

# @inline
function div!(c::Taylor1{T}, a::NumberNotSeries, b::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    c_coeffs = c.coeffs
    b_coeffs = b.coeffs
    kk = k+1
    @inbounds c_coeffs[kk] = zero(c_coeffs[kk])
    iszero(a) && !iszero(b) && return nothing
    if k == 0
        @inbounds c_coeffs[1] = a / b_coeffs[1]
        return nothing
    end
    @inbounds acc = c_coeffs[1] * b_coeffs[kk]
    @inbounds for i = 1:k-1
        acc += c_coeffs[i+1] * b_coeffs[k-i+1]
    end
    @inbounds c_coeffs[kk] = -acc / b_coeffs[1]
    return nothing
end

function div!(c::Taylor1, a::NumberNotSeries, b::Taylor1)
    @inbounds for k in eachindex(c)
        div!(c, a, b, k)
    end
    return nothing
end

@inline function div!(c::Taylor1{Taylor1{T}}, a::NumberNotSeries,
        b::Taylor1{Taylor1{T}}, k::Int) where {T<:NumberNotSeriesN}
    zero!(c, k)
    iszero(a) && !iszero(b) && return nothing
    c_coeffs = c.coeffs
    b_coeffs = b.coeffs
    if k == 0
        @inbounds div!(c_coeffs[1], a, b_coeffs[1])
        return nothing
    end
    kk = k + 1
    @inbounds mul!(c_coeffs[kk], c_coeffs[1], b_coeffs[kk])
    @inbounds for i = 1:k-1
        # c[k] += c[i] * b[k-i]
        muladd!(c_coeffs[kk], c_coeffs[i+1], b_coeffs[kk-i])
    end
    # @inbounds c[k] = -c[k]/b[0]
    @inbounds div_scalar!(c_coeffs[kk], -1, b_coeffs[1])
    return nothing
end

#
# @inline
function div!(c::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}},
        b::Taylor1{Taylor1{T}}, k::Int) where {T<:NumberNotSeriesN}
    zero!(c, k)
    iszero(a) && !iszero(b) && return nothing
    # order and coefficient of first factorized term
    ordfact = _orderfactorizedterm(a, b)
    if k == 0
        divfactorization!(c.coeffs[k+1], a, b, ordfact)
        return nothing
    end
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    b_order = order(b)
    kk = k + 1
    imin = max(0, k+ordfact-b_order)
    # c[k] = c[imin] * b[k+ordfact-imin]
    mul!(c_coeffs[kk], c_coeffs[imin+1], b_coeffs[kk+ordfact-imin])
    for i = imin+1:k-1
        # c[k] += c[i] * b[k+ordfact-i]
        for ord in eachindex(
                minlength(c_coeffs[kk], c_coeffs[i+1], b_coeffs[kk+ordfact-i]))
            muladd!(c_coeffs[kk], c_coeffs[i+1], b_coeffs[kk+ordfact-i], ord)
        end
    end
    aux = zero(c_coeffs[kk])
    if k+ordfact ≤ b_order
        # @inbounds aux <- a[k+ordfact] - c[k]
        for ord in eachindex(minlength(aux, a_coeffs[kk+ordfact]))
            subst!(aux, a_coeffs[kk+ordfact], c_coeffs[kk], ord)
        end
    else
        # @inbounds aux <- - c[k]
        for ord in eachindex(minlength(aux, a_coeffs[kk+ordfact]))
            subst!(aux, c_coeffs[kk], ord)
        end
    end
    # c[k] <- aux / b[ordfact]
    for ord in eachindex(c_coeffs[kk])
        div!(c_coeffs[kk], aux, b_coeffs[ordfact+1], ord)
    end
    return nothing
end

function div!(c::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}},
        b::Taylor1{T}, k::Int) where {T<:NumberNotSeriesN}
    zero!(c, k)
    iszero(a) && !iszero(b) && return nothing
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    kk = k + 1
    for j in eachindex(c_coeffs[kk])
        zero!(c_coeffs[kk], j)
        div!(c_coeffs[kk], a_coeffs[kk], b, j)
    end
    return nothing
end

function div!(c::Taylor1{Taylor1{T}}, a::Taylor1{T},
        b::Taylor1{Taylor1{T}}, k::Int) where {T<:NumberNotSeriesN}
    zero!(c.coeffs[k+1])
    iszero(a) && !iszero(b) && return nothing
    c_coeffs = c.coeffs
    b_coeffs = b.coeffs
    kk = k + 1
    # order and coefficient of first factorized term
    ordfact = _orderfactorizedterm(a, b_coeffs[kk])
    if k == 0
        c_coeffs[1] = a / b_coeffs[1]
        # divfactorization!(c_coeffs[kk], a, b_coeffs[kk], ordfact)
        return nothing
    end
    acc = zero(a)
    @inbounds for i = 0:k-1
        muladd!(acc, c_coeffs[i+1], b_coeffs[kk-i])
    end
    for j in eachindex(c_coeffs[1])
        subst!(acc, acc, j)
        div!(c_coeffs[kk], acc, b_coeffs[1], j)
    end
    return nothing
end

# @inline function div!(v::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}},
#         b::NumberNotSeries, k::Int) where {T<:NumberNotSeriesN}
#     # @inbounds v[k] = a[k] / b
#     for ord in eachindex(v)
#         div!(v, a, b, ord)
#     end
#     return nothing
# end

@inline function div!(c::Taylor1{TaylorN{T}}, a::NumberNotSeries,
        b::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
    zero!(c, k)
    iszero(a) && !iszero(b) && return nothing
    # order and coefficient of first factorized term
    # In this case, since a[k]=0 for k>0, we can simplify to:
    # ordfact, cdivfact = 0, a/b[0]
    if k == 0
        @inbounds div!(c.coeffs[1], a, b.coeffs[1])
        return nothing
    end
    c_coeffs = c.coeffs
    b_coeffs = b.coeffs
    kk = k + 1
    @inbounds mul!(c_coeffs[kk], c_coeffs[1], b_coeffs[kk])
    @inbounds for i = 1:k-1
        # c[k] += c[i] * b[k-i]
        mul!(c_coeffs[kk], c_coeffs[i+1], b_coeffs[kk-i])
    end
    # @inbounds c[k] = -c[k]/b[0]
    @inbounds div_scalar!(c_coeffs[kk], -1, b_coeffs[1])
    return nothing
end

# TODO: avoid allocations when T isa Taylor1
@inline function div!(v::HomogeneousPolynomial{T},
        a::HomogeneousPolynomial{T}, b::NumberNotSeriesN) where {T <: Number}
    _check_same_space(v, a)
    @inbounds for k in eachindex(v)
        v.coeffs[k] = a.coeffs[k] / b
    end
    return nothing
end

# NOTE: Due to the use of `zero!`, this `div!` method does *not* accumulate the result of a / b in c[k] (k > 0)
@inline function div!(c::TaylorN, a::TaylorN, b::TaylorN, k::Int)
    _check_same_space(c, a, b)
    if k==0
        @inbounds c.coeffs[1].coeffs[1] = constant_term(a) / constant_term(b)
        return nothing
    end
    zero!(c, k)
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    kk = k + 1
    @inbounds for i = 0:k-1
        mul!(c_coeffs[kk], c_coeffs[i+1], b_coeffs[kk-i])
    end
    @inbounds for i in eachindex(c_coeffs[kk])
        c_coeffs[kk].coeffs[i] =
            (a_coeffs[kk].coeffs[i] - c_coeffs[kk].coeffs[i]) / constant_term(b)
    end
    return nothing
end

# In-place division and assignment: c[k] = (c/a)[k]
# NOTE: Here `div!` *accumulates* the result of (c/a)[k] in c[k] (k > 0)
#
# Recursion algorithm:
#
# k = 0: c[0] <- c[0]/a[0]
# k = 1: c[1] <- c[1] - c[0]*a[1]
#        c[1] <- c[1]/a[0]
# k = 2: c[2] <- c[2] - c[0]*a[2] - c[1]*a[1]
#        c[2] <- c[2]/a[0]
# etc.
@inline function div!(c::TaylorN, a::TaylorN, k::Int)
    _check_same_space(c, a)
    if k==0
        @inbounds c.coeffs[1].coeffs[1] = constant_term(c) / constant_term(a)
        return nothing
    end
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    kk = k + 1
    @inbounds for i = 0:k-1
        mul_scalar!(c_coeffs[kk], -1, c_coeffs[i+1], a_coeffs[kk-i])
    end
    @inbounds for i in eachindex(c_coeffs[kk])
        c_coeffs[kk].coeffs[i] = c_coeffs[kk].coeffs[i] / constant_term(a)
    end
    return nothing
end

# In-place division and assignment: c[k] <- scalar * (c/a)[k]
# NOTE: Here `div!` *accumulates* the result of scalar * (c/a)[k] in c[k] (k > 0)
#
# Recursion algorithm:
#
# k = 0: c[0] <- scalar*c[0]/a[0]
# k = 1: c[1] <- scalar*c[1] - c[0]*a[1]
#        c[1] <- c[1]/a[0]
# k = 2: c[2] <- scalar*c[2] - c[0]*a[2] - c[1]*a[1]
#        c[2] <- c[2]/a[0]
# etc.
@inline function div_scalar!(c::TaylorN, scalar::NumberNotSeries,
        a::TaylorN, k::Int)
    _check_same_space(c, a)
    if k==0
        @inbounds c.coeffs[1].coeffs[1] =
            scalar * constant_term(c) / constant_term(a)
        return nothing
    end
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    kk = k + 1
    @inbounds mul!(c, scalar, c, k)
    @inbounds for i = 0:k-1
        mul_scalar!(c_coeffs[kk], -1, c_coeffs[i+1], a_coeffs[kk-i])
    end
    @inbounds for i in eachindex(c_coeffs[kk])
        c_coeffs[kk].coeffs[i] = c_coeffs[kk].coeffs[i] / constant_term(a)
    end
    return nothing
end

@inline function div_scalar!(c::Taylor1{T}, scalar::NumberNotSeries,
        a::Taylor1{T}, k::Int) where {T <: NumberNotSeries}
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    if k==0
        @inbounds c_coeffs[1] = scalar * c_coeffs[1] / a_coeffs[1]
        return nothing
    end
    kk = k+1
    @inbounds aux = scalar * c_coeffs[kk]
    @inbounds acc = zero(c_coeffs[kk])
    @inbounds for i = 0:k-1
        acc -= c_coeffs[i+1] * a_coeffs[k-i+1]
    end
    @inbounds c_coeffs[kk] = (acc + aux) / a_coeffs[1]
    return nothing
end

# NOTE: Here `div!` *accumulates* the result of a[k] / b[k] in c[k] (k > 0)
@inline function div!(c::TaylorN, a::NumberNotSeries, b::TaylorN, k::Int)
    _check_same_space(c, b)
    if k==0
        @inbounds c.coeffs[1].coeffs[1] = a / constant_term(b)
        return nothing
    end
    c_coeffs = c.coeffs
    b_coeffs = b.coeffs
    kk = k + 1
    @inbounds for i = 0:k-1
        mul!(c_coeffs[kk], c_coeffs[i+1], b_coeffs[kk-i])
    end
    @inbounds for i in eachindex(c_coeffs[kk])
        c_coeffs[kk].coeffs[i] = ( -c_coeffs[kk].coeffs[i] ) / constant_term(b)
    end
    return nothing
end

# c[k] <- a[k]/b, where b is a scalar
@inline function div!(c::TaylorN{T}, a::TaylorN{T}, b::NumberNotSeries,
        k::Int) where {T<:Number}
    _check_same_space(c, a)
    kk = k + 1
    @inbounds for i in eachindex(c.coeffs[kk])
        c.coeffs[kk].coeffs[i] = a.coeffs[kk].coeffs[i] / b
    end
    return nothing
end

# in-place division c <- c/a (assumes equal order among TaylorNs)
function div!(c::TaylorN, a::TaylorN)
    @inbounds for k in eachindex(c)
        div!(c, a, k)
    end
    return nothing
end

# in-place division c <- scalar*c/a (assumes equal order among TaylorNs)
function div_scalar!(c::TaylorN, scalar::NumberNotSeries, a::TaylorN)
    @inbounds for k in eachindex(c)
        div_scalar!(c, scalar, a, k)
    end
    return nothing
end

# in-place division c <- scalar*c/a (assumes equal order among TaylorNs)
function div_scalar!(c::Taylor1, scalar::NumberNotSeries, a::Taylor1)
    @inbounds for k in eachindex(c)
        div_scalar!(c, scalar, a, k)
    end
    return nothing
end

# c[k] <- (a/b)[k]
function div!(c::TaylorN, a::TaylorN, b::TaylorN)
    @inbounds for k in eachindex(c)
        div!(c, a, b, k)
    end
    return nothing
end

# c[k] <- (a/b)[k], where a is a scalar
function div!(c::TaylorN, a::NumberNotSeries, b::TaylorN)
    @inbounds for k in eachindex(c)
        div!(c, a, b, k)
    end
    return nothing
end

# c[k] <- a[k]/b, where b is a scalar
function div!(c::TaylorN{T}, a::TaylorN{T}, b::NumberNotSeries) where {T<:Number}
    @inbounds for k in eachindex(c)
        div!(c, a, b, k)
    end
    return nothing
end

# NOTE: Here `div!` *accumulates* the result of a / b in res[k] (k > 0)
@inline function div!(c::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        b::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeriesN}
    # order of first factorized term
    ordfact = _orderfactorizedterm(a, b)
    # Is the polynomial factorizable?
    _isthinzero(b.coeffs[ordfact+1]) && throw( ArgumentError(
        """Division does not define a Taylor1 polynomial;
        order k=$(ordfact) => leading coefficient of the
        denominator is zero.""") )
    zero!(c, k)
    if k == 0
        # @inbounds c[0] = a[ordfact]/b[ordfact]
        @inbounds div!(c.coeffs[1], a.coeffs[ordfact+1], b.coeffs[ordfact+1])
        return nothing
    end
    b_order = order(b)
    imin = max(0, k+ordfact-b_order)
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    kk = k + 1
    @inbounds mul!(c_coeffs[kk], c_coeffs[imin+1], b_coeffs[kk+ordfact-imin])
    @inbounds for i = imin+1:k-1
        mul!(c_coeffs[kk], c_coeffs[i+1], b_coeffs[kk+ordfact-i])
    end
        if k+ordfact ≤ b_order
        # @inbounds c[k] = (a[k+ordfact]-c[k]) / b[ordfact]
        @inbounds for l in eachindex(c_coeffs[kk])
            subst!(c_coeffs[kk], a_coeffs[kk+ordfact], c_coeffs[kk], l)
        end
        @inbounds div!(c_coeffs[kk], b_coeffs[ordfact+1])
    else
        # @inbounds c[k] = (-c[k]) / b[ordfact]
        @inbounds div_scalar!(c_coeffs[kk], -1, b_coeffs[ordfact+1])
    end
    return nothing
end

@inline function div!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        b::NumberNotSeries, k::Int) where {T<:NumberNotSeries}
    res_k = res.coeffs[k+1]
    a_k = a.coeffs[k+1]
    res_hps = res_k.coeffs
    a_hps = a_k.coeffs
    @inbounds for l in eachindex(res_hps)
        res_hp = res_hps[l].coeffs
        a_hp = a_hps[l].coeffs
        for m in eachindex(res_hp)
            res_hp[m] = a_hp[m]/b
        end
    end
    return nothing
end

