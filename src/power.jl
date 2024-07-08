# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#


function ^(a::HomogeneousPolynomial, n::Integer)
    n == 0 && return one(a)
    n == 1 && return copy(a)
    n == 2 && return square(a)
    n < 0 && throw(DomainError())
    return power_by_squaring(a, n)
end

#= The following method computes `a^float(n)` (except for cases like
Taylor1{Interval{T}}^n, where `power_by_squaring` is used), to
use internally `pow!`.
=#
^(a::Taylor1, n::Integer) = a^float(n)

function ^(a::TaylorN{T}, n::Integer) where {T<:Number}
    n == 0 && return one(a)
    n == 1 && return copy(a)
    n == 2 && return square(a)
    n < 0 && return inv( a^(-n) )
    return power_by_squaring(a, n)
end


for T in (:Taylor1, :TaylorN)
    @eval function ^(a::$T{T}, n::Integer) where {T<:Integer}
        n == 0 && return one(a)
        n == 1 && return copy(a)
        n == 2 && return square(a)
        n < 0 && throw(DomainError())
        return power_by_squaring(a, n)
    end

    @eval function ^(a::$T{Rational{T}}, n::Integer) where {T<:Integer}
        n == 0 && return one(a)
        n == 1 && return copy(a)
        n == 2 && return square(a)
        n < 0 && return inv( a^(-n) )
        return power_by_squaring(a, n)
    end

    @eval ^(a::$T, x::S) where {S<:Rational} = a^(x.num/x.den)

    @eval ^(a::$T, b::$T) = exp( b*log(a) )

    @eval ^(a::$T, x::T) where {T<:Complex} = exp( x*log(a) )
end

^(a::Taylor1{TaylorN{T}}, n::Integer) where {T<:NumberNotSeries} = a^float(n)

^(a::Taylor1{TaylorN{T}}, r::Rational) where {T<:NumberNotSeries} = a^(r.num/r.den)


# in-place form of power_by_squaring
function power_by_squaring!(y::TaylorN, x::TaylorN, aux1::TaylorN, aux2::TaylorN, p::Integer)
    # this method assumes `y`, `aux1` and `aux2` are initialized to zero(x)
    t = trailing_zeros(p) + 1
    p >>= t
    # aux1 = x
    for k in eachindex(aux1)
        identity!(aux1, x, k)
    end
    while (t -= 1) > 0
        # TODO: square x in-place (in appropriate method)
        # y = square(aux1)
        for k in eachindex(y)
            sqr!(y, aux1, k)
        end
        # aux1 = y
        for k in eachindex(aux1)
            identity!(aux1, y, k)
        end
    end
    # y = aux1
    for k in eachindex(y)
        identity!(y, aux1, k)
    end
    while p > 0
        t = trailing_zeros(p) + 1
        p >>= t
        while (t -= 1) ≥ 0
            # TODO: square x in-place (in appropriate method)
            # aux2 = square(aux1)
            for k in eachindex(aux2)
                sqr!(aux2, aux1, k)
            end
            # aux1 = aux2
            for k in eachindex(aux1)
                identity!(aux1, aux2, k)
            end
        end
        # aux2 = y
        for k in eachindex(aux2)
            identity!(aux2, y, k)
        end
        # y = aux2 * aux1
        # if $T == Taylor1
        #     for k in eachindex(y)
        #         mul!(y, aux2, aux1, k)
        #     end
        # end
        zero!(y)
        for k in eachindex(y)
            mul!(y, aux2, aux1, k)
        end
    end
    return nothing
end

# power_by_squaring; slightly modified from base/intfuncs.jl
# Licensed under MIT "Expat"
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval power_by_squaring(x::$T, p::Integer) = power_by_squaring(x, Val(p))
    @eval power_by_squaring(x::$T, ::Val{0}) = one(x)
    @eval power_by_squaring(x::$T, ::Val{1}) = copy(x)
    @eval power_by_squaring(x::$T, ::Val{2}) = square(x)
    @eval function power_by_squaring(x::$T, ::Val{P}) where P
        p = P # copy static parameter `P` into local variable `p`
        if $T == TaylorN
            y = zero(x)
            aux1 = zero(x)
            aux2 = zero(x)
            power_by_squaring!(y, x, aux1, aux2, p)
        else
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
        end
        return y
    end
end


## Real power ##
function ^(a::Taylor1{T}, r::S) where {T<:Number, S<:Real}
    a0 = constant_term(a)
    aux = one(a0)^r

    iszero(r) && return Taylor1(aux, a.order)
    aa = aux*a
    r == 1 && return aa
    r == 2 && return square(aa)
    r == 0.5 && return sqrt(aa)

    l0 = findfirst(a)
    lnull = trunc(Int, r*l0 )
    (lnull > a.order) && return Taylor1( zero(aux), a.order)

    c_order = l0 == 0 ? a.order : min(a.order, trunc(Int,r*a.order))
    c = Taylor1(zero(aux), c_order)
    for k in eachindex(c)
        pow!(c, aa, r, k)
    end

    return c
end

## Real power ##
# TODO: get rid of allocations
function ^(a::TaylorN, r::S) where {S<:Real}
    a0 = constant_term(a)
    aux = one(a0^r)

    iszero(r) && return TaylorN(aux, a.order)
    aa = aux*a
    r == 1 && return aa
    r == 2 && return square(aa)
    r == 0.5 && return sqrt(aa)
    isinteger(r) && return aa^round(Int,r) # uses power_by_squaring

    iszero(a0) && throw(DomainError(a,
        """The 0-th order TaylorN coefficient must be non-zero
        in order to expand `^` around 0."""))

    c = TaylorN( zero(aux), a.order)
    for ord in eachindex(a)
        pow!(c, aa, r, ord)
    end

    return c
end

function ^(a::Taylor1{TaylorN{T}}, r::S) where {T<:NumberNotSeries, S<:Real}
    a0 = constant_term(a)
    aux = one(a0)^r

    iszero(r) && return Taylor1(aux, a.order)
    aa = aux*a
    r == 1 && return aa
    r == 2 && return square(aa)
    r == 0.5 && return sqrt(aa)
    # Is the following needed?
    # isinteger(r) && r > 0 && iszero(constant_term(a[0])) &&
    #     return power_by_squaring(aa, round(Int,r))

    l0 = findfirst(a)
    lnull = trunc(Int, r*l0 )
    (lnull > a.order) && return Taylor1( zero(aux), a.order)

    c_order = l0 == 0 ? a.order : min(a.order, trunc(Int,r*a.order))
    c = Taylor1(zero(aux), c_order)
    for k in eachindex(c)
        pow!(c, aa, r, k)
    end

    return c
end


# Homogeneous coefficients for real power
@doc doc"""
    pow!(c, a, r::Real, k::Int)

Update the `k`-th expansion coefficient `c[k]` of `c = a^r`, for
both `c` and `a` either `Taylor1` or `TaylorN`.

The coefficients are given by

```math
c_k = \frac{1}{k a_0} \sum_{j=0}^{k-1} \big(r(k-j) -j\big)a_{k-j} c_j.
```

For `Taylor1` polynomials, a similar formula is implemented which
exploits `k_0`, the order of the first non-zero coefficient of `a`.

""" pow!

@inline function pow!(c::Taylor1{T}, a::Taylor1{T}, r::S, k::Int) where
        {T<:Number, S<:Real}

    if r == 0
        return one!(c, a, k)
    elseif r == 1
        return identity!(c, a, k)
    elseif r == 2
        return sqr!(c, a, k)
    elseif r == 0.5
        return sqrt!(c, a, k)
    end

    # Sanity
    zero!(c, k)

    # First non-zero coefficient
    l0 = findfirst(a)
    l0 < 0 && return nothing

    # The first non-zero coefficient of the result; must be integer
    !isinteger(r*l0) && throw(DomainError(a,
        """The 0-th order Taylor1 coefficient must be non-zero
        to raise the Taylor1 polynomial to a non-integer exponent."""))
    lnull = trunc(Int, r*l0 )
    kprime = k-lnull
    (kprime < 0 || lnull > a.order) && return nothing

    # Relevant for positive integer r, to avoid round-off errors
    isinteger(r) && r > 0 && (k > r*findlast(a)) && return nothing

    if k == lnull
        @inbounds c[k] = (a[l0])^r
        return nothing
    end

    # The recursion formula
    if l0+kprime ≤ a.order
        @inbounds c[k] = r * kprime * c[lnull] * a[l0+kprime]
    end
    for i = 1:k-lnull-1
        ((i+lnull) > a.order || (l0+kprime-i > a.order)) && continue
        aux = r*(kprime-i) - i
        @inbounds c[k] += aux * c[i+lnull] * a[l0+kprime-i]
    end
    @inbounds c[k] = c[k] / (kprime * a[l0])
    return nothing
end

@inline function pow!(c::TaylorN{T}, a::TaylorN{T}, r::S, k::Int) where
        {T<:NumberNotSeriesN, S<:Real}

    if r == 0
        return one!(c, a, k)
    elseif r == 1
        return identity!(c, a, k)
    elseif r == 2
        return sqr!(c, a, k)
    elseif r == 0.5
        return sqrt!(c, a, k)
    end

    if k == 0
        @inbounds c[0][1] = ( constant_term(a) )^r
        return nothing
    end

    # Sanity
    zero!(c, k)

    # The recursion formula
    @inbounds for i = 0:k-1
        aux = r*(k-i) - i
        # c[k] += a[k-i]*c[i]*aux
        mul_scalar!(c[k], aux, a[k-i], c[i])
    end
    # c[k] <- c[k]/(k * constant_term(a))
    @inbounds div!(c[k], c[k], k * constant_term(a))

    return nothing
end

@inline function pow!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}}, r::S,
        ordT::Int) where {T<:NumberNotSeries, S<:Real}

    if r == 0
        return one!(res, ordT)
    elseif r == 1
        return identity!(res, a, ordT)
    elseif r == 2
        return sqr!(res, a, ordT)
    elseif r == 0.5
        return sqrt!(res, a, ordT)
    end

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
    (kprime < 0 || lnull > a.order) && return nothing

    # Relevant for positive integer r, to avoid round-off errors
    isinteger(r) && r > 0 && (ordT > r*findlast(a)) && return nothing

    if ordT == lnull
        if isinteger(r)
            # TODO: get rid of allocations here
            @time res[ordT] = a[l0]^round(Int,r) # uses power_by_squaring
            return nothing
        end

        a0 = constant_term(a[l0])
        iszero(a0) && throw(DomainError(a[l0],
            """The 0-th order TaylorN coefficient must be non-zero
            in order to expand `^` around 0."""))

        for ordQ in eachindex(a[l0])
            pow!(res[ordT], a[l0], r, ordQ)
        end

        return nothing
    end

    # The recursion formula
    for i = 0:ordT-lnull-1
        ((i+lnull) > a.order || (l0+kprime-i > a.order)) && continue
        aux = r*(kprime-i) - i
        # res[ordT] += aux*res[i+lnull]*a[l0+kprime-i]
        @inbounds mul_scalar!(res[ordT], aux, res[i+lnull], a[l0+kprime-i])
    end
    # res[ordT] /= a[l0]*kprime
    @inbounds div_scalar!(res[ordT], 1/kprime, a[l0])

    return nothing
end


## Square ##
"""
    square(a::AbstractSeries) --> typeof(a)

Return `a^2`; see [`TaylorSeries.sqr!`](@ref).
""" square

for T in (:Taylor1, :TaylorN)
    @eval function square(a::$T)
        # c = $T( zero(constant_term(a)), a.order)
        # c = deepcopy(a)
        # zero!(c)
        c = zero(a)
        for k in eachindex(a)
            sqr!(c, a, k)
        end
        return c
    end
end

function square(a::HomogeneousPolynomial)
    order = 2*get_order(a)
    # NOTE: the following returns order 0, but could be get_order(), or get_order(a)
    order > get_order() && return HomogeneousPolynomial(zero(a[1]), 0)
    res = HomogeneousPolynomial(zero(a[1]), order)
    accsqr!(res, a)
    return res
end

function square(a::Taylor1{TaylorN{T}}) where {T<:NumberNotSeries}
    res = Taylor1(zero(a[0]), a.order)
    for ordT in eachindex(a)
        sqr!(res, a, ordT)
    end
    return res
end

#auxiliary function to avoid allocations
@inline function sqr_orderzero!(c::Taylor1{T}, a::Taylor1{T}) where {T<:NumberNotSeries}
    @inbounds c[0] = a[0]^2
    return nothing
end
@inline function sqr_orderzero!(c::TaylorN{T}, a::TaylorN{T}) where {T<:NumberNotSeries}
    @inbounds c[0][1] = a[0][1]^2
    return nothing
end
@inline function sqr_orderzero!(c::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}}) where {T<:NumberNotSeries}
    @inbounds for ord in eachindex(c[0])
        sqr!(c[0], a[0], ord)
    end
    return nothing
end
@inline function sqr_orderzero!(c::TaylorN{Taylor1{T}}, a::TaylorN{Taylor1{T}}) where {T<:NumberNotSeries}
    @inbounds for ord in eachindex(c[0][1])
        sqr!(c[0][1], a[0][1], ord)
    end
    return nothing
end
@inline function sqr_orderzero!(c::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}}) where {T<:NumberNotSeriesN}
    @inbounds for ord in eachindex(c[0])
        sqr!(c[0], a[0], ord)
    end
    return nothing
end

# Homogeneous coefficients for square
@doc doc"""
    sqr!(c, a, k::Int) --> nothing

Update the `k-th` expansion coefficient `c[k]` of `c = a^2`, for
both `c` and `a` either `Taylor1` or `TaylorN`.

The coefficients are given by

```math
\begin{aligned}
c_k &= 2 \sum_{j=0}^{(k-1)/2} a_{k-j} a_j,
    \text{ if $k$ is odd,} \\
c_k &= 2 \sum_{j=0}^{(k-2)/2} a_{k-j} a_j + (a_{k/2})^2,
    \text{ if $k$ is even.}
\end{aligned}
```

""" sqr!

for T = (:Taylor1, :TaylorN)
    @eval begin
        @inline function sqr!(c::$T{T}, a::$T{T}, k::Int) where {T<:Number}
            if k == 0
                sqr_orderzero!(c, a)
                return nothing
            end

            # Sanity
            zero!(c, k)

            # Recursion formula
            kodd = k%2
            kend = (k - 2 + kodd) >> 1
            if $T == Taylor1
                @inbounds for i = 0:kend
                    c[k] += a[i] * a[k-i]
                end
                @inbounds c[k] = 2 * c[k]
            else
                @inbounds for i = 0:kend
                    mul!(c[k], a[i], a[k-i])
                end
                @inbounds mul!(c, 2, c, k)
            end
            kodd == 1 && return nothing

            if $T == Taylor1
                @inbounds c[k] += a[k >> 1]^2
            else
                accsqr!(c[k], a[k >> 1])
            end

            return nothing
        end

        # in-place squaring: given `c`, compute expansion of `c^2` and save back into `c`
        @inline function sqr!(c::$T{T}, k::Int) where {T<:NumberNotSeries}
            if k == 0
                sqr_orderzero!(c, c)
                return nothing
            end

            # Recursion formula
            kodd = k%2
            kend = (k - 2 + kodd) >> 1
            if $T == Taylor1
                ck = deepcopy(c[k])
                (kodd == 0) && ( @inbounds c[k] = (c[k >> 1]^2)/2 )
                c[k] = c[0] * ck
                @inbounds for i = 1:kend
                    c[k] += c[i] * c[k-i]
                end
                @inbounds c[k] = 2 * c[k]
            else
                # zz = deepcopy(c[k])
                # zero!(c[k])
                # @show c[0][1], typeof(c[0][1])
                # _2c01 = 2c[0][1]
                # @inbounds mul!(c, _2c01, c, k)
                # mul!(c[k], zz, c[0])
                (kend ≥ 0) && mul!(c, c[0][1], c, k)
                @inbounds for i = 1:kend
                    mul!(c[k], c[i], c[k-i])
                end
                @inbounds mul!(c, 2, c, k)
                if (kodd == 0)
                    accsqr!(c[k], c[k >> 1])
                    # mul!(c, 0.5, c, k)
                end
            end

            return nothing
        end
    end
end

@inline function sqr!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        ordT::Int) where {T<:NumberNotSeries}
    # Sanity
    zero!(res, ordT)
    if ordT == 0
        @inbounds for ordQ in eachindex(a[0])
            @inbounds sqr!(res[0], a[0], ordQ)
        end
        return nothing
    end

    # Recursion formula
    kodd = ordT%2
    kend = (ordT - 2 + kodd) >> 1
    (kodd == 0) && @inbounds for ordQ in eachindex(a[0])
        sqr!(res[ordT], a[ordT >> 1], ordQ)
        mul!(res[ordT], 0.5, res[ordT], ordQ)
    end

    for i = 0:kend
        @inbounds for ordQ in eachindex(a[ordT])
            # mul! accumulates the result in res[ordT]
            mul!(res[ordT], a[i], a[ordT-i], ordQ)
        end
    end
    @inbounds for ordQ in eachindex(a[ordT])
        mul!(res[ordT], 2, res[ordT], ordQ)
    end

    return nothing
end


"""
    accsqr!(c, a)

Returns `c += a*a` with no allocation; all parameters are `HomogeneousPolynomial`.

"""
@inline function accsqr!(c::HomogeneousPolynomial{T}, a::HomogeneousPolynomial{T}) where
        {T<:NumberNotSeriesN}
    iszero(a) && return nothing

    @inbounds num_coeffs_a = size_table[a.order+1]

    @inbounds posTb = pos_table[c.order+1]
    @inbounds idxTb = index_table[a.order+1]

    @inbounds for na = 1:num_coeffs_a
        ca = a[na]
        iszero(ca) && continue
        inda = idxTb[na]
        pos = posTb[2*inda]
        c[pos] += ca^2
        @inbounds for nb = na+1:num_coeffs_a
            cb = a[nb]
            iszero(cb) && continue
            indb = idxTb[nb]
            pos = posTb[inda+indb]
            c[pos] += 2 * ca * cb
        end
    end

    return nothing
end



## Square root ##
function sqrt(a::Taylor1{T}) where {T<:Number}
    # First non-zero coefficient
    l0nz = findfirst(a)
    aux = zero(sqrt( constant_term(a) ))
    if l0nz < 0
        return Taylor1(aux, a.order)
    elseif isodd(l0nz) # l0nz must be pair
        throw(DomainError(a,
            """First non-vanishing Taylor1 coefficient must correspond
            to an **even power** in order to expand `sqrt` around 0."""))
    end

    # The last l0nz coefficients are dropped.
    lnull = l0nz >> 1 # integer division by 2
    c_order = l0nz == 0 ? a.order : a.order >> 1
    c = Taylor1( aux, c_order )
    aa = convert(Taylor1{eltype(aux)}, a)
    for k in eachindex(c)
        sqrt!(c, aa, k, lnull)
    end

    return c
end

function sqrt(a::TaylorN)
    @inbounds p0 = sqrt( constant_term(a) )
    if iszero(p0)
        throw(DomainError(a,
            """The 0-th order TaylorN coefficient must be non-zero
            in order to expand `sqrt` around 0."""))
    end

    c = TaylorN( p0, a.order)
    aa = convert(TaylorN{eltype(p0)}, a)
    for k in eachindex(c)
        sqrt!(c, aa, k)
    end

    return c
end

function sqrt(a::Taylor1{TaylorN{T}}) where {T<:NumberNotSeries}
    # First non-zero coefficient
    l0nz = findfirst(a)
    aux = TaylorN( zero(sqrt(constant_term(a[0]))), a[0].order )
    if l0nz < 0
        return Taylor1( aux, a.order )
    elseif isodd(l0nz) # l0nz must be pair
        throw(DomainError(a,
            """First non-vanishing Taylor1 coefficient must correspond
            to an **even power** in order to expand `sqrt` around 0."""))
    end

    # The last l0nz coefficients are dropped.
    lnull = l0nz >> 1 # integer division by 2
    c_order = l0nz == 0 ? a.order : a.order >> 1
    c = Taylor1( aux, c_order )
    aa = convert(Taylor1{eltype(aux)}, a)
    for k in eachindex(c)
        sqrt!(c, aa, k, lnull)
    end

    return c
end


# Homogeneous coefficients for the square-root
@doc doc"""
    sqrt!(c, a, k::Int, k0::Int=0)

Compute the `k-th` expansion coefficient `c[k]` of `c = sqrt(a)`
for both`c` and `a` either `Taylor1` or `TaylorN`.

The coefficients are given by

```math
\begin{aligned}
c_k &= \frac{1}{2 c_0} \big( a_k - 2 \sum_{j=1}^{(k-1)/2} c_{k-j}c_j\big),
    \text{ if $k$ is odd,} \\
c_k &= \frac{1}{2 c_0} \big( a_k - 2 \sum_{j=1}^{(k-2)/2} c_{k-j}c_j
    - (c_{k/2})^2\big), \text{ if $k$ is even.}
\end{aligned}
```

For `Taylor1` polynomials, `k0` is the order of the first non-zero
coefficient, which must be even.

""" sqrt!

@inline function sqrt!(c::Taylor1{T}, a::Taylor1{T}, k::Int, k0::Int=0) where {T<:Number}

    k < k0 && return nothing

    if k == k0
        @inbounds c[k] = sqrt(a[2*k0])
        return nothing
    end

    # Recursion formula
    kodd = (k - k0)%2
    # kend = div(k - k0 - 2 + kodd, 2)
    kend = (k - k0 - 2 + kodd) >> 1
    imax = min(k0+kend, a.order)
    imin = max(k0+1, k+k0-a.order)
    if k+k0 ≤ a.order
        @inbounds c[k] = a[k+k0]
    end
    if kodd == 0
        @inbounds c[k] -= (c[kend+k0+1])^2
    end
    imin ≤ imax && ( @inbounds c[k] -= 2 * c[imin] * c[k+k0-imin] )
    @inbounds for i = imin+1:imax
        c[k] -= 2 * c[i] * c[k+k0-i]
    end
    @inbounds c[k] = c[k] / (2*c[k0])

    return nothing
end

@inline function sqrt!(c::TaylorN{T}, a::TaylorN{T}, k::Int) where {T<:NumberNotSeriesN}

    if k == 0
        @inbounds c[0][1] = sqrt( constant_term(a) )
        return nothing
    end

    # Recursion formula
    kodd = k%2
    kend = (k - 2 + kodd) >> 1
    # c[k] <- a[k]
    @inbounds for i in eachindex(c[k])
        c[k][i] = a[k][i]
    end
    if kodd == 0
        # @inbounds c[k] <- c[k] - (c[kend+1])^2
        @inbounds mul_scalar!(c[k], -1, c[kend+1], c[kend+1])
    end
    @inbounds for i = 1:kend
        # c[k] <- c[k] - 2*c[i]*c[k-i]
        mul_scalar!(c[k], -2, c[i], c[k-i])
    end
    # @inbounds c[k] <- c[k] / (2*c[0])
    div!(c[k], c[k], 2*constant_term(c))

    return nothing
end

@inline function sqrt!(c::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}}, k::Int,
        k0::Int=0) where {T<:NumberNotSeries}

    k < k0 && return nothing

    if k == k0
        @inbounds for l in eachindex(c[k])
            sqrt!(c[k], a[2*k0], l)
        end
        return nothing
    end

    # Recursion formula
    kodd = (k - k0)%2
    # kend = div(k - k0 - 2 + kodd, 2)
    kend = (k - k0 - 2 + kodd) >> 1
    imax = min(k0+kend, a.order)
    imin = max(k0+1, k+k0-a.order)
    if k+k0 ≤ a.order
        # @inbounds c[k] += a[k+k0]
        ### TODO: add in-place add! method for Taylor1, TaylorN and mixtures: c[k] += a[k] -> add!(c, a, k)
        ###       and/or add identity! method such that each coeff is copied individually,
        ###       otherwise memory-mixing issues happen
        @inbounds for l in eachindex(c[k])
            for m in eachindex(c[k][l])
                c[k][l][m] = a[k+k0][l][m]
            end
        end
    end
    if kodd == 0
        # c[k] <- c[k] - c[kend+1]^2
        # TODO: use accsqr! here?
        @inbounds mul_scalar!(c[k], -1, c[kend+k0+1], c[kend+k0+1])
    end
    @inbounds for i = imin:imax
        # c[k] <- c[k] - 2 * c[i] * c[k+k0-i]
        mul_scalar!(c[k], -2, c[i], c[k+k0-i])
    end
    # @inbounds c[k] <- c[k] / (2*c[k0])
    @inbounds div_scalar!(c[k], 0.5, c[k0])

    return nothing
end
