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



# power_by_squaring; slightly modified from base/intfuncs.jl
# Licensed under MIT "Expat"
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval function power_by_squaring(x::$T, p::Integer)
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
    zero!(c, a, k)

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
        @inbounds c[0] = ( constant_term(a) )^r
        return nothing
    end

    # Sanity
    zero!(c, a, k)

    # The recursion formula
    for i = 0:k-1
        aux = r*(k-i) - i
        mul!(c[k], aux*a[k-i], c[i])
    end
    @inbounds c[k] = c[k] / (k * constant_term(a))

    return nothing
end

@inline function pow!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}}, r::S,
        ordT::Int) where {T<:NumberNotSeries, S<:Real}

    if r == 0
        return one!(res, a, ordT)
    elseif r == 1
        return identity!(res, a, ordT)
    elseif r == 2
        return sqr!(res, a, ordT)
    elseif r == 0.5
        return sqrt!(res, a, ordT)
    end

    # Sanity
    zero!(res, a, ordT)

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
        res[ordT] = a[l0]^r # if r is integer, uses power_by_squaring internally
        return nothing
    end

    # The recursion formula
    tmp =  TaylorN( zero(constant_term(a[ordT])), a[0].order)
    tmp1 = TaylorN( zero(constant_term(a[ordT])), a[0].order)
    for i = 0:ordT-lnull-1
        ((i+lnull) > a.order || (l0+kprime-i > a.order)) && continue
        aux = r*(kprime-i) - i
        @inbounds for ordQ in eachindex(tmp)
            tmp1[ordQ] = aux * res[i+lnull][ordQ]
            mul!(tmp, tmp1, a[l0+kprime-i], ordQ)
        end
    end
    @inbounds for ordQ in eachindex(tmp)
        tmp1[ordQ] = tmp[ordQ]/kprime
        div!(res[ordT], tmp1, a[l0], ordQ)
    end

    return nothing
end


## Square ##
"""
    square(a::AbstractSeries) --> typeof(a)

Return `a^2`; see [`TaylorSeries.sqr!`](@ref).
""" square

for T in (:Taylor1, :TaylorN)
    @eval function square(a::$T)
        c = $T( constant_term(a)^2, a.order)
        for k in 1:a.order
            sqr!(c, a, k)
        end
        return c
    end
end

function square(a::HomogeneousPolynomial)
    order = 2*a.order
    order > get_order() && return HomogeneousPolynomial(zero(a[1]), get_order())
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
                @inbounds c[0] = constant_term(a)^2
                return nothing
            end

            # Sanity
            zero!(c, a, k)

            # Recursion formula
            kodd = k%2
            kend = (k - 2 + kodd) >> 1
            @inbounds for i = 0:kend
                if $T == Taylor1
                    c[k] += a[i] * a[k-i]
                else
                    mul!(c[k], a[i], a[k-i])
                end
            end
            @inbounds c[k] = 2 * c[k]
            kodd == 1 && return nothing

            if $T == Taylor1
                @inbounds c[k] += a[k >> 1]^2
            else
                accsqr!(c[k], a[k >> 1])
            end

            return nothing
        end
    end
end

@inline function sqr!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        ordT::Int) where {T<:NumberNotSeries}
    # Sanity
    zero!(res, a, ordT)

    if ordT == 0
        @inbounds for ordQ in eachindex(a[0])
            @inbounds sqr!(res[0], a[0], ordQ)
        end
        return nothing
    end

    # Recursion formula
    kodd = ordT%2
    kend = (ordT - 2 + kodd) >> 1
    tmp = TaylorN( zero(constant_term(a[0])), a[0].order )
    for i = 0:kend
        @inbounds for ordQ in eachindex(a[ordT])
            mul!(tmp, a[i], a[ordT-i], ordQ)
        end
    end
    @inbounds for ordQ in eachindex(a[ordT])
        res[ordT][ordQ] = 2 * tmp[ordQ]
    end
    kodd == 1 && return nothing

    @inbounds for ordQ in eachindex(a[0])
        zero!(tmp, a[0], ordQ)
        sqr!(tmp, a[ordT >> 1], ordQ)
        add!(res[ordT], res[ordT], tmp, ordQ)
    end

    return nothing
end


"""
    sqr!(c, a)

Return `c += a*a` with no allocation; all parameters are `HomogeneousPolynomial`.

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
    # Sanity
    zero!(c, a, k)

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
    imin ≤ imax && ( @inbounds c[k] = c[imin] * c[k+k0-imin] )
    @inbounds for i = imin+1:imax
        c[k] += c[i] * c[k+k0-i]
    end
    if k+k0 ≤ a.order
        @inbounds aux = a[k+k0] - 2*c[k]
    else
        @inbounds aux = - 2*c[k]
    end
    if kodd == 0
        @inbounds aux = aux - (c[kend+k0+1])^2
    end
    @inbounds c[k] = aux / (2*c[k0])

    return nothing
end

@inline function sqrt!(c::TaylorN{T}, a::TaylorN{T}, k::Int) where {T<:NumberNotSeriesN}
    # Sanity
    zero!(c, a, k)

    if k == 0
        @inbounds c[0] = sqrt( constant_term(a) )
        return nothing
    end

    # Recursion formula
    kodd = k%2
    kend = (k - 2 + kodd) >> 1
    @inbounds for i = 1:kend
        mul!(c[k], c[i], c[k-i])
    end
    @inbounds aux = a[k] - 2*c[k]
    if kodd == 0
        @inbounds aux = aux - (c[kend+1])^2
    end
    @inbounds c[k] = aux / (2*constant_term(c))

    return nothing
end

@inline function sqrt!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}}, ordT::Int,
        ordT0::Int=0) where {T<:NumberNotSeries}
    # Sanity
    zero!(res, a, ordT)
    ordT < ordT0 && return nothing

    if ordT == ordT0
        for ordQ in eachindex(a[0])
            sqrt!(res[ordT], a[2*ordT0], ordQ)
        end
        return nothing
    end

    # Recursion formula
    @inbounds tmp = TaylorN( zero( constant_term(a[ordT])), a[0].order)
    @inbounds tmp1 = TaylorN( zero( constant_term(a[ordT])), a[0].order)
    ordT_odd = (ordT - ordT0)%2
    ordT_end = (ordT - ordT0 - 2 + ordT_odd) >> 1
    imax = min(ordT0+ordT_end, a.order)
    imin = max(ordT0+1, ordT+ordT0-a.order)
    if imin ≤ imax
        @inbounds for ordQ in eachindex(a[0])
            mul!(tmp, res[imin], res[ordT+ordT0-imin], ordQ)
        end
    end
    for i = imin+1:imax
        @inbounds for ordQ in eachindex(a[0])
            mul!(tmp, res[i], res[ordT+ordT0-i], ordQ)
        end
    end

    @inbounds for ordQ in eachindex(a[0])
        tmp[ordQ] = -2 * tmp[ordQ]
        if ordT+ordT0 ≤ a.order
            add!(tmp, a[ordT+ordT0], tmp, ordQ)
        end
        if ordT_odd == 0
            sqr!(tmp1, res[ordT_end+ordT0+1], ordQ)
            subst!(tmp, tmp, tmp1, ordQ)
        end
        tmp1[ordQ] = 2*res[ordT0][ordQ]
        div!(res[ordT], tmp, tmp1, ordQ)
    end

    return nothing
end
