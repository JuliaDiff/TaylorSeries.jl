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

for T in (:Taylor1, :TaylorN)
    @eval function ^(a::$T{T}, n::Integer) where {T<:Number}
        n == 0 && return one(a)
        n == 1 && return copy(a)
        n == 2 && return square(a)
        n < 0 && return inv( a^(-n) )
        return power_by_squaring(a, n)
    end

    @eval function ^(a::$T{T}, n::Integer) where {T<:Integer}
        n == 0 && return one(a)
        n == 1 && return copy(a)
        n == 2 && return square(a)
        n < 0 && throw(DomainError())
        return power_by_squaring(a, n)
    end

    @eval ^(a::$T, x::Rational) = a^(x.num/x.den)

    @eval ^(a::$T, b::$T) = exp( b*log(a) )

    @eval ^(a::$T, x::T) where {T<:Complex} = exp( x*log(a) )
end


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
function ^(a::Taylor1, r::S) where {S<:Real}
    r == zero(r) && return one(a)
    r == one(r)/2 && return sqrt(a)
    isinteger(r) && return a^round(Int,r)

    # First non-zero coefficient
    l0nz = findfirst(a)
    l0nz < 0 && return zero(a)

    # The first non-zero coefficient of the result; must be integer
    !isinteger(r*l0nz) && throw(ArgumentError(
        """The 0th order Taylor1 coefficient must be non-zero
        to raise the Taylor1 polynomial to a non-integer exponent."""))
    lnull = trunc(Int, r*l0nz )

    @inbounds aux = ( a[l0nz] )^r
    k0 = lnull+l0nz
    c = Taylor1( zero(aux), a.order)
    @inbounds c[lnull] = aux
    for k = k0+1:a.order+l0nz
        pow!(c, a, r, k, l0nz)
    end

    return c
end

## Real power ##
function ^(a::TaylorN, r::S) where {S<:Real}
    r == zero(r) && return TaylorN( one(eltype(a)), 0 )
    r == one(r)/2 && return sqrt(a)
    isinteger(r) && return a^round(Int,r)

    a0 = constant_term(a)
    @assert a0 != zero(a0)

    c = TaylorN( a0^r, a.order)
    for ord in 1:a.order
        pow!(c, a, r, ord)
    end

    return c
end


# Homogeneous coefficients for real power
doc"""
    pow!(c, a, r::Real, k::Int, k0::Int=0)

Update the `k-th` expansion coefficient `c[k]` of `c = a^r`, for
both `c` and `a` either `Taylor1` or `TaylorN`.

The coefficients are given by

```math
c_k = \frac{1}{k a_0} \sum_{j=0}^{k-1} \big(r(k-j) -j\big)a_{k-j} c_j.
```

For `Taylor1` polynomials, `k0` is the order of the first non-zero
coefficient of `a`.

"""
@inline function pow!(c::Taylor1, a::Taylor1, r::S, k::Int, l0::Int=0) where
        {S<:Real}

    if r == 0
        return one!(c, a, k)
    elseif r == 1
        return identity!(c, a, k)
    elseif r == 2
        return sqr!(c, a, k)
    elseif r == 0.5
        return sqrt!(c, a, k)
    end

    if k == l0
        @inbounds c[0] = ( a[l0] )^r
        return nothing
    end

    for i = 0:k-l0-1
        k-i > a.order && continue
        aux = r*(k-i) - i
        @inbounds c[k-l0] += aux * a[k-i] * c[i]
    end
    aux = k - l0*(r+1)
    @inbounds c[k-l0] = c[k-l0] / (aux * a[l0])

    return nothing
end

@inline function pow!(c::TaylorN, a::TaylorN, r::S, k::Int) where {S<:Real}

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

    for i = 0:k-1
        aux = r*(k-i) - i
        mul!(c[k], aux*a[k-i], c[i])
    end
    @inbounds c[k] = c[k] / (k * constant_term(a))

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
        for k = 1:a.order
            sqr!(c, a, k)
        end
        return c
    end
end

function square(a::HomogeneousPolynomial)
    T = eltype(a)

    order = 2*a.order
    order > get_order() && return HomogeneousPolynomial([zero(T)], get_order())

    res = HomogeneousPolynomial([zero(T)], order)
    sqr!(res, a)
    return res
end


# Homogeneous coefficients for square
doc"""
    sqr!(c, a, k::Int) --> nothing

Update the `k-th` expansion coefficient `c[k]` of `c = a^2`, for
both `c` and `a` either `Taylor1` or `TaylorN`.

The coefficients are given by

```math
\begin{eqnarray*}
c_k & = & 2 \sum_{j=0}^{(k-1)/2} a_{k-j} a_j,
    \text{ if k is odd,} \\
c_k & = & 2 \sum_{j=0}^{(k-2)/2} a_{k-j} a_j + (a_{k/2})^2,
    \text{ if k is even. }
\end{eqnarray*}
```

""" sqr!

for T = (:Taylor1, :TaylorN)
    @eval begin
        @inline function sqr!(c::$T, a::$T, k::Int)
            if k == 0
                @inbounds c[0] = constant_term(a)^2
                return nothing
            end

            kodd = k%2
            kend = div(k - 2 + kodd, 2)
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
                @inbounds c[k] += a[div(k,2)]^2
            else
                sqr!(c[k], a[div(k,2)])
            end

            return nothing
        end
    end
end


"""
    sqr!(c, a)

Return `c = a*a` with no allocation; all parameters are `HomogeneousPolynomial`.

"""
@inline function sqr!(c::HomogeneousPolynomial, a::HomogeneousPolynomial)
    iszero(a) && return nothing

    T = eltype(c)
    @inbounds num_coeffs_a = size_table[a.order+1]

    @inbounds posTb = pos_table[c.order+1]
    @inbounds idxTb = index_table[a.order+1]

    @inbounds for na = 1:num_coeffs_a
        ca = a[na]
        ca == zero(T) && continue
        inda = idxTb[na]
        pos = posTb[2*inda]
        c[pos] += ca * ca
        @inbounds for nb = na+1:num_coeffs_a
            cb = a[nb]
            cb == zero(T) && continue
            indb = idxTb[nb]
            pos = posTb[inda+indb]
            c[pos] += 2 * ca * cb
        end
    end

    return nothing
end



## Square root ##
function sqrt(a::Taylor1)

    # First non-zero coefficient
    l0nz = findfirst(a)
    if l0nz < 0
        return zero(a)
    elseif l0nz%2 == 1 # l0nz must be pair
        throw(ArgumentError(
        """First non-vanishing Taylor1 coefficient must correspond
        to an **even power** in order to expand `sqrt` around 0."""))
    end

    # The last l0nz coefficients are set to zero.
    lnull = div(l0nz, 2)
    @inbounds aux = sqrt( a[l0nz] )
    T = typeof(aux)

    c = Taylor1( zero(T), a.order )
    @inbounds c[lnull] = aux
    for k = lnull+1:a.order
        sqrt!(c, a, k, lnull)
    end

    return c
end

function sqrt(a::TaylorN)
    @inbounds p0 = sqrt( constant_term(a) )
    if p0 == zero(p0)
        throw(ArgumentError(
        """The 0-th order TaylorN coefficient must be non-zero
        in order to expand `sqrt` around 0."""))
    end

    c = TaylorN( p0, a.order)
    for k = 1:a.order
        sqrt!(c, a, k)
    end

    return c
end

# Homogeneous coefficients for the square-root
doc"""
    sqrt!(c, a, k::Int, k0::Int=0)

Compute the `k-th` expansion coefficient `c[k]` of `c = sqrt(a)`
for both`c` and `a` either `Taylor1` or `TaylorN`.

The coefficients are given by

```math
\begin{eqnarray*}
c_k &=& \frac{1}{2 c_0} \big( a_k - 2 \sum_{j=1}^{(k-1)/2} c_{k-j}c_j\big),
    \text{ if $k$ is odd,} \\
c_k &=& \frac{1}{2 c_0} \big( a_k - 2 \sum_{j=1}^{(k-2)/2} c_{k-j}c_j
    - (c_{k/2})^2\big), \text{ if $k$ is even.}
\end{eqnarray*}
```

For `Taylor1` polynomials, `k0` is the order of the first non-zero
coefficient, which must be even.

"""
@inline function sqrt!(c::Taylor1, a::Taylor1, k::Int, k0::Int=0)

    if k == k0
        @inbounds c[k] = sqrt(a[2*k0])
        return nothing
    end

    kodd = (k - k0)%2
    kend = div(k - k0 - 2 + kodd, 2)
    @inbounds for i = k0+1:k0+kend
        c[k] += c[i] * c[k+k0-i]
    end
    @inbounds aux = a[k+k0] - 2*c[k]
    if kodd == 0
        @inbounds aux = aux - (c[kend+k0+1])^2
    end
    @inbounds c[k] = aux / (2*c[k0])

    return nothing
end

@inline function sqrt!(c::TaylorN, a::TaylorN, k::Int)

    if k == 0
        @inbounds c[0] = sqrt( constant_term(a) )
        return nothing
    end

    kodd = k%2
    kend = div(k - 2 + kodd, 2)
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
