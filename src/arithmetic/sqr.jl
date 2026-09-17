# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#


## Square ##
"""
    square(a::AbstractSeries) --> typeof(a)

Return `a^2`; see [`TaylorSeries.sqr!`](@ref).
""" square

function square(a::Taylor1{T}) where {T}
    c = zero(a)
    aux = zero(a[0])
    for k in eachindex(a)
        sqr!(c, a, aux, k)
    end
    return c
end
function square(a::TaylorN{T}) where {T}
    c = zero(a)
    aux = zero(a[0][1])
    for k in eachindex(a)
        sqr!(c, a, aux, k)
    end
    return c
end

function square(a::HomogeneousPolynomial)
    order = 2*TS.order(a)
    # NOTE: the following returns order 0, but could be TS.order(), or TS.order(a)
    order > TS.order(a.space) && return HomogeneousPolynomial(a.space, zero(a[1]), 0)
    res = HomogeneousPolynomial(a.space, zero(a[1]), order)
    accsqr!(res, a)
    return res
end

#auxiliary function to avoid allocations
function sqr_orderzero!(c::Taylor1{T}, a::Taylor1{T}) where {T<:NumberNotSeries}
    @inbounds c[0] = a[0]^2
    return nothing
end
function sqr_orderzero!(c::TaylorN{T}, a::TaylorN{T}) where {T<:NumberNotSeries}
    @inbounds c[0][1] = a[0][1]^2
    return nothing
end
function sqr_orderzero!(c::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}}) where
        {T<:NumberNotSeries}
    aux = zero(a[0])
    @inbounds for ord in eachindex(c[0])
        sqr!(c[0], a[0], aux, ord)
    end
    return nothing
end
function sqr_orderzero!(c::TaylorN{Taylor1{T}}, a::TaylorN{Taylor1{T}}) where
        {T<:NumberNotSeries}
    aux = zero(a[0][1][0])
    @inbounds for ord in eachindex(c[0][1])
        sqr!(c[0][1], a[0][1], aux, ord)
    end
    return nothing
end
function sqr_orderzero!(c::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}}) where
        {T<:Number}
    aux = zero(a[0][0])
    @inbounds for ord in eachindex(c[0])
        sqr!(c[0], a[0], aux, ord)
    end
    return nothing
end
# function sqr_orderzero!(c::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}}) where
#         {T<:NumberNotSeries}
#     aux = zero(a[0][0])
#     @inbounds for ord in eachindex(c[0])
#         sqr!(c[0], a[0], aux, ord)
#     end
#     return nothing
# end

# Homogeneous coefficients for square
@doc doc"""
    sqr!(c, a, aux, k::Int) --> nothing

Update the `k-th` expansion coefficient `c[k]` of `c = a^2`, for
both `c` and `a` either `Taylor1{T}` or `TaylorN{T}`; `aux::T`
is an auxiliary.

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

function sqr!(c::Taylor1{T}, a::Taylor1{T}, ::T, k::Int) where {T<:Number}
    if k == 0
        sqr_orderzero!(c, a)
        return nothing
    end
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    kk = k+1
    @inbounds acc = zero(c_coeffs[kk])
    # Recursion formula
    kodd = k%2
    kend = (k - 2 + kodd) >> 1
    @inbounds for i = 1:kend+1
        acc += a_coeffs[i] * a_coeffs[kk-i+1]
    end
    acc = 2 * acc
    if kodd == 0
        @inbounds acc += a_coeffs[(k >> 1)+1]^2
    end
    @inbounds c_coeffs[kk] = acc
    return nothing
end

function sqr!(c::TaylorN{T}, a::TaylorN{T}, ::T, k::Int) where {T<:Number}
    if k == 0
        sqr_orderzero!(c, a)
        return nothing
    end
    # Sanity
    zero!(c, k)
    # Recursion formula
    kodd = k%2
    kend = (k - 2 + kodd) >> 1
    @inbounds for i = 0:kend
        mul!(c[k], a[i], a[k-i])
    end
    @inbounds mul!(c, 2, c, k)
    kodd == 1 && return nothing
    accsqr!(c[k], a[k >> 1])
    return nothing
end

# in-place squaring: given `c`, compute expansion of `c^2` and save back into `c`
function sqr!(c::Taylor1{T}, k::Int) where {T<:NumberNotSeries}
    if k == 0
        sqr_orderzero!(c, c)
        return nothing
    end
    c_coeffs = c.coeffs
    kk = k+1
    # Recursion formula
    kodd = k%2
    kend = (k - 2 + kodd) >> 1
    @inbounds acc = zero(c_coeffs[kk])
    (kend >= 0) && ( @inbounds acc = c_coeffs[1] * c_coeffs[kk] )
    @inbounds for i = 2:kend+1
        acc += c_coeffs[i] * c_coeffs[kk-i+1]
    end
    acc = 2 * acc
    (kodd == 0) && ( @inbounds acc += c_coeffs[(k >> 1)+1]^2 )
    @inbounds c_coeffs[kk] = acc
    return nothing
end

function sqr!(c::TaylorN{T}, k::Int) where {T<:NumberNotSeries}
    if k == 0
        sqr_orderzero!(c, c)
        return nothing
    end
    # Recursion formula
    kodd = k%2
    kend = (k - 2 + kodd) >> 1
    (kend >= 0) && ( @inbounds mul!(c, c[0][1], c, k) )
    @inbounds for i = 1:kend
        mul!(c[k], c[i], c[k-i])
    end
    @inbounds mul!(c, 2, c, k)
    if kodd == 0
        accsqr!(c[k], c[k >> 1])
    end
    return nothing
end

function sqr!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}}, aux::TaylorN{T},
        ordT::Int) where {T<:NumberNotSeries}
    # Sanity
    zero!(res, ordT)
    if ordT == 0
        @inbounds for ordQ in eachindex(a[0])
            @inbounds sqr!(res[0], a[0], aux[0][1], ordQ)
        end
        return nothing
    end
    # Recursion formula
    kodd = ordT%2
    kend = (ordT - 2 + kodd) >> 1
    zero!(aux)
    (kodd == 0) && @inbounds for ordQ in eachindex(a[0])
        sqr!(res[ordT], a[ordT >> 1], aux[0][1], ordQ)
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

function sqr!(c::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}}, aux::Taylor1{T},
        k::Int) where {T<:NumberNotSeries}
    if k == 0
        sqr_orderzero!(c, a)
        return nothing
    end
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    kk = k+1
    @inbounds c_k = c_coeffs[kk]
    zero!(c_k)
    # Recursion formula
    kodd = k%2
    kend = (k - 2 + kodd) >> 1
    @inbounds for i = 0:kend
        mul_scalar!(aux, 2, a_coeffs[i+1], a_coeffs[kk-i])
        add!(c_k, c_k, aux)
    end
    kodd == 1 && return nothing
    # c[k] += a[k >> 1]^2
    aaux = zero(aux[0])
    zero!(aux)
    @inbounds a_mid = a_coeffs[(k >> 1)+1]
    @inbounds for j in eachindex(a_mid)
        sqr!(aux, a_mid, aaux, j)
    end
    add!(c_k, c_k, aux)
    return nothing
end

function sqr!(c::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}}, aux::Taylor1{T},
        k::Int) where {T<:Number}
    if k == 0
        sqr_orderzero!(c, a)
        return nothing
    end
    # Sanity
    zero!(c[k])
    zero!(aux)
    # Recursion formula
    kodd = k%2
    kend = (k - 2 + kodd) >> 1
    @inbounds for i = 0:kend
        for j in eachindex(a[k])
            # c[k] += 2 * a[i] * a[k-i]
            mul_scalar!(aux, 2, a[i], a[k-i], j)
            add!(c[k], c[k], aux, j)
        end
    end
    kodd == 1 && return nothing
    # @inbounds c[k] += a[k >> 1]^2
    aaux = zero(aux[0])
    for j in eachindex(a[k])
        zero!(aux, j)
        sqr!(aux, a[k >> 1], aaux, j)
        add!(c[k], c[k], aux, j)
    end
    return nothing
end


"""
    accsqr!(c, a)

Returns `c += a*a` with no allocation; all parameters are `HomogeneousPolynomial`.

"""
function accsqr!(c::HomogeneousPolynomial{T}, a::HomogeneousPolynomial{T}) where
        {T<:NumberNotSeriesN}
    _check_same_space(c, a)
    iszero(a) && return nothing

    sp = c.space
    degree_a = order(a)
    _check_homogeneous_product_order(c, a, a)
    degree_a == 0 && return _muladd_scalar_unchecked!(c, a[1], a)
    order_a = degree_a+1
    @inbounds num_coeffs_a = sp.size_table[order_a]
    input_positions = _product_table(sp, degree_a, degree_a).input_positions

    @inbounds for na = 1:num_coeffs_a
        ca = a[na]
        _isthinzero(ca) && continue
        pos = input_positions[(na-1)*num_coeffs_a + na]
        c[pos] += ca^2
        @inbounds for nb = na+1:num_coeffs_a
            cb = a[nb]
            _isthinzero(cb) && continue
            pos = input_positions[(na-1)*num_coeffs_a + nb]
            c[pos] += 2 * ca * cb
        end
    end

    return nothing
end
