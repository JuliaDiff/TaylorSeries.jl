# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

## Multiplication ##
for T in (:Taylor1, :TaylorN)
    @eval begin
        function *(a::T, b::$T{S}) where {T<:NumberNotSeries, S<:NumberNotSeries}
            v = $T( a * b[0], order(b))
            @inbounds for k in eachindex(v)
                mul!(v, b, a, k)
            end
            return v
        end
        *(b::$T{S}, a::T) where {T<:NumberNotSeries, S<:NumberNotSeries} = a * b
        function *(a::T, b::$T{T}) where {T<:NumberNotSeries}
            v = $T( a * b[0], order(b))
            @inbounds for k in eachindex(v)
                mul!(v, b, a, k)
            end
            return v
        end
        *(b::$T{T}, a::T) where {T<:NumberNotSeries} = a * b
    end
end

*(a::T, b::HomogeneousPolynomial{S}) where {T<:NumberNotSeries,
    S<:NumberNotSeries} = HomogeneousPolynomial(b.space, a * b.coeffs, order(b))
*(b::HomogeneousPolynomial{S}, a::T) where {T<:NumberNotSeries,
    S<:NumberNotSeries} = a * b
*(a::T, b::HomogeneousPolynomial{T}) where {T<:NumberNotSeries} =
    HomogeneousPolynomial(b.space, a * b.coeffs, order(b))
*(b::HomogeneousPolynomial{T}, a::T) where {T<:NumberNotSeries} = a * b

for T in (:HomogeneousPolynomial, :TaylorN)
    @eval begin
        *(a::T, b::$T{Taylor1{S}}) where {T<:NumberNotSeries,
            S<:NumberNotSeries} = $T( a .* b.coeffs, order(b))
        *(b::$T{Taylor1{S}}, a::T) where {T<:NumberNotSeries,
            S<:NumberNotSeries} = a * b
        *(a::T, b::Taylor1{$T{S}}) where {T<:NumberNotSeries,
            S<:NumberNotSeries} = Taylor1(a .* b.coeffs)
        *(b::Taylor1{$T{S}}, a::T) where
            {T<:NumberNotSeries, S<:NumberNotSeries} = a * b
        *(a::Taylor1{T}, b::$T{Taylor1{S}}) where
            {T<:NumberNotSeries, S<:NumberNotSeries} = $T(a .* b.coeffs, order(b))
        *(b::$T{Taylor1{R}}, a::Taylor1{T}) where
            {T<:NumberNotSeries, R<:NumberNotSeries} = a * b
        *(a::$T{T}, b::Taylor1{$T{S}}) where {T<:NumberNotSeries,
            S<:NumberNotSeries} = Taylor1(a .* b.coeffs)
        *(b::Taylor1{$T{S}}, a::$T{T}) where {T<:NumberNotSeries,
            S<:NumberNotSeries} = a * b
    end
end

function *(a::Taylor1{T}, b::Taylor1{T}) where {T<:Number}
    if order(a) != order(b)
        a, b = fixorder(a, b)
    end
    c = zero(a)
    for ord in eachindex(c)
        _muladd_unchecked!(c, a, b, ord) # updates c[ord]
    end
    return c
end

function *(a::TaylorN{T}, b::TaylorN{T}) where {T<:NumberNotSeriesN}
    _check_same_space(a, b)
    if order(a) != order(b)
        a, b = fixorder(a, b)
    end
    c = zero(a)
    for ord in eachindex(c)
        _muladd_unchecked!(c, a, b, ord) # updates c[ord]
    end
    return c
end

function *(a::T, b::Taylor1{Taylor1{T}}) where {T<:NumberNotSeriesN}
    v = Taylor1( a * b[0], order(b))
    @inbounds for k in eachindex(v)
        mul!(v, b, a, k)
    end
    return v
end
*(b::Taylor1{Taylor1{T}}, a::T) where {T<:NumberNotSeriesN} = a * b

function *(a::Taylor1{T}, b::Taylor1{Taylor1{T}}) where {T<:NumberNotSeriesN}
    v = Taylor1( a * b[0], order(b))
    @inbounds for k in eachindex(v)
        mul!(v, b, a, k)
    end
    return v
end
*(b::Taylor1{Taylor1{T}}, a::Taylor1{T}) where {T<:NumberNotSeriesN} = a * b


*(a::HomogeneousPolynomial{T}, b::HomogeneousPolynomial{S}) where
    {T<:NumberNotSeriesN,S<:NumberNotSeriesN} = *(promote(a,b)...)

function *(a::HomogeneousPolynomial{T}, b::HomogeneousPolynomial{T}) where
        {T<:NumberNotSeriesN}
    _check_same_space(a, b)
    order = TS.order(a) + TS.order(b)
    # NOTE: the following returns order 0, but could be TS.order(), or TS.order(a)
    order > TS.order(a.space) && return HomogeneousPolynomial(a.space, zero(a[1]), TS.order(a))
    res = HomogeneousPolynomial(a.space, zero(a[1]), order)
    mul!(res, a, b)
    return res
end

function *(a::Taylor1{TaylorN{T}}, b::Taylor1{TaylorN{S}}) where
        {T<:NumberNotSeries, S<:NumberNotSeries}
    R = promote_type(T,S)
    return *(convert(Taylor1{TaylorN{R}}, a), convert(Taylor1{TaylorN{R}}, b))
end

function *(a::Taylor1{TaylorN{T}}, b::Taylor1{TaylorN{T}}) where {T<:NumberNotSeries}
    _check_same_space(a[0], b[0])
    if (order(a) != order(b)) || any(order.(a.coeffs) .!= order.(b.coeffs))
        a, b = fixorder(a, b)
    end
    res = zero(a)
    for ordT in eachindex(a)
        _mul_unchecked!(res, a, b, ordT)
    end
    return res
end


# Internal multiplication functions
function mul!(c::Taylor1{T}, a::Taylor1{T}, b::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    kk = k+1
    @inbounds acc = zero(c_coeffs[kk])
    @inbounds for i = 1:kk
        acc += a_coeffs[i] * b_coeffs[kk-i+1]
    end
    @inbounds c_coeffs[kk] = acc
    return nothing
end
function mul!(v::Taylor1{T}, a::Taylor1{S}, b::NumberNotSeries, k::Int) where
        {T<:NumberNotSeries, S<:NumberNotSeries}
    @inbounds v.coeffs[k+1] = a.coeffs[k+1] * b
    return nothing
end
mul!(v::Taylor1{T}, a::NumberNotSeries, b::Taylor1{S}, k::Int) where
        {T<:NumberNotSeries, S<:NumberNotSeries} = mul!(v, b, a, k)

function mul!(v::Taylor1{T}, a::Taylor1{S}, b::NumberNotSeries) where
        {T<:NumberNotSeries, S<:NumberNotSeries}
    v_coeffs = v.coeffs
    a_coeffs = a.coeffs
    @inbounds for i in eachindex(v_coeffs)
        v_coeffs[i] = a_coeffs[i] * b
    end
    return nothing
end
mul!(v::Taylor1{T}, a::NumberNotSeries, b::Taylor1{S}) where
        {T<:NumberNotSeries, S<:NumberNotSeries} = mul!(v, b, a)
#
function muladd!(c::Taylor1{T}, a::Taylor1{T}, b::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    kk = k+1
    @inbounds acc = c_coeffs[kk]
    @inbounds for i = 1:kk
        acc += a_coeffs[i] * b_coeffs[kk-i+1]
    end
    @inbounds c_coeffs[kk] = acc
    return nothing
end

function muladd!(c::Taylor1{T}, a::Taylor1{T}, b::Taylor1{T}) where
        {T<:NumberNotSeries}
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    @inbounds for kk in eachindex(c_coeffs)
        acc = c_coeffs[kk]
        for i = 1:kk
            acc += a_coeffs[i] * b_coeffs[kk-i+1]
        end
        c_coeffs[kk] = acc
    end
    return nothing
end

@inline function _muladd_unchecked!(c::Taylor1{T}, a::Taylor1{T},
        b::Taylor1{T}, k::Int) where {T<:Number}
    mul!(c, a, b, k)
    return nothing
end

# function muladd!(v::Taylor1{T}, a::Taylor1{T}, b::NumberNotSeries, k::Int) where
#         {T<:NumberNotSeries}
#     @inbounds v[k] += a[k] * b
#     return nothing
# end
# muladd!(v::Taylor1{T}, a::NumberNotSeries, b::Taylor1{T}, k::Int) where
#         {T<:NumberNotSeries} = muladd!(v, b, a, k)
# Implements c[k] = scalar \sum_i a[i] b[k-i]
function mul_scalar!(c::Taylor1{T}, scalar::NumberNotSeries, a::Taylor1{T},
        b::Taylor1{T}, k::Int) where {T<:NumberNotSeries}
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    kk = k+1
    @inbounds acc = zero(c_coeffs[kk])
    @inbounds for i = 1:kk
        acc += a_coeffs[i] * b_coeffs[kk-i+1]
    end
    @inbounds c_coeffs[kk] = scalar * acc
    return nothing
end

function mul_scalar!(c::Taylor1{T}, scalar::NumberNotSeries, a::Taylor1{T},
        b::Taylor1{T}) where {T<:NumberNotSeries}
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    @inbounds for kk in eachindex(c_coeffs)
        acc = zero(c_coeffs[kk])
        for i = 1:kk
            acc += a_coeffs[i] * b_coeffs[kk-i+1]
        end
        c_coeffs[kk] = scalar * acc
    end
    return nothing
end

# NOTE: For TaylorN, `mul!` (`muladd!`) *accumulates* the result of a * b in c[k]
mul!(c::TaylorN{T}, a::TaylorN{T}, b::TaylorN{T}, k::Int) where
        {T<:Number} = muladd!(c, a, b, k)
function mul!(v::TaylorN, a::TaylorN, b::NumberNotSeries, k::Int)
    _check_same_space(v, a)
    @inbounds for i in eachindex(v[k])
        v[k][i] = a[k][i] * b
    end
    return nothing
end
mul!(v::TaylorN{T}, a::NumberNotSeries, b::TaylorN{T}, k::Int) where
    {T<:Number} = mul!(v, b, a, k)
#
function muladd!(c::TaylorN{T}, a::TaylorN{T}, b::TaylorN{T},
        k::Int) where {T<:Number}
    _check_same_space(c, a, b)
    _muladd_unchecked!(c, a, b, k)
    return nothing
end
function muladd!(v::TaylorN, a::TaylorN, b::NumberNotSeries, k::Int)
    _check_same_space(v, a)
    @inbounds for i in eachindex(v[k])
        v[k][i] += a[k][i] * b
    end
    return nothing
end
muladd!(v::TaylorN, a::NumberNotSeries, b::TaylorN, k::Int) = muladd!(v, b, a, k)
function mul_scalar!(c::TaylorN{T}, scalar::NumberNotSeries, a::TaylorN{T},
        b::TaylorN{T}, k::Int) where {T<:Number}
    _check_same_space(c, a, b)
    _mul_scalar_unchecked!(c, scalar, a, b, k)
    return nothing
end

# Nested Taylor1s
function mul!(c::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}}, b::Taylor1{Taylor1{T}},
        k::Int) where {T<:NumberNotSeries}
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    kk = k+1
    @inbounds c_k = c_coeffs[kk]
    zero!(c_k)
    @inbounds for i = 1:kk
        muladd!(c_k, a_coeffs[i], b_coeffs[kk-i+1])
    end
    return nothing
end

function mul!(c::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}}, b::Taylor1{Taylor1{T}},
        k::Int) where {T<:NumberNotSeriesN}
    @inbounds for j in eachindex(c[k])
        zero!(c[k], j)
        for i = 0:k
            muladd!(c[k], a[i], b[k-i], j)
        end
    end
    return nothing
end

mul!(v::Taylor1{Taylor1{T}}, a::Taylor1{T}, b::Taylor1{Taylor1{T}}, k::Int) where
        {T<:NumberNotSeriesN} = mul!(v, b, a, k)
function mul!(v::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}}, b::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    @inbounds v_k = v.coeffs[k+1]
    @inbounds a_k = a.coeffs[k+1]
    if v_k === a_k
        mul!(v_k, b)
    else
        mul!(v_k, a_k, b)
    end
    return nothing
end

function mul!(v::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}}, b::Taylor1{T}, k::Int) where
        {T<:NumberNotSeriesN}
    @inbounds for i in eachindex(v[k])
        mul!(v[k], a[k], b, i)
    end
    return nothing
end

mul!(v::Taylor1{Taylor1{T}}, a::NumberNotSeries, b::Taylor1{Taylor1{T}}, k::Int) where
        {T<:NumberNotSeriesN} = mul!(v, b, a, k)
function mul!(v::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}}, b::NumberNotSeries,
        k::Int) where {T<:NumberNotSeries}
    @inbounds mul!(v.coeffs[k+1], a.coeffs[k+1], b)
    return nothing
end

function mul!(v::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}}, b::NumberNotSeries,
        k::Int) where {T<:NumberNotSeriesN}
    @inbounds for i in eachindex(v[k])
        mul!(v[k], a[k], b, i)
    end
    return nothing
end

function muladd!(c::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}},
        b::Taylor1{Taylor1{T}}, k::Int) where {T<:NumberNotSeries}
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    kk = k+1
    @inbounds c_k = c_coeffs[kk]
    @inbounds for i = 1:kk
        muladd!(c_k, a_coeffs[i], b_coeffs[kk-i+1])
    end
    return nothing
end

function muladd!(c::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}},
        b::Taylor1{Taylor1{T}}, k::Int) where {T<:NumberNotSeriesN}
    @inbounds for j in eachindex(c[k])
        for i = 0:k
            muladd!(c[k], a[i], b[k-i], j)
        end
    end
    return nothing
end

# muladd!(v::Taylor1{Taylor1{T}}, a::Taylor1{T}, b::Taylor1{Taylor1{T}}, k::Int) where
#         {T<:NumberNotSeriesN} = muladd!(v, b, a, k)
# function muladd!(v::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}}, b::Taylor1{T}, k::Int) where
#         {T<:NumberNotSeriesN}
#     @inbounds for i in eachindex(v[k])
#         muladd!(v[k], a[k], b, i)
#     end
#     return nothing
# end
# function muladd!(v::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}}, b::T, k::Int) where
#         {T<:NumberNotSeriesN}
#     @inbounds for i in eachindex(v[k])
#         muladd!(v[k], a[k], b, i)
#     end
#     return nothing
# end
# muladd!(v::Taylor1{Taylor1{T}}, a::T, b::Taylor1{Taylor1{T}}, k::Int) where
#     {T<:NumberNotSeriesN} = muladd!(v, b, a, k)
function mul_scalar!(c::Taylor1{Taylor1{T}}, scalar::NumberNotSeries,
        a::Taylor1{Taylor1{T}}, b::Taylor1{Taylor1{T}}, k::Int) where
        {T<:NumberNotSeries}
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    kk = k+1
    @inbounds c_k = c_coeffs[kk]
    zero!(c_k)
    @inbounds for i = 1:kk
        muladd!(c_k, a_coeffs[i], b_coeffs[kk-i+1])
    end
    c_k_coeffs = c_k.coeffs
    @inbounds for i in eachindex(c_k_coeffs)
        c_k_coeffs[i] *= scalar
    end
    return nothing
end

function mul_scalar!(c::Taylor1{Taylor1{T}}, scalar::NumberNotSeries,
        a::Taylor1{Taylor1{T}}, b::Taylor1{Taylor1{T}}, k::Int) where {T<:Number}
    mul!(c, a, b, k)
    # c[k] <- scalar * c[k]
    for ord in eachindex(c[k])
        mul!(c[k], c[k], scalar, ord)
    end
    return nothing
end

function mul_scalar!(c::Taylor1{Taylor1{T}}, scalar::NumberNotSeries,
        a::Taylor1{Taylor1{T}}, b::Taylor1{Taylor1{T}}) where
        {T<:NumberNotSeries}
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    @inbounds for kk in eachindex(c_coeffs)
        c_k = c_coeffs[kk]
        zero!(c_k)
        for i = 1:kk
            muladd!(c_k, a_coeffs[i], b_coeffs[kk-i+1])
        end
        c_k_coeffs = c_k.coeffs
        for i in eachindex(c_k_coeffs)
            c_k_coeffs[i] *= scalar
        end
    end
    return nothing
end


# for T in (:Taylor1, :TaylorN)
#     @eval begin
#         function mul!(v::$T{T}, a::$T{T}, b::NumberNotSeries) where {T<:Number}
#             for k in eachindex(v)
#                 mul!(v, a, b, k)
#             end
#             return nothing
#         end
#         mul!(v::$T{T}, a::NumberNotSeries, b::$T{T}) where {T<:Number} = mul!(v, b, a)
#     end
# end

# in-place product: `a` <- `a*b`
# this method computes the product `a*b` and saves it back into `a`
# assumes `a` and `b` are of same order
function mul!(a::TaylorN{T}, b::TaylorN{T}) where {T<:Number}
    @inbounds for k in reverse(eachindex(a))
        mul!(a, a, b[0][1], k)
        for l in 1:k
            mul!(a[k], a[k-l], b[l])
        end
    end
    return nothing
end
function mul!(a::Taylor1{T}, b::Taylor1{T}) where {T<:NumberNotSeries}
    @inbounds for k in reverse(eachindex(a))
        # a[k] <- a[k]*b[0]
        mul!(a, a, b[0], k)
        for l in 1:k
            # a[k] <- a[k] + a[k-l] * b[l]
            a[k] += a[k-l] * b[l]
        end
    end
    return nothing
end
function mul!(a::Taylor1{TaylorN{T}}, b::Taylor1{TaylorN{T}}) where
        {T<:NumberNotSeries}
    @inbounds for k in reverse(eachindex(a))
        mul!(a, a, b[0], k)
        for l in 1:k
            # a[k] += a[k-l] * b[l]
            for m in eachindex(a[k])
                mul!(a[k], a[k-l], b[l], m)
            end
        end
    end
    return nothing
end
function mul!(a::Taylor1{Taylor1{T}}, b::Taylor1{Taylor1{T}}) where
        {T<:NumberNotSeries}
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    @inbounds for kk in reverse(eachindex(a_coeffs))
        a_k = a_coeffs[kk]
        mul!(a_k, b_coeffs[1])
        for l = 2:kk
            muladd!(a_k, a_coeffs[kk-l+1], b_coeffs[l])
        end
    end
    return nothing
end

function mul!(a::Taylor1{Taylor1{T}}, b::Taylor1{Taylor1{T}}) where
        {T<:NumberNotSeriesN}
    @inbounds for k in reverse(eachindex(a))
        # a[k] <- a[k]*b[0]
        mul!(a, a, b[0], k)
        for l in 1:k
            # a[k] <- a[k] + a[k-l] * b[l]
            for m in eachindex(a[k])
                muladd!(a[k], a[k-l], b[l], m)
            end
        end
    end
    return nothing
end

function _mul_unchecked!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        b::Taylor1{TaylorN{T}}, ordT::Int) where {T<:NumberNotSeries}
    zero!(res, ordT)
    for k in 0:ordT
        @inbounds for ordQ in eachindex(a[ordT])
            _muladd_unchecked!(res[ordT], a[k], b[ordT-k], ordQ)
        end
    end
    return nothing
end

function mul!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        b::Taylor1{TaylorN{T}}, ordT::Int) where {T<:NumberNotSeries}
    _check_same_space(res[0], a[0], b[0])
    _mul_unchecked!(res, a, b, ordT)
    return nothing
end

function mul!(res::Taylor1{TaylorN{T}}, a::NumberNotSeries,
        b::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
    res_k = res.coeffs[k+1]
    b_k = b.coeffs[k+1]
    res_hps = res_k.coeffs
    b_hps = b_k.coeffs
    @inbounds for l in eachindex(res_hps)
        res_hp = res_hps[l].coeffs
        b_hp = b_hps[l].coeffs
        for m in eachindex(res_hp)
            res_hp[m] = a*b_hp[m]
        end
    end
    return nothing
end
mul!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}}, b::NumberNotSeries,
    k::Int) where {T<:NumberNotSeries} = mul!(res, b, a, k)


# in-place product (assumes equal order)
function mul!(c::Taylor1{T}, a::Taylor1{T}, b::Taylor1{T}) where
        {T<:NumberNotSeries}
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    @inbounds for kk in eachindex(c_coeffs)
        acc = zero(c_coeffs[kk])
        for i = 1:kk
            acc += a_coeffs[i] * b_coeffs[kk-i+1]
        end
        c_coeffs[kk] = acc
    end
    return nothing
end

# Fallback for nested coefficient types; scalar Taylor1 has a coefficient-vector
# specialization above.
function mul!(c::Taylor1{T}, a::Taylor1{T}, b::Taylor1{T}) where {T<:Number}
    for k in eachindex(c)
        mul!(c, a, b, k)
    end
end

function mul!(c::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}},
        b::Taylor1{Taylor1{T}}) where {T<:NumberNotSeries}
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    @inbounds for kk in eachindex(c_coeffs)
        c_k = c_coeffs[kk]
        zero!(c_k)
        for i = 1:kk
            muladd!(c_k, a_coeffs[i], b_coeffs[kk-i+1])
        end
    end
    return nothing
end

function mul!(c::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        b::Taylor1{TaylorN{T}}) where {T<:NumberNotSeries}
    _check_same_space(c[0], a[0], b[0])
    for k in eachindex(c)
        _mul_unchecked!(c, a, b, k)
    end
    return nothing
end

function mul!(c::TaylorN{T}, a::TaylorN{T}, b::TaylorN{T}) where {T<:NumberNotSeriesN}
    _check_same_space(c, a, b)
    for k in eachindex(c)
        _muladd_unchecked!(c, a, b, k)
    end
end

# function muladd!(c::Taylor1{T}, a::Taylor1{T}, b::Taylor1{T}) where {T<:Number}
#     for k in eachindex(c)
#         muladd!(c, a, b, k)
#     end
# end
#
# function mul_scalar!(c::Taylor1{T}, scalar::NumberNotSeries, a::Taylor1{T},
#         b::Taylor1{T}) where {T<:Number}
#     for k in eachindex(c)
#         mul_scalar!(c, scalar, a, b, k)
#     end
# end

function mul_scalar!(c::TaylorN{T}, scalar::NumberNotSeries, a::TaylorN{T},
        b::TaylorN{T}) where {T<:NumberNotSeriesN}
    _check_same_space(c, a, b)
    for k in eachindex(c)
        _mul_scalar_unchecked!(c, scalar, a, b, k)
    end
end


@doc doc"""
    mul!(c, a, b, k::Int) --> nothing

Update the `k`-th expansion coefficient `c[k]` of `c = a * b`,
where all `c`, `a`, and `b` are either `Taylor1` or `TaylorN`.
Note that for `TaylorN` the result of `a * b` is accumulated in `c[k]`.

The coefficients are given by

```math
c_k = \sum_{j=0}^k a_j b_{k-j}.
```

""" mul!


@inline function _check_homogeneous_product_order(c::HomogeneousPolynomial,
        a::HomogeneousPolynomial, b::HomogeneousPolynomial)
    order(c) == order(a) + order(b) ||
        throw(DimensionMismatch("result homogeneous degree must equal the sum of input degrees"))
    return nothing
end

@inline function _muladd_scalar_unchecked!(c::HomogeneousPolynomial, scalar,
        a::HomogeneousPolynomial)
    _isthinzero(scalar) && return nothing
    @inbounds for i in eachindex(c)
        ai = a[i]
        _isthinzero(ai) && continue
        c[i] += scalar * ai
    end
    return nothing
end

@inline function _mul_unchecked!(c::HomogeneousPolynomial, a::HomogeneousPolynomial,
        b::HomogeneousPolynomial)
    (_isthinzero(b) || _isthinzero(a)) && return nothing
    degree_a = order(a)
    degree_b = order(b)
    degree_a == 0 && return _muladd_scalar_unchecked!(c, a[1], b)
    degree_b == 0 && return _muladd_scalar_unchecked!(c, b[1], a)

    sp = c.space
    order_a = degree_a+1
    order_b = degree_b+1
    @inbounds num_coeffs_a = sp.size_table[order_a]
    @inbounds num_coeffs_b = sp.size_table[order_b]
    input_positions = _product_table(sp, degree_a, degree_b).input_positions
    pair = 1
    @inbounds for na in 1:num_coeffs_a
        ca = a[na]
        if _isthinzero(ca)
            pair += num_coeffs_b
            continue
        end
        @inbounds for nb in 1:num_coeffs_b
            cb = b[nb]
            if !_isthinzero(cb)
                pos = input_positions[pair]
                c[pos] += ca * cb
            end
            pair += 1
        end
    end
    return nothing
end

@inline function _mul_output_major_unchecked!(c::HomogeneousPolynomial,
        a::HomogeneousPolynomial, b::HomogeneousPolynomial)
    (_isthinzero(b) || _isthinzero(a)) && return nothing
    degree_a = order(a)
    degree_b = order(b)
    degree_a == 0 && return _muladd_scalar_unchecked!(c, a[1], b)
    degree_b == 0 && return _muladd_scalar_unchecked!(c, b[1], a)

    table = _init_output_major_product_table!(c.space, degree_a, degree_b)
    offsets = table.output_offsets
    output_pairs = table.output_pairs
    num_right = table.num_right
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    @inbounds for pos in 1:length(offsets)-1
        acc = c_coeffs[pos]
        for csr_pos in offsets[pos]:(offsets[pos+1]-1)
            pair = Int(output_pairs[csr_pos]) - 1
            na = pair ÷ num_right + 1
            nb = pair - (na-1) * num_right + 1
            acc += a_coeffs[na] * b_coeffs[nb]
        end
        c_coeffs[pos] = acc
    end
    return nothing
end

@inline function _mul_scalar_unchecked!(c::HomogeneousPolynomial,
        scalar::NumberNotSeries, a::HomogeneousPolynomial,
        b::HomogeneousPolynomial)
    (_isthinzero(scalar) || _isthinzero(b) || _isthinzero(a)) && return nothing
    degree_a = order(a)
    degree_b = order(b)
    degree_a == 0 && return _muladd_scalar_unchecked!(c, scalar * a[1], b)
    degree_b == 0 && return _muladd_scalar_unchecked!(c, scalar * b[1], a)

    sp = c.space
    order_a = degree_a+1
    order_b = degree_b+1
    @inbounds num_coeffs_a = sp.size_table[order_a]
    @inbounds num_coeffs_b = sp.size_table[order_b]
    input_positions = _product_table(sp, degree_a, degree_b).input_positions
    pair = 1
    @inbounds for na in 1:num_coeffs_a
        ca = a[na]
        if _isthinzero(ca)
            pair += num_coeffs_b
            continue
        end
        sca = scalar * ca
        @inbounds for nb in 1:num_coeffs_b
            cb = b[nb]
            if !_isthinzero(cb)
                pos = input_positions[pair]
                c[pos] += sca * cb
            end
            pair += 1
        end
    end
    return nothing
end

@inline function _muladd_unchecked!(c::TaylorN{T}, a::TaylorN{T},
        b::TaylorN{T}, k::Int) where {T<:Number}
    @inbounds _mul_output_major_unchecked!(c[k], a[0], b[k])
    @inbounds for i = 1:k
        _mul_output_major_unchecked!(c[k], a[i], b[k-i])
    end
    return nothing
end

@inline function _mul_scalar_unchecked!(c::TaylorN{T}, scalar::NumberNotSeries,
        a::TaylorN{T}, b::TaylorN{T}, k::Int) where {T<:Number}
    @inbounds _mul_scalar_unchecked!(c[k], scalar, a[0], b[k])
    @inbounds for i = 1:k
        _mul_scalar_unchecked!(c[k], scalar, a[i], b[k-i])
    end
    return nothing
end


"""
    mul!(c, a, b) --> nothing

Accumulates in `c` the result of `a*b` with minimum allocation. Arguments
c, a and b are `HomogeneousPolynomial`.

"""
@inline function mul!(c::HomogeneousPolynomial, a::HomogeneousPolynomial,
        b::HomogeneousPolynomial)
    _check_same_space(c, a, b)
    _check_homogeneous_product_order(c, a, b)
    _mul_output_major_unchecked!(c, a, b)
    return nothing
end


"""
    mul_scalar!(c, scalar, a, b) --> nothing

Accumulates in `c` the result of `scalar*a*b` with minimum allocation. Arguments
c, a and b are `HomogeneousPolynomial`; `scalar` is a NumberNotSeries.

"""
@inline function mul_scalar!(c::HomogeneousPolynomial, scalar::NumberNotSeries, a::HomogeneousPolynomial,
        b::HomogeneousPolynomial)
    _check_same_space(c, a, b)
    _check_homogeneous_product_order(c, a, b)
    _mul_scalar_unchecked!(c, scalar, a, b)
    return nothing
end
