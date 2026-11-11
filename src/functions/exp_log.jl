# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

# Exponential and logarithmic functions: exp, expm1, log, log1p


# Out-of-place functions
for T in (:Taylor1, :TaylorN)
    @eval begin
        function exp(a::$T)
            aux = exp(constant_term(a))
            aa = convert($T{typeof(aux)}, a)
            c = _constant_series_like(a, aux, order(a))
            for k in eachindex(a)
                exp!(c, aa, k)
            end
            return c
        end

        function expm1(a::$T)
            aux = expm1(constant_term(a))
            aa = convert($T{typeof(aux)}, a)
            c = _constant_series_like(a, aux, order(a))
            for k in eachindex(a)
                expm1!(c, aa, k)
            end
            return c
        end

        function log(a::$T)
            iszero(constant_term(a)) && throw(DomainError(a,
                """The 0-th order coefficient must be non-zero in order to expand `log` around 0."""))
            aux = log(constant_term(a))
            aa = convert($T{typeof(aux)}, a)
            c = _constant_series_like(a, aux, order(a))
            for k in eachindex(a)
                log!(c, aa, k)
            end
            return c
        end

        function log1p(a::$T)
            # constant_term(a) < -one(constant_term(a)) && throw(DomainError(a,
            #         """The 0-th order coefficient must be larger than -1 in order to expand `log1`."""))
            aux = log1p(constant_term(a))
            # aa = one(aux) * a
            aa = convert($T{typeof(aux)}, a)
            c = _constant_series_like(a, aux, order(a))
            for k in eachindex(a)
                log1p!(c, aa, k)
            end
            return c
        end
    end
end

function log1p(a::Taylor1{TaylorN{T}}) where {T<:NumberNotSeries}
    aux = log1p(constant_term(a))
    aa = convert(Taylor1{typeof(aux)}, a)
    c = _constant_series_like(a, aux, order(a))
    denom = _constant_series_like(aux, zero(aux[0][1]), order(aux))
    _log1p_denominator!(denom, aa.coeffs[1])
    for k in eachindex(a)
        log1p!(c, aa, k, denom)
    end
    return c
end


# In-place functions (recursions) for Taylor1
function exp!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    if k == 0
        @inbounds c_coeffs[1] = exp(a_coeffs[1])
        return nothing
    end
    kk = k+1
    @inbounds acc = zero(c_coeffs[kk])
    @inbounds for i = 1:k
        acc += (k-i+1) * a_coeffs[k-i+2] * c_coeffs[i]
    end
    @inbounds c_coeffs[kk] = acc / k
    return nothing
end
function exp!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:Number}
    if k == 0
        a0 = constant_term(a)
        for j in eachindex(a[k])
            exp!(c[k], a0, j)
        end
        return nothing
    end
    zero!(c[k])
    aux = zero(c[k])
    @inbounds for i = 0:k-1
        # c[k] += (k-i) * a[k-i] * c[i]
        for j in eachindex(c[k])
            mul_scalar!(aux, k-i, a[k-i], c[i], j)
            add!(c[k], c[k], aux, j)
            zero!(aux, j)
        end
    end
    # div!(c, c, k, k)
    for i in eachindex(c[k])
        identity!(aux, c[k], i)
    end
    for i in eachindex(c[k])
        div!(c[k], aux, k, i)
    end
    return nothing
end

function expm1!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    if k == 0
        @inbounds c.coeffs[1] = expm1(a.coeffs[1])
        return nothing
    end
    c_coeffs = c.coeffs
    @inbounds cc = c_coeffs[1]
    @inbounds c_coeffs[1] += one(cc)
    exp!(c, a, k)
    @inbounds c_coeffs[1] = cc
    return nothing
end
function expm1!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:Number}
    if k == 0
        for j in eachindex(a[k])
            expm1!(c[k], constant_term(a), j)
        end
        return nothing
    end
    cc = c[0]
    c[0] += one(c[0])
    exp!(c, a, k)
    c[0] = cc
    return nothing
end

function log!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    @inbounds a0 = a_coeffs[1]
    if k == 0
        @inbounds c_coeffs[1] = log(a0)
        return nothing
    elseif k == 1
        @inbounds c_coeffs[2] = a_coeffs[2] / a0
        return nothing
    end
    kk = k+1
    @inbounds acc = zero(c_coeffs[kk])
    @inbounds for i = 1:k-1
        acc += (k-i) * a_coeffs[i+1] * c_coeffs[k-i+1]
    end
    @inbounds c_coeffs[kk] = (a_coeffs[kk] - acc/k) / a0
    return nothing
end
function log!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:Number}
    if k == 0
        aux = constant_term(a)
        # @inbounds c[0] = log(constant_term(a))
        for j in eachindex(a[k])
            log!(c[0], aux, j)
        end
        return nothing
    elseif k == 1
        aux = constant_term(a)
        # @inbounds c[1] = a[1] / constant_term(a)
        for j in eachindex(a[1])
            div!(c[1], a[1], aux, j)
        end
        return nothing
    end
    zero!(c[k])
    aux = zero(c[k])
    @inbounds for i = 1:k-1
        # c[k] += (k-i) * a[i] * c[k-i]
        for j in eachindex(a[k])
            mul_scalar!(aux, k-i, a[i], c[k-i], j)
            add!(c[k], c[k], aux, j)
        end
    end
    # @inbounds c[k] = (a[k] - c[k]/k) / constant_term(a)
    for j in eachindex(a[k])
        zero!(aux, j)
        div!(aux, c[k], k, j)
        subst!(aux, a[k], aux, j)
    end
    for j in eachindex(a[k])
        div!(c[k], aux, constant_term(a), j)
    end
    return nothing
end

function log1p!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    @inbounds a0 = a_coeffs[1]
    a0p1 = a0+one(a0)
    if k == 0
        @inbounds c_coeffs[1] = log1p(a0)
        return nothing
    elseif k == 1
        @inbounds c_coeffs[2] = a_coeffs[2] / a0p1
        return nothing
    end
    kk = k+1
    @inbounds acc = zero(c_coeffs[kk])
    @inbounds for i = 1:k-1
        acc += (k-i) * a_coeffs[i+1] * c_coeffs[k-i+1]
    end
    @inbounds c_coeffs[kk] = (a_coeffs[kk] - acc/k) / a0p1
    return nothing
end
function log1p!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:Number}
    a0 = constant_term(a)
    a0p1 = a0+one(a0)
    if k == 0
        for j in eachindex(a[0])
            @inbounds log1p!(c[0], a0, j)
        end
        return nothing
    elseif k == 1
        # @inbounds c[1] = a[1] / a0p1
        for j in eachindex(a[1])
            div!(c[1], a[1], a0p1, j)
        end
        return nothing
    end
    zero!(c[k])
    a0 = zero(c[k])
    @inbounds for i = 1:k-1
        # c[k] += (k-i) * a[i] * c[k-i]
        for j in eachindex(a[k])
            mul_scalar!(a0, k-i, a[i], c[k-i], j)
            add!(c[k], c[k], a0, j)
        end
    end
    # @inbounds c[k] = (a[k] - c[k]/k) / a0p1
    for j in eachindex(a[k])
        zero!(a0, j)
        div!(a0, c[k], k, j)
        subst!(a0, a[k], a0, j)
    end
    for j in eachindex(c[k])
        div!(c[k], a0, a0p1, j)
    end
    return nothing
end


# In-place functions (recursions) for TaylorN
function exp!(c::TaylorN{T}, a::TaylorN{T}, k::Int) where {T<:Number}
    if k == 0
        @inbounds c[0] = exp(constant_term(a))
        return nothing
    end
    zero!(c, k)
    @inbounds for i = 0:k-1
        mul_scalar!(c[k], k-i, a[k-i], c[i])
    end
    @inbounds div!(c[k], c[k], k)
    return nothing
end

function expm1!(c::TaylorN{T}, a::TaylorN{T}, k::Int) where {T<:Number}
    if k == 0
        @inbounds c[0] = expm1(constant_term(a))
        return nothing
    end
    zero!(c, k)
    # c[k] += k * a[k] * (c[0] + 1): the factor c[0] + 1 is a constant (degree 0), so this is
    # the degree-0 branch of `mul_scalar!`, without building that polynomial
    @inbounds c00 = c[0][1]
    @inbounds _muladd_scalar_unchecked!(c[k], k * (c00 + one(c00)), a[k])
    @inbounds for i = 1:k-1
        mul_scalar!(c[k], k-i, a[k-i], c[i])
    end
    @inbounds div!(c[k], c[k], k)
    return nothing
end

function log!(c::TaylorN{T}, a::TaylorN{T}, k::Int) where {T<:Number}
    if k == 0
        @inbounds c[0] = log(constant_term(a))
        return nothing
    elseif k == 1
        a0 = constant_term(a)
        @inbounds c_1 = c[1].coeffs
        @inbounds a_1 = a[1].coeffs
        @inbounds for j in eachindex(c_1)
            c_1[j] = a_1[j] / a0                    # c[1] <- a[1] / a0
        end
        return nothing
    end
    zero!(c, k)
    @inbounds for i = 1:k-1
        mul_scalar!(c[k], k-i, a[i], c[k-i])
    end
    # c[k] <- (a[k] - c[k]/k) / a0
    a0 = constant_term(a)
    @inbounds c_k = c[k].coeffs
    @inbounds a_k = a[k].coeffs
    @inbounds for j in eachindex(c_k)
        c_k[j] = (a_k[j] - c_k[j]/k) / a0
    end
    return nothing
end

function log1p!(c::TaylorN{T}, a::TaylorN{T}, k::Int) where {T<:Number}
    if k == 0
        a0 = constant_term(a)
        @inbounds c[0] = log1p(a0)
        return nothing
    elseif k == 1
        a0 = constant_term(a)
        a0p1 = a0+one(a0)
        @inbounds c_1 = c[1].coeffs
        @inbounds a_1 = a[1].coeffs
        @inbounds for j in eachindex(c_1)
            c_1[j] = a_1[j] / a0p1                  # c[1] <- a[1] / a0p1
        end
        return nothing
    end
    a0 = constant_term(a)
    a0p1 = a0+one(a0)
    zero!(c, k)
    @inbounds for i = 1:k-1
        mul_scalar!(c[k], k-i, a[i], c[k-i])
    end
    # c[k] <- (a[k] - c[k]/k) / a0p1
    @inbounds c_k = c[k].coeffs
    @inbounds a_k = a[k].coeffs
    @inbounds for j in eachindex(c_k)
        c_k[j] = (a_k[j] - c_k[j]/k) / a0p1
    end
    return nothing
end


# In-place functions (recursions) for Taylor1{TaylorN}
function exp!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        k::Int) where {T<:NumberNotSeries}
    if k == 0
        @inbounds for ordQ in eachindex(a[0])
            # zero!(res[0], a[0], ordQ)
            exp!(res[0], a[0], ordQ)
        end
        return nothing
    end
    # The recursion formula
    zero!(res[k])
    for i = 0:k-1
        @inbounds for ordQ in eachindex(a[0])
            mul_scalar!(res[k], k-i, res[i], a[k-i], ordQ)
        end
    end
    div!(res, res, k, k)
    return nothing
end

function expm1!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        k::Int) where {T<:NumberNotSeries}
    if k == 0
        @inbounds for ordQ in eachindex(a[0])
            # zero!(res[0], a[0], ordQ)
            expm1!(res[0], a[0], ordQ)
        end
        return nothing
    end
    # The recursion formula
    tmp = _constant_series_like(a[k], zero(a[k][0][1]), order(a[0]))
    zero!(res[k])
    # i=0 term of sum
    # tmp[ordQ] is scaled in place
    @inbounds for ordQ in eachindex(a[0])
        one!(tmp, a[0], ordQ)
        add!(tmp, res[0], tmp, ordQ)
        tmp_hp = tmp[ordQ].coeffs
        for j in eachindex(tmp_hp)
            tmp_hp[j] = k * tmp_hp[j]
        end
        # zero!(res[k], a[0], ordQ)
        mul!(res[k], a[k], tmp, ordQ)
    end
    for i = 1:k-1
        @inbounds for ordQ in eachindex(a[0])
            # tmp[ordQ] <- (k-i) * res[i][ordQ], in place
            tmp_hp = tmp[ordQ].coeffs
            res_hp = res[i][ordQ].coeffs
            for j in eachindex(tmp_hp)
                tmp_hp[j] = (k-i) * res_hp[j]
            end
            mul!(res[k], tmp, a[k-i], ordQ)
        end
    end
    div!(res, res, k, k)
    return nothing
end

function log!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        k::Int) where {T<:NumberNotSeries}
    res_coeffs = res.coeffs
    a_coeffs = a.coeffs
    res_k = res_coeffs[k+1]
    a0 = a_coeffs[1]
    if k == 0
        @inbounds for ordQ in eachindex(a0)
            # zero!(res[k], a[0], ordQ)
            log!(res_k, a0, ordQ)
        end
        return nothing
    elseif k == 1
        a1 = a_coeffs[2]
        @inbounds for ordQ in eachindex(a0)
            zero!(res_k[ordQ])
            div!(res_k, a1, a0, ordQ)
        end
        return nothing
    end
    # The recursion formula
    zero!(res_k)
    for i = 1:k-1
        res_ki = res_coeffs[k-i+1]
        a_i = a_coeffs[i+1]
        @inbounds for ordQ in eachindex(a0)
            mul_scalar!(res_k, k-i, res_ki, a_i, ordQ)
        end
    end
    div!(res, res, k, k)
    a_k = a_coeffs[k+1]
    @inbounds for ordQ in eachindex(a0)
        subst!(res_k, a_k, res_k, ordQ)
        div!(res_k, a0, ordQ)
    end
    return nothing
end

function log1p!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        k::Int) where {T<:NumberNotSeries}
    res_coeffs = res.coeffs
    a_coeffs = a.coeffs
    res_k = res_coeffs[k+1]
    a0 = a_coeffs[1]
    if k == 0
        @inbounds for ordQ in eachindex(a0)
            # zero!(res[k], a[0], ordQ)
            log1p!(res_k, a0, ordQ)
        end
        return nothing
    end
    tmp1 = _constant_series_like(a_coeffs[k+1], zero(a_coeffs[k+1][0][1]), order(a0))
    _log1p_denominator!(tmp1, a0)
    log1p!(res, a, k, tmp1)
    return nothing
end


function _log1p_denominator!(tmp::TaylorN{T}, a0::TaylorN{T}) where {T<:Number}
    zero!(tmp)
    @inbounds for ordQ in eachindex(a0)
        one!(tmp, a0, ordQ)
        add!(tmp, tmp, a0, ordQ)
    end
    return tmp
end

function log1p!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        k::Int, tmp1::TaylorN{T}) where {T<:NumberNotSeries}
    res_coeffs = res.coeffs
    a_coeffs = a.coeffs
    res_k = res_coeffs[k+1]
    a0 = a_coeffs[1]
    if k == 0
        @inbounds for ordQ in eachindex(a0)
            log1p!(res_k, a0, ordQ)
        end
        return nothing
    end
    zero!(res_k)
    if k == 1
        a1 = a_coeffs[2]
        @inbounds for ordQ in eachindex(a0)
            div!(res_k, a1, tmp1, ordQ)
        end
        return nothing
    end
    # The recursion formula
    for i = 1:k-1
        res_ki = res_coeffs[k-i+1]
        a_i = a_coeffs[i+1]
        @inbounds for ordQ in eachindex(a0)
            mul_scalar!(res_k, k-i, res_ki, a_i, ordQ)
        end
    end
    div!(res, res, k, k)
    a_k = a_coeffs[k+1]
    @inbounds for ordQ in eachindex(a0)
        subst!(res_k, a_k, res_k, ordQ)
        div!(res_k, tmp1, ordQ)
    end
    return nothing
end


# Documentation of the recursion relations
@doc doc"""
    exp!(c, a, k) --> nothing

Update the `k-th` expansion coefficient `c[k+1]` of `c = exp(a)`
for both `c` and `a` either `Taylor1` or `TaylorN`.

The coefficients are given by

```math
\begin{equation*}
c_k = \frac{1}{k} \sum_{j=0}^{k-1} (k-j) a_{k-j} c_j.
\end{equation*}
```

""" exp!


@doc doc"""
    log!(c, a, k) --> nothing

Update the `k-th` expansion coefficient `c[k+1]` of `c = log(a)`
for both `c` and `a` either `Taylor1` or `TaylorN`.

The coefficients are given by

```math
\begin{equation*}
c_k = \frac{1}{a_0} \big(a_k - \frac{1}{k} \sum_{j=0}^{k-1} j a_{k-j} c_j \big).
\end{equation*}
```

""" log!