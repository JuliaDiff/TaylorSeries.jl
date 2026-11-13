# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

# Hyperbolic functions: sinh, cosh, sinhcosh, tanh
#
# `sinhcosh!` for `Taylor1{Taylor1}` uses `_sincos_nested_recursion!`, defined in
# trigonometric.jl


# Out-of-place functions
for T in (:Taylor1, :TaylorN)
    @eval begin
        function sinhcosh(a::$T)
            order = TS.order(a)
            aux = sinh(constant_term(a))
            # aa = one(aux) * a
            aa = convert($T{typeof(aux)}, a)
            s = _constant_series_like(a, aux, order)
            c = _constant_series_like(a, cosh(constant_term(a)), order)
            for k in eachindex(a)
                sinhcosh!(s, c, aa, k)
            end
            return s, c
        end
        sinh(a::$T) = sinhcosh(a)[1]
        cosh(a::$T) = sinhcosh(a)[2]

        function tanh(a::$T)
            order = TS.order(a)
            aux = tanh( constant_term(a) )
            # aa = one(aux) * a
            aa = convert($T{typeof(aux)}, a)
            c = _constant_series_like(a, aux, order)
            c2 = _constant_series_like(a, aux^2, order)
            for k in eachindex(a)
                tanh!(c, aa, c2, k)
            end
            return c
        end
    end
end

# In-place functions (recursions) for Taylor1
function sinhcosh!(s::Taylor1{T}, c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    s_coeffs = s.coeffs
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    if k == 0
        @inbounds a0 = a_coeffs[1]
        @inbounds s_coeffs[1] = sinh(a0)
        @inbounds c_coeffs[1] = cosh(a0)
        return nothing
    end
    kk = k+1
    @inbounds s_acc = zero(s_coeffs[kk])
    @inbounds c_acc = zero(c_coeffs[kk])
    @inbounds for i = 1:k
        x = i * a_coeffs[i+1]
        s_acc += x * c_coeffs[k-i+1]
        c_acc += x * s_coeffs[k-i+1]
    end
    @inbounds s_coeffs[kk] = s_acc / k
    @inbounds c_coeffs[kk] = c_acc / k
    return nothing
end
function sinhcosh!(s::Taylor1{T}, c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:Number}
    if k == 0
        for j in eachindex(a[k])
            sinhcosh!(s[k], c[k], a[k], j)
        end
        return nothing
    end
    zero!(s[k])
    zero!(c[k])
    x = zero(a[k])
    @inbounds for i = 1:k
        for j in eachindex(a[k])
            mul_scalar!(x, i, a[i], c[k-i], j)
            add!(s[k], s[k], x, j)
            mul_scalar!(x, i, a[i], s[k-i], j)
            add!(c[k], c[k], x, j)
        end
    end
    for j in eachindex(a[k])
        div!(s[k], s[k], k, j)
        div!(c[k], c[k], k, j)
    end
    return nothing
end

function sinhcosh!(s::Taylor1{Taylor1{T}}, c::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}},
        k::Int) where {T<:NumberNotSeries}
    if k == 0
        for j in eachindex(a[k])
            sinhcosh!(s[k], c[k], a[k], j)
        end
        return nothing
    end
    # x = i * a[i]:  s[k] += x * c[k-i],  c[k] += x * s[k-i]
    _sincos_nested_recursion!(s, c, a, k, identity, identity)
    return nothing
end

function tanh!(c::Taylor1{T}, a::Taylor1{T}, c2::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    c2_coeffs = c2.coeffs
    if k == 0
        @inbounds aux = tanh(a_coeffs[1])
        @inbounds c_coeffs[1] = aux
        @inbounds c2_coeffs[1] = aux^2
        return nothing
    end
    kk = k+1
    @inbounds acc = zero(c_coeffs[kk])
    @inbounds for i = 1:k
        acc += (k-i+1) * a_coeffs[k-i+2] * c2_coeffs[i]
    end
    @inbounds c_coeffs[kk] = a_coeffs[kk] - acc/k
    @inbounds sqr!(c2, c, zero(c_coeffs[1]), k)
    return nothing
end
function tanh!(c::Taylor1{T}, a::Taylor1{T}, c2::Taylor1{T}, k::Int) where
        {T<:Number}
    if k == 0
        for j in eachindex(a[k])
            tanh!(c[0], a[0], c2[0], j)
        end
        return nothing
    end
    zero!(c, k)
    aux = zero(c[k])
    @inbounds for i = 0:k-1
        # c[k] += (k-i) * a[k-i] * c2[i]
        for j in eachindex(a[k])
            mul_scalar!(aux, k-i, a[k-i], c2[i], j)
            add!(c[k], c[k], aux, j)
        end
    end
    # @inbounds c[k] = a[k] - c[k]/k
    for j in eachindex(a[k])
        div!(c[k], c[k], k, j)     # c[k] <- c[k]/k
        subst!(c[k], a[k], c[k], j) # c[k] <- a[k] - c[k]
    end
    zero!(aux)
    sqr!(c2, c, aux, k)
    return nothing
end


# In-place functions (recursions) for TaylorN
function sinhcosh!(s::TaylorN{T}, c::TaylorN{T}, a::TaylorN{T}, k::Int) where
        {T<:Number}
    if k == 0
        @inbounds s[0] = sinh( constant_term(a) )
        @inbounds c[0] = cosh( constant_term(a) )
        return nothing
    end
    # Check spaces
    _check_same_space(s, c, a)
    # The recursion formula
    zero!(s, k)
    zero!(c, k)
    @inbounds for i = 1:k
        _mul_scalar_unchecked!(s[k], i, a[i], c[k-i])
        _mul_scalar_unchecked!(c[k], i, a[i], s[k-i])
    end
    @inbounds div!(s[k], s[k], k)
    @inbounds div!(c[k], c[k], k)
    return nothing
end

function tanh!(c::TaylorN{T}, a::TaylorN{T}, c2::TaylorN{T}, k::Int) where
        {T<:Number}
    if k == 0
        @inbounds aux = tanh( constant_term(a) )
        @inbounds c[0] = aux
        @inbounds c2[0] = aux^2
        return nothing
    end
    # Check spaces
    _check_same_space(c, a, c2)
    # The recursion formula
    zero!(c, k)
    @inbounds for i = 0:k-1
        _mul_scalar_unchecked!(c[k], k-i, a[k-i], c2[i])
    end
    @inbounds for l in eachindex(c[k])
        c[k][l] = a[k][l] - c[k][l]/k
    end
    sqr!(c2, c, zero(c[0][1]), k)
    return nothing
end


# In-place functions (recursions) for Taylor1{TaylorN}
function sinhcosh!(s::Taylor1{TaylorN{T}}, c::Taylor1{TaylorN{T}},
        a::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
    if k == 0
        @inbounds for ordQ in eachindex(a[0])
            sinhcosh!(s[0], c[0], a[0], ordQ)
        end
        return nothing
    end
    # Check spaces
    _check_same_space(s, c, a)
    # The recursion formula
    zero!(s[k])
    zero!(c[k])
    @inbounds for i = 1:k
        for ordQ in eachindex(a[0])
            # s[k] += i * a[i] * c[k-i];  c[k] += i * a[i] * s[k-i], in one pass
            _mul_scalar2_unchecked!(s[k], c[k], i, a[i], c[k-i], s[k-i], ordQ)
        end
    end
    div!(s, s, k, k)
    div!(c, c, k, k)
    return nothing
end

function tanh!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        res2::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
    if k == 0
        @inbounds res[0] = tanh( a[0] )
        # zero!(res2, res, 0)
        sqr!(res2, res, 0)
        return nothing
    end
    # Check spaces
    _check_same_space(res, a, res2)
    # The recursion formula
    zero!(res[k])
    for i = 0:k-1
        @inbounds for ordQ in eachindex(a[0])
            _mul_scalar_unchecked!(res[k], k-i, res2[i], a[k-i], ordQ)
        end
    end
    aux = res2[k]
    @inbounds for ordQ in eachindex(a[0])
        # zero!(tmp, res[k], ordQ)
        # aux[ordQ] <- res[k][ordQ] / k, in place
        aux_hp = aux[ordQ].coeffs
        res_hp = res[k][ordQ].coeffs
        for j in eachindex(aux_hp)
            aux_hp[j] = res_hp[j] / k
        end
        subst!(res[k], a[k], aux, ordQ)
    end
    zero!(res2[k])
    sqr!(res2, res, k)
    return nothing
end


# Documentation of the recursion relations
@doc doc"""
    sinhcosh!(s, c, a, k)

Update the `k-th` expansion coefficients `s[k+1]` and `c[k+1]`
of `s = sinh(a)` and `c = cosh(a)` simultaneously, for `s`, `c` and `a`
either `Taylor1` or `TaylorN`.

The coefficients are given by

```math
\begin{aligned}
s_k = \frac{1}{k} \sum_{j=0}^{k-1} (k-j) a_{k-j} c_j, \\
c_k = \frac{1}{k} \sum_{j=0}^{k-1} (k-j) a_{k-j} s_j.
\end{aligned}
```

""" sinhcosh!


@doc doc"""
    tanh!(c, a, p, k)

Update the `k-th` expansion coefficients `c[k+1]` of `c = tanh(a)`,
for `c` and `a` either `Taylor1` or `TaylorN`; `p = a^2` and
is passed as an argument for efficiency.

```math
\begin{equation*}
c_k = a_k - \frac{1}{k} \sum_{j=0}^{k-1} (k-j) a_{k-j} p_j.
\end{equation*}
```

""" tanh!