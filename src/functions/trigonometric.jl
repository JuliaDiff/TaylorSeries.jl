# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

# Trigonometric functions: sin, cos, sincos, sinpi, cospi, sincospi, cis, cispi, tan

# Construct a complex series depending on the coefficients type, use for `cis` and `cispi`
_complex_series(c::T, s::T) where {T<:Real} = complex(c, s)
_complex_series(c::T, s::T) where {T<:Complex} = c + im * s
function _complex_series(c::Taylor1, s::Taylor1)
    c_coeffs = c.coeffs
    s_coeffs = s.coeffs
    z1 = _complex_series(c_coeffs[1], s_coeffs[1])
    v = FixedSizeVectorDefault{typeof(z1)}(undef, length(c_coeffs))
    v[1] = z1
    @inbounds for i in 2:length(c_coeffs)
        v[i] = _complex_series(c_coeffs[i], s_coeffs[i])
    end
    return Taylor1(v)
end
function _complex_series(c::HomogeneousPolynomial, s::HomogeneousPolynomial)
    c_coeffs = c.coeffs
    s_coeffs = s.coeffs
    # one allocation (the coefficient array), then filled in place
    z = HomogeneousPolynomial(c.space, _complex_series(c_coeffs[1], s_coeffs[1]), order(c))
    z_coeffs = z.coeffs
    @inbounds for i in 2:length(c_coeffs)
        z_coeffs[i] = _complex_series(c_coeffs[i], s_coeffs[i])
    end
    return z
end
function _complex_series(c::TaylorN, s::TaylorN)
    c_hps = c.coeffs
    s_hps = s.coeffs
    z1 = _complex_series(c_hps[1], s_hps[1])
    v = FixedSizeVectorDefault{typeof(z1)}(undef, length(c_hps))
    v[1] = z1
    @inbounds for i in 2:length(c_hps)
        v[i] = _complex_series(c_hps[i], s_hps[i])
    end
    return TaylorN{numtype(z1)}(v, c.space) # trusted constructor: same space and orders as `c`
end


# Out-of-place functions
for T in (:Taylor1, :TaylorN)
    @eval begin
        function sincos(a::$T)
            order = TS.order(a)
            aux, auxc = sincos(constant_term(a))
            # aa = one(aux) * a
            aa = convert($T{typeof(aux)}, a)
            s = _constant_series_like(a, aux, order)
            c = _constant_series_like(a, auxc, order)
            for k in eachindex(a)
                sincos!(s, c, aa, k)
            end
            return s, c
        end
        sin(a::$T) = sincos(a)[1]
        cos(a::$T) = sincos(a)[2]

        function sincospi(a::$T)
            order = TS.order(a)
            aux, auxc = sincospi(constant_term(a))
            aa = one(aux) * a
            # aa = convert($T{typeof(aux)}, a)
            s = _constant_series_like(a, aux, order)
            c = _constant_series_like(a, auxc, order)
            for k in eachindex(a)
                sincospi!(s, c, aa, k)
            end
            return s, c
        end
        sinpi(a::$T) = sincospi(a)[1]
        cospi(a::$T) = sincospi(a)[2]

        function cis(a::$T)
            s, c = sincos(a)
            return _complex_series(c, s)
        end

        function cispi(a::$T)
            s, c = sincospi(a)
            return _complex_series(c, s)
        end

        function tan(a::$T)
            order = TS.order(a)
            aux = tan(constant_term(a))
            # aa = one(aux) * a
            aa = convert($T{typeof(aux)}, a)
            c = _constant_series_like(a, aux, order)
            c2 = _constant_series_like(a, aux^2, order)
            for k in eachindex(a)
                tan!(c, aa, c2, k)
            end
            return c
        end
    end
end


# In-place functions (recursions) for Taylor1
function sincos!(s::Taylor1{T}, c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    s_coeffs = s.coeffs
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    if k == 0
        @inbounds s_coeffs[1], c_coeffs[1] = sincos(a_coeffs[1])
        return nothing
    end
    kk = k+1
    @inbounds s_acc = zero(s_coeffs[kk])
    @inbounds c_acc = zero(c_coeffs[kk])
    @inbounds for i = 1:k
        x = i * a_coeffs[i+1]
        s_acc += x * c_coeffs[k-i+1]
        c_acc -= x * s_coeffs[k-i+1]
    end
    @inbounds s_coeffs[kk] = s_acc / k
    @inbounds c_coeffs[kk] = c_acc / k
    return nothing
end
function sincos!(s::Taylor1{T}, c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:Number}
    if k == 0
        for j in eachindex(a[k])
            @inbounds sincos!(s[0], c[0], constant_term(a), j)
        end
        return nothing
    end
    x = zero(c[k])
    zero!(s[k])
    zero!(c[k])
    @inbounds for i = 1:k
        # x = i * a[i]
        # s[k] += x * c[k-i]
        # c[k] -= x * s[k-i]
        for j in eachindex(a[k])
            mul_scalar!(x,  i, a[i], c[k-i], j)
            add!(s[k], s[k], x, j)
            mul_scalar!(x, -i, a[i], s[k-i], j)
            add!(c[k], c[k], x, j)
        end
    end
    for j in eachindex(a[k])
        div!(s[k], s[k], k, j)
        div!(c[k], c[k], k, j)
    end
    return nothing
end

# `Taylor1{Taylor1}`: the recursion (k ≥ 1) shared by `sincos!`, `sincospi!` and `sinhcosh!`
@inline function _sincos_nested_recursion!(s::Taylor1{Taylor1{T}}, c::Taylor1{Taylor1{T}},
        a::Taylor1{Taylor1{T}}, k::Int, fs::FS, fc::FC) where {T<:NumberNotSeries, FS, FC}
    @inbounds s_k = s[k]
    @inbounds c_k = c[k]
    zero!(s_k)
    zero!(c_k)
    s_kc = s_k.coeffs
    c_kc = c_k.coeffs
    @inbounds for i = 1:k
        sc_s = fs(i)
        sc_c = fc(i)
        a_ic = a[i].coeffs
        c_kic = c[k-i].coeffs
        s_kic = s[k-i].coeffs
        for j in eachindex(a[k])
            jj = j + 1
            s_acc = zero(T)
            c_acc = zero(T)
            for l = 1:jj
                a_il = a_ic[l]
                s_acc += a_il * c_kic[jj-l+1]
                c_acc += a_il * s_kic[jj-l+1]
            end
            s_kc[jj] = s_kc[jj] + convert(T, s_acc * sc_s)
            c_kc[jj] = c_kc[jj] + convert(T, c_acc * sc_c)
        end
    end
    @inbounds for jj in eachindex(s_kc)
        s_kc[jj] = s_kc[jj] / k
        c_kc[jj] = c_kc[jj] / k
    end
    return nothing
end

function sincos!(s::Taylor1{Taylor1{T}}, c::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}},
        k::Int) where {T<:NumberNotSeries}
    if k == 0
        for j in eachindex(a[k])
            @inbounds sincos!(s[0], c[0], constant_term(a), j)
        end
        return nothing
    end
    # x = i * a[i]:  s[k] += x * c[k-i],  c[k] -= x * s[k-i]
    _sincos_nested_recursion!(s, c, a, k, identity, -)
    return nothing
end


# The factor `pi` of the recursions of `sincospi!`; by default is `pi`, but
# for e.g. intervals returns a proper enclosure.
_pi_factor(::Type{T}) where {T<:Number} = pi
_pi_factor(::Type{<:AbstractSeries{T}}) where {T} = _pi_factor(T)


function sincospi!(s::Taylor1{T}, c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    s_coeffs = s.coeffs
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    if k == 0
        @inbounds s_coeffs[1], c_coeffs[1] = sincospi(a_coeffs[1])
        return nothing
    end
    kk = k+1
    pi_T = _pi_factor(T)
    @inbounds s_acc = zero(s_coeffs[kk])
    @inbounds c_acc = zero(c_coeffs[kk])
    @inbounds for i = 1:k
        x = (i * pi_T) * a_coeffs[i+1]
        s_acc += x * c_coeffs[k-i+1]
        c_acc -= x * s_coeffs[k-i+1]
    end
    @inbounds s_coeffs[kk] = s_acc / k
    @inbounds c_coeffs[kk] = c_acc / k
    return nothing
end
function sincospi!(s::Taylor1{T}, c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:Number}
    if k == 0
        a0 = constant_term(a)
        for j in eachindex(a[k])
            @inbounds sincospi!(s[0], c[0], a0, j)
        end
        return nothing
    end
    x = zero(c[k])
    zero!(s[k])
    zero!(c[k])
    pi_T = _pi_factor(T)
    @inbounds for i = 1:k
        # x = i * pi * a[i]
        # s[k] += x * c[k-i]
        # c[k] -= x * s[k-i]
        for j in eachindex(a[k])
            mul_scalar!(x,  i*pi_T, a[i], c[k-i], j)
            add!(s[k], s[k], x, j)
            mul_scalar!(x, -i*pi_T, a[i], s[k-i], j)
            add!(c[k], c[k], x, j)
        end
    end
    for j in eachindex(a[k])
        div!(s[k], s[k], k, j)
        div!(c[k], c[k], k, j)
    end
    return nothing
end
function sincospi!(s::Taylor1{Taylor1{T}}, c::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}},
        k::Int) where {T<:NumberNotSeries}
    if k == 0
        a0 = constant_term(a)
        for j in eachindex(a[k])
            @inbounds sincospi!(s[0], c[0], a0, j)
        end
        return nothing
    end
    # x = i * pi_T * a[i]:  s[k] += x * c[k-i],  c[k] -= x * s[k-i]
    pi_T = _pi_factor(T)
    _sincos_nested_recursion!(s, c, a, k, i -> i*pi_T, i -> -i*pi_T)
    return nothing
end


function tan!(c::Taylor1{T}, a::Taylor1{T}, c2::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    c2_coeffs = c2.coeffs
    if k == 0
        @inbounds aux = tan(a_coeffs[1])
        @inbounds c_coeffs[1] = aux
        @inbounds c2_coeffs[1] = aux^2
        return nothing
    end
    kk = k+1
    @inbounds acc = zero(c_coeffs[kk])
    @inbounds for i = 1:k
        acc += (k-i+1) * a_coeffs[k-i+2] * c2_coeffs[i]
    end
    @inbounds c_coeffs[kk] = acc / k + a_coeffs[kk]
    @inbounds sqr!(c2, c, zero(c_coeffs[1]), k)
    return nothing
end
function tan!(c::Taylor1{T}, a::Taylor1{T}, c2::Taylor1{T}, k::Int) where
        {T<:Number}
    if k == 0
        aux = constant_term(a)
        @inbounds for j in eachindex(a[k])
            tan!(c[0], aux, c2[0], j)
        end
        return nothing
    end
    zero!(c[k])
    aux = zero(c[k])
    @inbounds for i = 0:k-1
        # c[k] += (k-i) * a[k-i] * c2[i]
        for j in eachindex(a[k])
            mul_scalar!(aux, k-i, a[k-i], c2[i], j)
            add!(c[k], c[k], aux, j)
        end
    end
    for j in eachindex(a[k])
        div!(c[k], c[k], k, j) # c[k] <- c[k]/k
        add!(c[k], c[k], a[k], j) # c[k] <- c[k] + a[k]
        zero!(aux, j)
    end
    sqr!(c2, c, aux, k)
    return nothing
end


# In-place functions (recursions) for TaylorN
function sincos!(s::TaylorN{T}, c::TaylorN{T}, a::TaylorN{T}, k::Int) where
        {T<:Number}
    if k == 0
        a0 = constant_term(a)
        @inbounds s[0][1], c[0][1] = sincos( a0 )
        return nothing
    end
    zero!(s, k)
    zero!(c, k)
    @inbounds for i = 1:k
        mul_scalar!(s[k],  i, a[i], c[k-i])
        mul_scalar!(c[k], -i, a[i], s[k-i])
    end
    @inbounds div!(s[k], s[k], k)
    @inbounds div!(c[k], c[k], k)
    return nothing
end

function sincospi!(s::TaylorN{T}, c::TaylorN{T}, a::TaylorN{T}, k::Int) where
        {T<:Number}
    if k == 0
        a0 = constant_term(a)
        @inbounds s[0][1], c[0][1] = sincospi( a0 )
        return nothing
    end
    # Check spaces
    _check_same_space(s, c, a)
    zero!(s, k)
    zero!(c, k)
    @inbounds for i = 1:k
        _mul_scalar_unchecked!(s[k],  i, a[i], c[k-i], pi)
        _mul_scalar_unchecked!(c[k], -i, a[i], s[k-i], pi)
    end
    @inbounds div!(s[k], s[k], k)
    @inbounds div!(c[k], c[k], k)
    return nothing
end

function tan!(c::TaylorN{T}, a::TaylorN{T}, c2::TaylorN{T}, k::Int) where
        {T<:Number}
    a0 = constant_term(a)
    if k == 0
        @inbounds aux = tan(a0)
        @inbounds c[0][1] = aux
        @inbounds c2[0][1] = aux^2
        return nothing
    end
    zero!(c, k)
    @inbounds for i = 0:k-1
        mul_scalar!(c[k], k-i, a[k-i], c2[i])
    end
    # c[k] <- c[k]/k
    div!(c[k], c[k], k)
    # c[k] <- c[k] + a[k]
    add!(c, a, c, k)
    sqr!(c2, c, zero(a0), k)
    return nothing
end


# In-place functions (recursions) for Taylor1{TaylorN}
function sincos!(s::Taylor1{TaylorN{T}}, c::Taylor1{TaylorN{T}},
        a::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
    if k == 0
        @inbounds for ordQ in eachindex(a[0])
            sincos!(s[0], c[0], a[0], ordQ)
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
            # s[k] += i * a[i] * c[k-i];  c[k] += i * a[i] * s[k-i]
            _mul_scalar2_unchecked!(s[k], c[k], i, a[i], c[k-i], s[k-i], ordQ)
        end
    end
    div!(s, s, k, k)
    subst!(c, c, k)
    div!(c, c, k, k)
    return nothing
end

function sincospi!(s::Taylor1{TaylorN{T}}, c::Taylor1{TaylorN{T}},
        a::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
    if k == 0
        @inbounds for ordQ in eachindex(a[0])
            sincospi!(s[0], c[0], a[0], ordQ)
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
            # s[k] += i * (pi*a[i]) * c[k-i];  c[k] += i * (pi*a[i]) * s[k-i]
            _mul_scalar2_unchecked!(s[k], c[k], i, a[i], c[k-i], s[k-i], ordQ, pi)
        end
    end
    div!(s, s, k, k)
    subst!(c, c, k)
    div!(c, c, k, k)
    return nothing
end

function tan!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        res2::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
    if k == 0
        @inbounds res[0] = tan( a[0] )
        # zero!(res2, res, 0)
        sqr!(res2, res, 0)
        return nothing
    end
    # The recursion formula
    zero!(res[k])
    for i = 0:k-1
        @inbounds for ordQ in eachindex(a[0])
            mul_scalar!(res[k], k-i, res2[i], a[k-i], ordQ)
        end
    end
    @inbounds for ordQ in eachindex(a[0])
        div!(res[k][ordQ], res[k][ordQ], k)
        add!(res[k], a[k], res[k], ordQ)
    end
    sqr!(res2, res, k)
    return nothing
end


# Documentation of the recursion relations
@doc doc"""
    sincos!(s, c, a, k) --> nothing

Update the `k-th` expansion coefficients `s[k+1]` and `c[k+1]`
of `s = sin(a)` and `c = cos(a)` simultaneously, for `s`, `c` and `a`
either `Taylor1` or `TaylorN`.

The coefficients are given by

```math
\begin{aligned}
s_k &=  \frac{1}{k}\sum_{j=0}^{k-1} (k-j) a_{k-j} c_j ,\\
c_k &= -\frac{1}{k}\sum_{j=0}^{k-1} (k-j) a_{k-j} s_j.
\end{aligned}
```

""" sincos!


@doc doc"""
    tan!(c, a, p, k::Int) --> nothing

Update the `k-th` expansion coefficients `c[k+1]` of `c = tan(a)`,
for `c` and `a` either `Taylor1` or `TaylorN`; `p = c^2` and
is passed as an argument for efficiency.

The coefficients are given by

```math
\begin{equation*}
c_k = a_k + \frac{1}{k} \sum_{j=0}^{k-1} (k-j) a_{k-j} p_j.
\end{equation*}
```

""" tan!