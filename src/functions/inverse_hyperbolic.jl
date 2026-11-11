# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

# Inverse hyperbolic functions: asinh, acosh, atanh


# Out-of-place functions
for T in (:Taylor1, :TaylorN)
    @eval begin
        function asinh(a::$T)
            order = TS.order(a)
            a0 = constant_term(a)
            aux = asinh(a0)
            # aa = one(aux) * a
            aa = convert($T{typeof(aux)}, a)
            c = _constant_series_like(a, aux, order)
            r = _constant_series_like(a, sqrt(a0^2 + 1), order)
            iszero(constant_term(r)) && throw(DomainError(a,
                """Series expansion of asinh(x) diverges at x = ±im."""))
            for k in eachindex(a)
                asinh!(c, aa, r, k)
            end
            return c
        end

        function acosh(a::$T)
            a0 = constant_term(a)
            a0^2 == one(a0) && throw(DomainError(a,
                """Series expansion of acosh(x) diverges at x = ±1."""))
            order = TS.order(a)
            aux = acosh(a0)
            # aa = one(aux) * a
            aa = convert($T{typeof(aux)}, a)
            c = _constant_series_like(a, aux, order)
            r = _constant_series_like(a, sqrt(a0^2 - 1), order)
            for k in eachindex(a)
                acosh!(c, aa, r, k)
            end
            return c
        end

        function atanh(a::$T)
            order = TS.order(a)
            a0 = constant_term(a)
            aux = atanh(a0)
            # aa = one(aux) * a
            aa = convert($T{typeof(aux)}, a)
            c = _constant_series_like(a, aux, order)
            r = _constant_series_like(a, one(aux) - a0^2, order)
            iszero(constant_term(r)) && throw(DomainError(a,
                """Series expansion of atanh(x) diverges at x = ±1."""))
            for k in eachindex(a)
                atanh!(c, aa, r, k)
            end
            return c
        end
    end
end


# In-place functions (recursions) for Taylor1
function asinh!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    if k == 0
        @inbounds c[0] = asinh( a[0] )
        @inbounds r[0] = sqrt( a[0]^2 + 1 )
        return nothing
    end
    zero!(c, k)
    @inbounds for i in 1:k-1
        c[k] += (k-i) * r[i] * c[k-i]
    end
    # Update r[k] <- (sqrt(a^2+1))[k]
    sqr!(r, a, zero(constant_term(a)), k)
    sqrt!(r, r, r, k)
    @inbounds c[k] = (a[k] - c[k]/k) / constant_term(r)
    return nothing
end
function asinh!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
        {T<:Number}
    if k == 0
        a0 = constant_term(a)
        for j in eachindex(a0)
            asinh!(c[0], a0, r[0], j)
        end
        return nothing
    end
    zero!(c[k])
    aux = zero(c)
    @inbounds for i in 1:k-1
        # c[k] += (k-i) * r[i] * c[k-i]
        # mul_scalar!(c[k], k-i, r[i], c[k-i])
        for j in eachindex(a[k])
            mul_scalar!(aux[k], k-i, r[i], c[k-i], j)
            add!(c[k], c[k], aux[k], j)
        end
    end
    zero!(aux, k)
    sqrt!(r, a^2+1, aux, k)
    # @inbounds c[k] = (a[k] - c[k]/k) / constant_term(r)
    @inbounds for j in eachindex(a[k])
        zero!(aux[k], j)
        div!(c[k], c[k], k, j)     # c[k] <- c[k]/k
        subst!(aux[k], a[k], c[k], j) # aux <- a[k] - c[k]
    end
    @inbounds for j in eachindex(a[k])
        div!(c[k], aux[k], constant_term(r), j) # c[k] <- aux / constant_term(r)
    end
    return nothing
end

function acosh!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    if k == 0
        @inbounds c[0] = acosh( a[0] )
        @inbounds r[0] = sqrt( a[0]^2 - 1 )
        return nothing
    end
    zero!(c, k)
    @inbounds for i in 1:k-1
        c[k] += (k-i) * r[i] * c[k-i]
    end
    # Update r[k] <- (sqrt(a^2-1))[k]
    sqr!(r, a, zero(constant_term(a)), k)
    sqrt!(r, r, r, k)
    @inbounds c[k] = (a[k] - c[k]/k) / constant_term(r)
    return nothing
end
function acosh!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
        {T<:Number}
    if k == 0
        a0 = constant_term(a)
        for j in eachindex(a0)
            acosh!(c[k], a0, r[k], j)
        end
        return nothing
    end
    zero!(c[k])
    aux = zero(c)
    @inbounds for i in 1:k-1
        # c[k] += (k-i) * r[i] * c[k-i]
        for j in eachindex(a[k])
            mul_scalar!(aux[k], k-i, r[i], c[k-i], j)
            add!(c[k], c[k], aux[k], j)
        end
    end
    zero!(aux, k)
    sqrt!(r, a^2-1, aux, k)
    # @inbounds c[k] = (a[k] - c[k]/k) / constant_term(r)
    @inbounds for j in eachindex(a[k])
        zero!(aux[k], j)
        div!(c[k], c[k], k, j)     # c[k] <- c[k]/k
        subst!(aux[k], a[k], c[k], j) # aux <- a[k] - c[k]
    end
    @inbounds for j in eachindex(a[k])
        div!(c[k], aux[k], constant_term(r), j) # c[k] <- aux / constant_term(r)
    end
    return nothing
end

function atanh!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    r_coeffs = r.coeffs
    @inbounds a0 = a_coeffs[1]
    if k == 0
        @inbounds c_coeffs[1] = atanh(a0)
        @inbounds r_coeffs[1] = 1 - a0^2
        return nothing
    end
    kk = k+1
    @inbounds acc = zero(c_coeffs[kk])
    @inbounds for i in 1:k-1
        acc += (k-i) * r_coeffs[i+1] * c_coeffs[k-i+1]
    end
    @inbounds sqr!(r, a, zero(a0), k)
    @inbounds c_coeffs[kk] = (a_coeffs[kk] + acc/k) / r_coeffs[1]
    return nothing
end
function atanh!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
        {T<:Number}
    if k == 0
        a0 = constant_term(a)
        @inbounds for j in eachindex(a[k])
            atanh!(c[0], a0, r[0], j)
        end
        return nothing
    end
    zero!(c, k)
    aux = zero(c[k])
    @inbounds for i in 1:k-1
        # c[k] += (k-i) * r[i] * c[k-i]
        for j in eachindex(a[k])
            mul_scalar!(aux, k-i, r[i], c[k-i], j)
            add!(c[k], c[k], aux, j)
        end
    end
    zero!(aux)
    @inbounds sqr!(r, a, aux, k)
    # @inbounds c[k] = (a[k] + c[k]/k) / constant_term(r)
    @inbounds for j in eachindex(a[k])
        zero!(aux, j)
        div!(c[k], c[k], k, j)     # c[k] <- c[k]/k
        add!(aux, a[k], c[k], j) # aux <- a[k] + c[k]
    end
    @inbounds for j in eachindex(a[k])
        div!(c[k], aux, constant_term(r), j) # c[k] <- aux / constant_term(r)
    end
    return nothing
end

# In-place functions (recursions) for TaylorN
function asinh!(c::TaylorN{T}, a::TaylorN{T}, r::TaylorN{T}, k::Int) where
        {T<:Number}
    if k == 0
        @inbounds c[0] = asinh( a[0][1] )
        @inbounds r[0] = sqrt( a[0][1]^2 + 1 )
        return nothing
    end
    zero!(c, k)
    @inbounds for i in 1:k-1
        mul_scalar!(c[k], k-i, r[i], c[k-i])
    end
    # Update r[k] <- (sqrt(a^2+1))[k]
    sqr!(r, a, zero(constant_term(a)), k)
    sqrt!(r, r, r, k)
    # c[k] <- (a[k] - c[k]/k) / r[0], coefficient by coefficient
    r0 = constant_term(r)
    @inbounds c_k, a_k = c[k].coeffs, a[k].coeffs
    @inbounds for j in eachindex(c_k)
        c_k[j] = (a_k[j] - c_k[j]/k) / r0
    end
    return nothing
end

function acosh!(c::TaylorN{T}, a::TaylorN{T}, r::TaylorN{T}, k::Int) where
        {T<:Number}
    if k == 0
        @inbounds c[0] = acosh( a[0][1] )
        @inbounds r[0] = sqrt( a[0][1]^2 - 1 )
        return nothing
    end
    zero!(c, k)
    @inbounds for i in 1:k-1
        mul_scalar!(c[k], k-i, r[i], c[k-i])
    end
    # Update r[k] <- (sqrt(a^2-1))[k]
    sqr!(r, a, zero(constant_term(a)), k)
    sqrt!(r, r, r, k)
    # c[k] <- (a[k] - c[k]/k) / r[0], coefficient by coefficient
    r0 = constant_term(r)
    @inbounds c_k, a_k = c[k].coeffs, a[k].coeffs
    @inbounds for j in eachindex(c_k)
        c_k[j] = (a_k[j] - c_k[j]/k) / r0
    end
    return nothing
end

function atanh!(c::TaylorN{T}, a::TaylorN{T}, r::TaylorN{T}, k::Int) where
        {T<:Number}
    if k == 0
        @inbounds c[0] = atanh( a[0][1] )
        @inbounds r[0] = 1 - a[0][1]^2
        return nothing
    end
    zero!(c, k)
    @inbounds for i in 1:k-1
        mul_scalar!(c[k], k-i, r[i], c[k-i])
    end
    sqr!(r, a, zero(a[0][1]), k)
    # c[k] <- (a[k] + c[k]/k) / r[0], coefficient by coefficient
    r0 = constant_term(r)
    @inbounds c_k, a_k = c[k].coeffs, a[k].coeffs
    @inbounds for j in eachindex(c_k)
        c_k[j] = (a_k[j] + c_k[j]/k) / r0
    end
    return nothing
end


# In-place functions (recursions) for Taylor1{TaylorN}
function asinh!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        r::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
    tmp = _constant_series_like(a[0], zero(a[0][0][1]), order(a[0]))
    if k == 0
        tmp0 = zero(tmp)
        @inbounds res[0] = asinh( a[0] )
        # r[0] = sqrt(1+a[0]^2)
        for ord in eachindex(a[0])
            sqr!(r[0], a[0], tmp0[0][1], ord)
        end
        zero!(tmp0, 0)
        for ordQ in eachindex(a[0])
            one!(tmp, a[0], ordQ)
            add!(tmp, tmp, r[0], ordQ)
            # zero!(r[0], tmp, ordQ)
            sqrt!(r[0], tmp, tmp0, ordQ)
        end
        return nothing
    end
    # The recursion formula
    zero!(res[k])
    for i in 1:k-1
        @inbounds for ordQ in eachindex(a[0])
            mul_scalar!(res[k], k-i, res[k-i], r[i], ordQ)
        end
    end
    div!(res, res, k, k)
    @inbounds for ordQ in eachindex(a[0])
        subst!(tmp, a[k], res[k], ordQ)
        zero!(res[k][ordQ])
        div!(res[k], tmp, r[0], ordQ)
        zero!(tmp, ordQ)
    end
    # Update aux term r = sqrt(s), s = 1+a^2
    sqr!(r, a, tmp, k)  # r[k] <- (a^2)[k]
    sqrt!(r, r, r, k)   # r[k] <- (sqrt(s))[k]
    return nothing
end

function acosh!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        r::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
    tmp = _constant_series_like(a[0], zero(a[0][0][1]), order(a[0]))
    if k == 0
        tmp0 = zero(tmp)
        @inbounds res[0] = acosh( a[0] )
        # r[0] = sqrt(a[0]^2-1)
        # r[0] = square(a[0])
        for ord in eachindex(a[0])
            sqr!(r[0], a[0], tmp0[0][1], ord)
        end
        zero!(tmp0, 0)
        for ordQ in eachindex(a[0])
            one!(tmp, a[0], ordQ)
            subst!(tmp, r[0], tmp, ordQ)
            # zero!(r[0], tmp, ordQ)
            sqrt!(r[0], tmp, tmp0, ordQ)
        end
        return nothing
    end
    # The recursion formula
    zero!(res[k])
    for i in 1:k-1
        @inbounds for ordQ in eachindex(a[0])
            mul_scalar!(res[k], k-i, res[k-i], r[i], ordQ)
        end
    end
    div!(res, res, k, k)
    @inbounds for ordQ in eachindex(a[0])
        subst!(tmp, a[k], res[k], ordQ)
        zero!(res[k][ordQ])
        div!(res[k], tmp, r[0], ordQ)
    end
    # Update aux term r = sqrt(s), s = 1+a^2
    sqr!(r, a, tmp, k)  # r[k] <- (a^2)[k]
    sqrt!(r, r, r, k)   # r[k] <- (sqrt(s))[k]
    return nothing
end

function atanh!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        r::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
    tmp = _constant_series_like(a[0], zero(a[0][0][1]), order(a[0]))
    if k == 0
        res[0] = atanh( a[0] )
        # zero!(r, a, 0)
        sqr!(r, a, tmp, 0)
        subst!(r, one(a[0][0][1]), r, 0)
        return nothing
    end
    # The recursion formula
    zero!(res[k])
    for i in 1:k-1
        @inbounds for ordQ in eachindex(a[0])
            mul_scalar!(res[k], k-i, res[k-i], r[i], ordQ)
        end
    end
    @inbounds for ordQ in eachindex(a[0])
        # zero!(tmp, res[k], ordQ)
        # tmp[ordQ] <- res[k][ordQ] / k, in place
        tmp_hp = tmp[ordQ].coeffs
        res_hp = res[k][ordQ].coeffs
        for j in eachindex(tmp_hp)
            tmp_hp[j] = res_hp[j] / k
        end
        add!(tmp, a[k], tmp, ordQ)
        zero!(res[k][ordQ])
        div!(res[k], tmp, r[0], ordQ)
    end
    zero!(tmp)
    zero!(r[k])
    sqr!(r, a, tmp, k)
    return nothing
end


# Documentation of the recursion relations
@doc doc"""
    asinh!(c, a, r, k)

Update the `k-th` expansion coefficients `c[k+1]` of `c = asinh(a)`,
for `c` and `a` either `Taylor1` or `TaylorN`; `r = sqrt(1-c^2)` and
is passed as an argument for efficiency.

```math
\begin{equation*}
c_k = \frac{1}{ \sqrt{r_0} }
    \big( a_k - \frac{1}{k} \sum_{j=1}^{k-1} j r_{k-j} c_j \big).
\end{equation*}
```

""" asinh!


@doc doc"""
    acosh!(c, a, r, k)

Update the `k-th` expansion coefficients `c[k+1]` of `c = acosh(a)`,
for `c` and `a` either `Taylor1` or `TaylorN`; `r = sqrt(c^2-1)` and
is passed as an argument for efficiency.

```math
\begin{equation*}
c_k = \frac{1}{ r_0 }
    \big( a_k - \frac{1}{k} \sum_{j=1}^{k-1} j r_{k-j} c_j \big).
\end{equation*}
```

""" acosh!


@doc doc"""
    atanh!(c, a, r, k)

Update the `k-th` expansion coefficients `c[k+1]` of `c = atanh(a)`,
for `c` and `a` either `Taylor1` or `TaylorN`; `r = 1-a^2` and
is passed as an argument for efficiency.

```math
\begin{equation*}
c_k = \frac{1}{r_0}\big(a_k + \frac{1}{k} \sum_{j=1}^{k-1} j r_{k-j} c_j\big).
\end{equation*}
```

""" atanh!