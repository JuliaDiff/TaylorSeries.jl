# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

# Inverse trigonometric functions: asin, acos, atan


# Out-of-place functions
for T in (:Taylor1, :TaylorN)
    @eval begin
        function asin(a::$T)
            a0 = constant_term(a)
            a0^2 == one(a0) && throw(DomainError(a,
                """Series expansion of asin(x) diverges at x = ±1."""))
            order = TS.order(a)
            aux = asin(a0)
            # aa = one(aux) * a
            aa = convert($T{typeof(aux)}, a)
            c = _constant_series_like(a, aux, order)
            r = _constant_series_like(a, sqrt(1 - a0^2), order)
            for k in eachindex(a)
                asin!(c, aa, r, k)
            end
            return c
        end

        function acos(a::$T)
            a0 = constant_term(a)
            a0^2 == one(a0) && throw(DomainError(a,
                """Series expansion of asin(x) diverges at x = ±1."""))
            order = TS.order(a)
            aux = acos(a0)
            # aa = one(aux) * a
            aa = convert($T{typeof(aux)}, a)
            c = _constant_series_like(a, aux, order)
            r = _constant_series_like(a, sqrt(1 - a0^2), order)
            for k in eachindex(a)
                acos!(c, aa, r, k)
            end
            return c
        end

        function atan(a::$T)
            order = TS.order(a)
            a0 = constant_term(a)
            aux = atan(a0)
            # aa = one(aux) * a
            aa = convert($T{typeof(aux)}, a)
            c = _constant_series_like(a, aux, order)
            r = _constant_series_like(a, one(aux) + a0^2, order)
            iszero(constant_term(r)) && throw(DomainError(a,
                """Series expansion of atan(x) diverges at x = ±im."""))
            for k in eachindex(a)
                atan!(c, aa, r, k)
            end
            return c
        end

        function atan(a::$T, b::$T)
            c = atan(a/b)
            constant_term!(c, atan(constant_term(a), constant_term(b)))
            return c
        end
    end
end


# In-place functions (recursions) for Taylor1
function asin!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    if k == 0
        @inbounds c[0] = asin( a[0] )
        @inbounds r[0] = sqrt( 1 - a[0]^2 )
        return nothing
    end
    zero!(c, k)
    @inbounds for i in 1:k-1
        c[k] += (k-i) * r[i] * c[k-i]
    end
    # Compute k-th coefficient of auxiliary term s=1-a^2
    aux = zero(a)
    zero!(r, k) # r[k] <- 0
    sqr!(r, a, aux[0], k) # r[k] <- (a^2)[k]
    subst!(r, r, k) # r[k] <- -r[k]
    zero!(aux, 0)
    sqrt!(r, r, aux, k) # r[k] <- (sqrt(r))[k]
    @inbounds c[k] = (a[k] - c[k]/k) / constant_term(r)
    return nothing
end
function asin!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
        {T<:Number}
    if k == 0
        aux = constant_term(a)
        @inbounds for j in eachindex(a[k])
            asin!(c[0], aux, r[0], j) # c[0] = asin( aux )
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
    # Compute k-th coefficient of auxiliary term s=1-a^2
    zero!(r, k) # r[k] <- 0
    zero!(aux)
    sqr!(r, a, aux[0], k) # r[k] <- (a^2)[k]
    subst!(r, r, k) # r[k] <- -r[k]
    zero!(aux, 0)
    sqrt!(r, r, aux, k) # r[k] <- (sqrt(r))[k]
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

function acos!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    if k == 0
        @inbounds c[0] = acos( a[0] )
        @inbounds r[0] = sqrt( 1 - a[0]^2 )
        return nothing
    end
    zero!(c, k)
    @inbounds for i in 1:k-1
        c[k] += (k-i) * r[i] * c[k-i]
    end
    # Compute k-th coefficient of auxiliary term s=1-a^2
    aux = zero(a)
    zero!(r, k) # r[k] <- 0
    sqr!(r, a, aux[0], k) # r[k] <- (a^2)[k]
    subst!(r, r, k) # r[k] <- -r[k]
    zero!(aux)
    sqrt!(r, r, aux, k) # r[k] <- (sqrt(r))[k]
    @inbounds c[k] = -(a[k] + c[k]/k) / constant_term(r)
    return nothing
end
function acos!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
        {T<:Number}
    if k == 0
        aux = constant_term(a)
        @inbounds for j in eachindex(a[k])
            acos!(c[0], aux, r[0], j) # c[0] = acos( aux )
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
    # Compute k-th coefficient of auxiliary term s=1-a^2
    zero!(r, k) # r[k] <- 0
    zero!(aux)
    sqr!(r, a, aux[0], k) # r[k] <- (a^2)[k]
    subst!(r, r, k) # r[k] <- -r[k]
    zero!(aux, 0)
    sqrt!(r, r, aux, k) # r[k] <- (sqrt(r))[k]
    # @inbounds c[k] = - (a[k] + c[k]/k) / constant_term(r)
    @inbounds for j in eachindex(a[k])
        zero!(aux[k], j)
        div!(c[k], c[k], k, j)     # c[k] <- c[k]/k
        add!(aux[k], a[k], c[k], j) # aux <- a[k] + c[k]
        subst!(aux[k], aux[k], j)
    end
    @inbounds for j in eachindex(a[k])
        div!(c[k], aux[k], constant_term(r), j) # c[k] <- aux / constant_term(r)
    end
    return nothing
end

function atan!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    r_coeffs = r.coeffs
    if k == 0
        @inbounds a0 = a_coeffs[1]
        @inbounds c_coeffs[1] = atan(a0)
        @inbounds r_coeffs[1] = 1 + a0^2
        return nothing
    end
    kk = k+1
    @inbounds acc = zero(c_coeffs[kk])
    @inbounds for i in 1:k-1
        acc += (k-i) * r_coeffs[i+1] * c_coeffs[k-i+1]
    end
    @inbounds sqr!(r, a, zero(a_coeffs[1]), k)
    @inbounds c_coeffs[kk] = (a_coeffs[kk] - acc/k) / r_coeffs[1]
    return nothing
end
function atan!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
        {T<:Number}
    if k == 0
        aux = constant_term(a)
        for j in eachindex(a[k])
            atan!(c[0], aux, r[0], j)
        end
        return nothing
    end
    zero!(c[k])
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
    # @inbounds c[k] = (a[k] - c[k]/k) / constant_term(r)
    for j in eachindex(a[k])
        zero!(aux, j)
        div!(c[k], c[k], k, j)     # c[k] <- c[k]/k
        subst!(aux, a[k], c[k], j) # aux <- a[k] - c[k]
    end
    @inbounds for j in eachindex(a[k])
        div!(c[k], aux, constant_term(r), j) # c[k] <- aux / constant_term(r)
    end
    return nothing
end


# In-place functions (recursions) for TaylorN
function asin!(c::TaylorN{T}, a::TaylorN{T}, r::TaylorN{T}, k::Int) where
        {T<:Number}
    if k == 0
        @inbounds c[0][1] = asin( a[0][1] )
        @inbounds r[0][1] = sqrt( 1 - a[0][1]^2 )
        return nothing
    end
    zero!(c, k)
    @inbounds for i in 1:k-1
        mul_scalar!(c[k], k-i, r[i], c[k-i])
    end
    # Compute k-th coefficient of auxiliary term s=1-a^2
    zero!(r, k) # r[k] <- 0
    sqr!(r, a, zero(constant_term(a)), k) # r[k] <- (a^2)[k]
    subst!(r, r, k) # r[k] <- -r[k]
    sqrt!(r, r, r, k) # r[k] <- (sqrt(r))[k]
    for l in eachindex(c[k])
        @inbounds c[k][l] = (a[k][l] - c[k][l]/k) / constant_term(r)
    end
    return nothing
end

function acos!(c::TaylorN{T}, a::TaylorN{T}, r::TaylorN{T}, k::Int) where
        {T<:Number}
    if k == 0
        @inbounds c[0][1] = acos( a[0][1] )
        @inbounds r[0][1] = sqrt( 1 - a[0][1]^2 )
        return nothing
    end
    zero!(c, k)
    @inbounds for i in 1:k-1
        mul_scalar!(c[k], k-i, r[i], c[k-i])
    end
    # Compute k-th coefficient of auxiliary term s=1-a^2
    zero!(r, k) # r[k] <- 0
    sqr!(r, a, zero(constant_term(a)), k) # r[k] <- (a^2)[k]
    subst!(r, r, k) # r[k] <- -r[k]
    sqrt!(r, r, r, k) # r[k] <- (sqrt(r))[k]
    for l in eachindex(c[k])
        @inbounds c[k][l] = -(a[k][l] + c[k][l]/k) / constant_term(r)
    end
    return nothing
end

function atan!(c::TaylorN{T}, a::TaylorN{T}, r::TaylorN{T}, k::Int) where
        {T<:Number}
    a0 = constant_term(a)
    if k == 0
        @inbounds c[0] = atan( a0 )
        @inbounds r[0] = 1 + a0^2
        return nothing
    end
    zero!(c, k)
    @inbounds for i in 1:k-1
        mul_scalar!(c[k], k-i, r[i], c[k-i])
    end
    @inbounds sqr!(r, a, zero(a0), k)
    for l in eachindex(c[k])
        @inbounds c[k][l] = (a[k][l] - c[k][l]/k) / constant_term(r)
    end
    return nothing
end


# In-place functions (recursions) for Taylor1{TaylorN}
function asin!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        r::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
    if k == 0
        @inbounds for ordQ in eachindex(res[0])
            asin!(res[0], a[0], r[0], ordQ)
        end
        return nothing
    end
    # The recursion formula
    @inbounds zero!(res[k])
    @inbounds for i in 1:k-1
        for ordQ in eachindex(a[0])
            mul_scalar!(res[k], k-i, res[k-i], r[i], ordQ)
        end
    end
    @inbounds div!(res, res, k, k)
    @inbounds for ordQ in eachindex(a[0])
        subst!(res[k], a[k], res[k], ordQ)
        div!(res[k], r[0], ordQ)
    end
    # k-th coefficient of auxiliary term s=1-a^2
    zero!(r, k) # r[k] <- 0
    sqr!(r, a, k) # r[k] <- (a^2)[k]
    subst!(r, r, k) # r[k] <- -r[k]
    sqrt!(r, r, r, k) # r[k] <- (sqrt(r))[k]
    return nothing
end

function acos!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        r::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
    if k == 0
        @inbounds for ordQ in eachindex(res[0])
            acos!(res[0], a[0], r[0], ordQ)
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
        add!(res[k], a[k], res[k], ordQ)
        subst!(res[k], res[k], ordQ)
        div!(res[k], r[0], ordQ)
    end
    # k-th coefficient of auxiliary term s=1-a^2
    zero!(r, k) # r[k] <- 0
    sqr!(r, a, k) # r[k] <- (a^2)[k]
    subst!(r, r, k) # r[k] <- -r[k]
    sqrt!(r, r, r, k) # r[k] <- (sqrt(r))[k]
    return nothing
end

function atan!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        r::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
    aux = _constant_series_like(a[0], zero(a[0][0][1]), order(a[0]))
    if k == 0
        res[0] = atan( a[0] )
        # zero!(r, a, 0)
        sqr!(r, a, aux, 0)
        add!(r, r, one(a[0][0][1]), 0)
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
        # zero!(aux, res[k], ordQ)
        # aux[ordQ] <- - res[k][ordQ] / k, in place
        aux_hp = aux[ordQ].coeffs
        res_hp = res[k][ordQ].coeffs
        for j in eachindex(aux_hp)
            aux_hp[j] = (- res_hp[j]) / k
        end
        add!(aux, a[k], aux, ordQ)
        zero!(res[k][ordQ])
        div!(res[k], aux, r[0], ordQ)
    end
    zero!(aux)
    zero!(r[k])
    sqr!(r, a, aux, k)
    return nothing
end


# Documentation of the recursion relations
@doc doc"""
    asin!(c, a, r, k)

Update the `k-th` expansion coefficients `c[k+1]` of `c = asin(a)`,
for `c` and `a` either `Taylor1` or `TaylorN`; `r = sqrt(1-c^2)` and
is passed as an argument for efficiency.

```math
\begin{equation*}
c_k = \frac{1}{ \sqrt{r_0} }
    \big( a_k - \frac{1}{k} \sum_{j=1}^{k-1} j r_{k-j} c_j \big).
\end{equation*}
```

""" asin!


@doc doc"""
    acos!(c, a, r, k)

Update the `k-th` expansion coefficients `c[k+1]` of `c = acos(a)`,
for `c` and `a` either `Taylor1` or `TaylorN`; `r = sqrt(1-c^2)` and
is passed as an argument for efficiency.

```math
\begin{equation*}
c_k = - \frac{1}{ r_0 }
    \big( a_k - \frac{1}{k} \sum_{j=1}^{k-1} j r_{k-j} c_j \big).
\end{equation*}
```

""" acos!


@doc doc"""
    atan!(c, a, r, k)

Update the `k-th` expansion coefficients `c[k+1]` of `c = atan(a)`,
for `c` and `a` either `Taylor1` or `TaylorN`; `r = 1+a^2` and
is passed as an argument for efficiency.

```math
\begin{equation*}
c_k = \frac{1}{r_0}\big(a_k - \frac{1}{k} \sum_{j=1}^{k-1} j r_{k-j} c_j\big).
\end{equation*}
```

""" atan!