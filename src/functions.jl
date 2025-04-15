# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#


# Functions
for T in (:Taylor1, :TaylorN)
    @eval begin
        function exp(a::$T)
            order = a.order
            aux = exp(constant_term(a))
            aa = convert($T{typeof(aux)}, a)
            c = $T( aux, order )
            for k in eachindex(a)
                exp!(c, aa, k)
            end
            return c
        end

        function expm1(a::$T)
            order = a.order
            aux = expm1(constant_term(a))
            aa = convert($T{typeof(aux)}, a)
            c = $T( aux, order )
            for k in eachindex(a)
                expm1!(c, aa, k)
            end
            return c
        end

        function log(a::$T)
            iszero(constant_term(a)) && throw(DomainError(a,
                """The 0-th order coefficient must be non-zero in order to expand `log` around 0."""))
            order = a.order
            aux = log(constant_term(a))
            aa = convert($T{typeof(aux)}, a)
            c = $T( aux, order )
            for k in eachindex(a)
                log!(c, aa, k)
            end
            return c
        end
        function log1p(a::$T)
            # constant_term(a) < -one(constant_term(a)) && throw(DomainError(a,
            #         """The 0-th order coefficient must be larger than -1 in order to expand `log1`."""))
            order = a.order
            aux = log1p(constant_term(a))
            # aa = one(aux) * a
            aa = convert($T{typeof(aux)}, a)
            c = $T( aux, order )
            for k in eachindex(a)
                log1p!(c, aa, k)
            end
            return c
        end

        sin(a::$T) = sincos(a)[1]
        cos(a::$T) = sincos(a)[2]
        function sincos(a::$T)
            order = a.order
            aux, auxc = sincos(constant_term(a))
            # aa = one(aux) * a
            aa = convert($T{typeof(aux)}, a)
            s = $T( aux, order )
            c = $T( auxc, order )
            for k in eachindex(a)
                sincos!(s, c, aa, k)
            end
            return s, c
        end

        sinpi(a::$T) = sincospi(a)[1]
        cospi(a::$T) = sincospi(a)[2]
        function sincospi(a::$T)
            order = a.order
            aux, auxc = sincospi(constant_term(a))
            aa = one(aux) * a
            # aa = convert($T{typeof(aux)}, a)
            s = $T( aux, order )
            c = $T( auxc, order )
            for k in eachindex(a)
                sincospi!(s, c, aa, k)
            end
            return s, c
        end

        function tan(a::$T)
            order = a.order
            aux = tan(constant_term(a))
            # aa = one(aux) * a
            aa = convert($T{typeof(aux)}, a)
            c = $T(aux, order)
            c2 = $T(aux^2, order)
            for k in eachindex(a)
                tan!(c, aa, c2, k)
            end
            return c
        end

        function asin(a::$T)
            a0 = constant_term(a)
            a0^2 == one(a0) && throw(DomainError(a,
                """Series expansion of asin(x) diverges at x = ±1."""))
            order = a.order
            aux = asin(a0)
            # aa = one(aux) * a
            aa = convert($T{typeof(aux)}, a)
            c = $T( aux, order )
            r = $T( sqrt(1 - a0^2), order )
            for k in eachindex(a)
                asin!(c, aa, r, k)
            end
            return c
        end

        function acos(a::$T)
            a0 = constant_term(a)
            a0^2 == one(a0) && throw(DomainError(a,
                """Series expansion of asin(x) diverges at x = ±1."""))
            order = a.order
            aux = acos(a0)
            # aa = one(aux) * a
            aa = convert($T{typeof(aux)}, a)
            c = $T( aux, order )
            r = $T( sqrt(1 - a0^2), order )
            for k in eachindex(a)
                acos!(c, aa, r, k)
            end
            return c
        end

        function atan(a::$T)
            order = a.order
            a0 = constant_term(a)
            aux = atan(a0)
            # aa = one(aux) * a
            aa = convert($T{typeof(aux)}, a)
            c = $T( aux, order)
            r = $T(one(aux) + a0^2, order)
            iszero(constant_term(r)) && throw(DomainError(a,
                """Series expansion of atan(x) diverges at x = ±im."""))
            for k in eachindex(a)
                atan!(c, aa, r, k)
            end
            return c
        end

        function atan(a::$T, b::$T)
            c = atan(a/b)
            c[0] = atan(constant_term(a), constant_term(b))
            return c
        end

        sinh(a::$T) = sinhcosh(a)[1]
        cosh(a::$T) = sinhcosh(a)[2]
        function sinhcosh(a::$T)
            order = a.order
            aux = sinh(constant_term(a))
            # aa = one(aux) * a
            aa = convert($T{typeof(aux)}, a)
            s = $T( aux, order)
            c = $T( cosh(constant_term(a)), order)
            for k in eachindex(a)
                sinhcosh!(s, c, aa, k)
            end
            return s, c
        end

        function tanh(a::$T)
            order = a.order
            aux = tanh( constant_term(a) )
            # aa = one(aux) * a
            aa = convert($T{typeof(aux)}, a)
            c = $T( aux, order)
            c2 = $T( aux^2, order)
            for k in eachindex(a)
                tanh!(c, aa, c2, k)
            end
            return c
        end


        function asinh(a::$T)
            order = a.order
            a0 = constant_term(a)
            aux = asinh(a0)
            # aa = one(aux) * a
            aa = convert($T{typeof(aux)}, a)
            c = $T( aux, order )
            r = $T( sqrt(a0^2 + 1), order )
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
            order = a.order
            aux = acosh(a0)
            # aa = one(aux) * a
            aa = convert($T{typeof(aux)}, a)
            c = $T( aux, order )
            r = $T( sqrt(a0^2 - 1), order )
            for k in eachindex(a)
                acosh!(c, aa, r, k)
            end
            return c
        end

        function atanh(a::$T)
            order = a.order
            a0 = constant_term(a)
            aux = atanh(a0)
            # aa = one(aux) * a
            aa = convert($T{typeof(aux)}, a)
            c = $T( aux, order)
            r = $T(one(aux) - a0^2, order)
            iszero(constant_term(r)) && throw(DomainError(a,
                """Series expansion of atanh(x) diverges at x = ±1."""))
            for k in eachindex(a)
                atanh!(c, aa, r, k)
            end
            return c
        end

    end
end


# Recursive functions (homogeneous coefficients)
@inline function zero!(a::Taylor1{T}, k::Int) where {T<:NumberNotSeries}
    a[k] = zero(a[k])
    return nothing
end

@inline function zero!(a::Taylor1{T}, k::Int) where {T<:Number}
    for l in eachindex(a[k])
        zero!(a[k], l)
    end
    return nothing
end

@inline function zero!(a::Taylor1{T}) where {T<:Number}
    for k in eachindex(a)
        zero!(a, k)
    end
    return nothing
end

@inline function zero!(a::HomogeneousPolynomial{T}, k::Int) where {T<:Number}
    a[k] = zero(a[k])
    return nothing
end

@inline function zero!(a::HomogeneousPolynomial{T}) where {T<:Number}
    for k in eachindex(a)
        zero!(a, k)
    end
    return nothing
end

@inline function zero!(a::TaylorN{T}, k::Int) where {T<:Number}
    zero!(a[k])
    return nothing
end

@inline function zero!(a::TaylorN{T}) where {T<:Number}
    for k in eachindex(a)
        zero!(a, k)
    end
    return nothing
end

@inline function one!(a::Taylor1{T}, k::Int) where {T<:NumberNotSeries}
    k == 0 ? a[k] = one(a[k]) : a[k] = zero(a[k])
    return nothing
end

@inline function one!(a::Taylor1{T}, k::Int) where {T<:Number}
    for j in eachindex(a[k])
        one!(a[k], j)
    end
    return nothing
end

@inline function one!(a::Taylor1{T}) where {T<:Number}
    for k in eachindex(a)
        one!(a, k)
    end
    return nothing
end

@inline function one!(c::TaylorN{T}, k::Int) where {T<:Number}
    zero!(c, k)
    if k == 0
        @inbounds c[0][1] = one(constant_term(c[0][1]))
    end
    return nothing
end

@inline function identity!(c::HomogeneousPolynomial{T}, a::HomogeneousPolynomial{T},
        k::Int) where {T<:Number}
    @inbounds c[k] = identity(a[k])
    return nothing
end

@inline function identity!(c::HomogeneousPolynomial{Taylor1{T}},
        a::HomogeneousPolynomial{Taylor1{T}}, k::Int) where {T<:NumberNotSeries}
    @inbounds for l in eachindex(c[k])
        identity!(c[k], a[k], l)
    end
    return nothing
end

# Taylor1 (including nested Taylor1s)
@inline function one!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where {T<:NumberNotSeries}
    zero!(c, k)
    (k == 0) && (@inbounds c[0] = one(constant_term(a)))
    return nothing
end
@inline function one!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where {T<:Number}
    for i in eachindex(a[k])
        one!(c[k], a[k], i)
    end
    return nothing
end

@inline function identity!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    @inbounds c[k] = identity(a[k])
    return nothing
end
@inline function identity!(c::Taylor1{T}, a::T, k::Int) where {T<:Number}
    zero!(c[k])
    if k == 0
        @inbounds c[k] = copy(a)
    end
    return nothing
end

@inline function exp!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    if k == 0
        @inbounds c[0] = exp(constant_term(a))
        return nothing
    end
    zero!(c, k)
    @inbounds for i = 0:k-1
        c[k] += (k-i) * a[k-i] * c[i]
    end
    div!(c, c, k, k)
    return nothing
end
@inline function exp!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:Number}
    aux = constant_term(a)
    if k == 0
        for j in eachindex(a[k])
            exp!(c[k], aux, j)
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
    div!(c, c, k, k)
    return nothing
end

@inline function expm1!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    if k == 0
        @inbounds c[0] = expm1(constant_term(a))
        return nothing
    end
    exp!(c, a, k)
    return nothing
end
@inline function expm1!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:Number}
    if k == 0
        for j in eachindex(a[k])
            expm1!(c[k], constant_term(a), j)
        end
        return nothing
    end
    exp!(c, a, k)
    return nothing
end

@inline function log!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    if k == 0
        @inbounds c[0] = log(constant_term(a))
        return nothing
    elseif k == 1
        @inbounds c[1] = a[1] / constant_term(a)
        return nothing
    end
    zero!(c, k)
    @inbounds for i = 1:k-1
        c[k] += (k-i) * a[i] * c[k-i]
    end
    @inbounds c[k] = (a[k] - c[k]/k) / constant_term(a)
    return nothing
end
@inline function log!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:Number}
    aux = constant_term(a)
    if k == 0
        # @inbounds c[0] = log(constant_term(a))
        for j in eachindex(a[k])
            log!(c[0], aux, j)
        end
        return nothing
    elseif k == 1
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

@inline function log1p!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    a0 = constant_term(a)
    a0p1 = a0+one(a0)
    if k == 0
        @inbounds c[0] = log1p(a0)
        return nothing
    elseif k == 1
        @inbounds c[1] = a[1] / a0p1
        return nothing
    end
    zero!(c, k)
    @inbounds for i = 1:k-1
        c[k] += (k-i) * a[i] * c[k-i]
    end
    @inbounds c[k] = (a[k] - c[k]/k) / a0p1
    return nothing
end
@inline function log1p!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where
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

@inline function sincos!(s::Taylor1{T}, c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    if k == 0
        @inbounds s[0], c[0] = sincos( constant_term(a) )
        return nothing
    end
    zero!(s, k)
    zero!(c, k)
    @inbounds for i = 1:k
        x = i * a[i]
        s[k] += x * c[k-i]
        c[k] -= x * s[k-i]
    end
    s[k] = s[k] / k
    c[k] = c[k] / k
    return nothing
end
@inline function sincos!(s::Taylor1{T}, c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:Number}
    if k == 0
        for j in eachindex(a[k])
            @inbounds sincos!(s[0], c[0], constant_term(a), j)
        end
        return nothing
    end
    x = zero(c[k])
    for j in eachindex(a[k])
        zero!(s[k], j)
        zero!(c[k], j)
        @inbounds for i = 1:k
            # x = i * a[i]
            mul!(x, a[i], i, j)
            # s[k] += x * c[k-i]
            muladd!(s[k], x, c[k-i], j)
            subst!(x, x, j)
            # c[k] -= x * s[k-i]
            muladd!(c[k], x, s[k-i], j)
            subst!(x, x, j)
        end
    end
    for j in eachindex(a[k])
        div!(s[k], s[k], k, j)
        div!(c[k], c[k], k, j)
    end
    return nothing
end

@inline function sincospi!(s::Taylor1{T}, c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    if k == 0
        a0 = constant_term(a)
        @inbounds s[0], c[0] = sincospi( a0 )
        return nothing
    end
    mul!(a, pi, a, k)
    sincos!(s, c, a, k)
    return nothing
end
@inline function sincospi!(s::Taylor1{T}, c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:Number}
    if k == 0
        a0 = constant_term(a)
        for j in eachindex(a[k])
            @inbounds sincospi!(s[0], c[0], a0, j)
        end
        return nothing
    end
    aux = zero(a)
    mul!(aux, a, pi)
    sincos!(s, c, aux, k)
    return nothing
end

@inline function tan!(c::Taylor1{T}, a::Taylor1{T}, c2::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    if k == 0
        @inbounds aux = tan( constant_term(a) )
        @inbounds c[0] = aux
        @inbounds c2[0] = aux^2
        return nothing
    end
    zero!(c, k)
    @inbounds for i = 0:k-1
        c[k] += (k-i) * a[k-i] * c2[i]
    end
    # c[k] <- c[k]/k
    div!(c, c, k, k)
    # c[k] <- c[k] + a[k]
    add!(c, a, c, k)
    sqr!(c2, c, k)
    return nothing
end
@inline function tan!(c::Taylor1{T}, a::Taylor1{T}, c2::Taylor1{T}, k::Int) where
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
    end
    sqr!(c2, c, k)
    return nothing
end

@inline function asin!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    if k == 0
        a0 = constant_term(a)
        @inbounds c[0] = asin( a0 )
        @inbounds r[0] = sqrt( 1 - a0^2 )
        return nothing
    end
    zero!(c, k)
    @inbounds for i in 1:k-1
        c[k] += (k-i) * r[i] * c[k-i]
    end
    # Compute k-th coefficient of auxiliary term s=1-a^2
    zero!(r, k) # r[k] <- 0
    sqr!(r, a, k) # r[k] <- (a^2)[k]
    subst!(r, r, k) # r[k] <- -r[k]
    sqrt!(r, r, k) # r[k] <- (sqrt(r))[k]
    @inbounds c[k] = (a[k] - c[k]/k) / constant_term(r)
    return nothing
end
@inline function asin!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
        {T<:Number}
    if k == 0
        aux = constant_term(a)
        @inbounds for j in eachindex(a[k])
            asin!(c[0], aux, r[0], j) # c[0] = asin( aux )
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
    # Compute k-th coefficient of auxiliary term s=1-a^2
    zero!(r, k) # r[k] <- 0
    sqr!(r, a, k) # r[k] <- (a^2)[k]
    subst!(r, r, k) # r[k] <- -r[k]
    sqrt!(r, r, k) # r[k] <- (sqrt(r))[k]
    # @inbounds c[k] = (a[k] - c[k]/k) / constant_term(r)
    @inbounds for j in eachindex(a[k])
        zero!(aux, j)
        div!(c[k], c[k], k, j)     # c[k] <- c[k]/k
        subst!(aux, a[k], c[k], j) # aux <- a[k] - c[k]
    end
    @inbounds for j in eachindex(a[k])
        div!(c[k], aux, constant_term(r), j) # c[k] <- aux / constant_term(r)
    end
    return nothing
end

@inline function acos!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    if k == 0
        a0 = constant_term(a)
        @inbounds c[0] = acos( a0 )
        @inbounds r[0] = sqrt( 1 - a0^2 )
        return nothing
    end
    zero!(c, k)
    @inbounds for i in 1:k-1
        c[k] += (k-i) * r[i] * c[k-i]
    end
    # Compute k-th coefficient of auxiliary term s=1-a^2
    zero!(r, k) # r[k] <- 0
    sqr!(r, a, k) # r[k] <- (a^2)[k]
    subst!(r, r, k) # r[k] <- -r[k]
    sqrt!(r, r, k) # r[k] <- (sqrt(r))[k]
    @inbounds c[k] = -(a[k] + c[k]/k) / constant_term(r)
    # end
    return nothing
end
@inline function acos!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
        {T<:Number}
    if k == 0
        aux = constant_term(a)
        @inbounds for j in eachindex(a[k])
            acos!(c[0], aux, r[0], j) # c[0] = acos( aux )
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
    # Compute k-th coefficient of auxiliary term s=1-a^2
    zero!(r, k) # r[k] <- 0
    sqr!(r, a, k) # r[k] <- (a^2)[k]
    subst!(r, r, k) # r[k] <- -r[k]
    sqrt!(r, r, k) # r[k] <- (sqrt(r))[k]
    # @inbounds c[k] = - (a[k] + c[k]/k) / constant_term(r)
    @inbounds for j in eachindex(a[k])
        zero!(aux, j)
        div!(c[k], c[k], k, j)     # c[k] <- c[k]/k
        add!(aux, a[k], c[k], j) # aux <- a[k] + c[k]
        subst!(aux, aux, j)
    end
    @inbounds for j in eachindex(a[k])
        div!(c[k], aux, constant_term(r), j) # c[k] <- aux / constant_term(r)
    end
    return nothing
end

@inline function atan!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    if k == 0
        a0 = constant_term(a)
        @inbounds c[0] = atan( a0 )
        @inbounds r[0] = 1 + a0^2
        return nothing
    end
    zero!(c, k)
    @inbounds for i in 1:k-1
        c[k] += (k-i) * r[i] * c[k-i]
    end
    @inbounds sqr!(r, a, k)
    @inbounds c[k] = (a[k] - c[k]/k) / constant_term(r)
    return nothing
end
@inline function atan!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
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
    @inbounds sqr!(r, a, k)
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

@inline function sinhcosh!(s::Taylor1{T}, c::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    if k == 0
        @inbounds s[0] = sinh( constant_term(a) )
        @inbounds c[0] = cosh( constant_term(a) )
        return nothing
    end
    x = a[1]
    zero!(s, k)
    zero!(c, k)
    @inbounds for i = 1:k
        x = i * a[i]
        s[k] += x * c[k-i]
        c[k] += x * s[k-i]
    end
    @inbounds div!(s, s, k, k)
    @inbounds div!(c, c, k, k)
    return nothing
end
@inline function sinhcosh!(s::Taylor1{T}, c::Taylor1{T}, a::Taylor1{T}, k::Int) where
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
            mul!(x, a[i], i, j)
            muladd!(s[k], x, c[k-i], j)
            muladd!(c[k], x, s[k-i], j)
        end
    end
    for j in eachindex(a[k])
        div!(s[k], s[k], k, j)
        div!(c[k], c[k], k, j)
    end
    return nothing
end

@inline function tanh!(c::Taylor1{T}, a::Taylor1{T}, c2::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    if k == 0
        @inbounds aux = tanh( constant_term(a) )
        @inbounds c[0] = aux
        @inbounds c2[0] = aux^2
        return nothing
    end
    zero!(c, k)
    @inbounds for i = 0:k-1
        c[k] += (k-i) * a[k-i] * c2[i]
    end
    @inbounds c[k] = a[k] - c[k]/k
    sqr!(c2, c, k)
    return nothing
end
@inline function tanh!(c::Taylor1{T}, a::Taylor1{T}, c2::Taylor1{T}, k::Int) where
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
    sqr!(c2, c, k)
    return nothing
end

@inline function asinh!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    if k == 0
        a0 = constant_term(a)
        @inbounds c[0] = asinh( a0 )
        @inbounds r[0] = sqrt( a0^2 + 1 )
        return nothing
    end
    zero!(c, k)
    @inbounds for i in 1:k-1
        c[k] += (k-i) * r[i] * c[k-i]
    end
    sqrt!(r, a^2+1, k)
    @inbounds c[k] = (a[k] - c[k]/k) / constant_term(r)
    return nothing
end
@inline function asinh!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
        {T<:Number}
    if k == 0
        a0 = constant_term(a)
        for j in eachindex(a0)
            asinh!(c[0], a0, r[0], j)
        end
        return nothing
    end
    zero!(c[k])
    aux = zero(c[k])
    @inbounds for i in 1:k-1
        # c[k] += (k-i) * r[i] * c[k-i]
        # mul_scalar!(c[k], k-i, r[i], c[k-i])
        for j in eachindex(a[k])
            mul_scalar!(aux, k-i, r[i], c[k-i], j)
            add!(c[k], c[k], aux, j)
        end
    end
    sqrt!(r, a^2+1, k)
    # @inbounds c[k] = (a[k] - c[k]/k) / constant_term(r)
    @inbounds for j in eachindex(a[k])
        zero!(aux, j)
        div!(c[k], c[k], k, j)     # c[k] <- c[k]/k
        subst!(aux, a[k], c[k], j) # aux <- a[k] - c[k]
    end
    @inbounds for j in eachindex(a[k])
        div!(c[k], aux, constant_term(r), j) # c[k] <- aux / constant_term(r)
    end
    return nothing
end

@inline function acosh!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    if k == 0
        a0 = constant_term(a)
        @inbounds c[0] = acosh( a0 )
        @inbounds r[0] = sqrt( a0^2 - 1 )
        return nothing
    end
    zero!(c, k)
    @inbounds for i in 1:k-1
        c[k] += (k-i) * r[i] * c[k-i]
    end
    sqrt!(r, a^2-1, k)
    @inbounds c[k] = (a[k] - c[k]/k) / constant_term(r)
    return nothing
end
@inline function acosh!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
        {T<:Number}
    if k == 0
        a0 = constant_term(a)
        for j in eachindex(a0)
            acosh!(c[k], a0, r[k], j)
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
    sqrt!(r, a^2-1, k)
    # @inbounds c[k] = (a[k] - c[k]/k) / constant_term(r)
    @inbounds for j in eachindex(a[k])
        zero!(aux, j)
        div!(c[k], c[k], k, j)     # c[k] <- c[k]/k
        subst!(aux, a[k], c[k], j) # aux <- a[k] - c[k]
    end
    @inbounds for j in eachindex(a[k])
        div!(c[k], aux, constant_term(r), j) # c[k] <- aux / constant_term(r)
    end
    return nothing
end

@inline function atanh!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    if k == 0
        a0 = constant_term(a)
        @inbounds c[0] = atanh( a0 )
        @inbounds r[0] = 1 - a0^2
        return nothing
    end
    zero!(c, k)
    @inbounds for i in 1:k-1
        c[k] += (k-i) * r[i] * c[k-i]
    end
    @inbounds sqr!(r, a, k)
    @inbounds c[k] = (a[k] + c[k]/k) / constant_term(r)
    return nothing
end
@inline function atanh!(c::Taylor1{T}, a::Taylor1{T}, r::Taylor1{T}, k::Int) where
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
    @inbounds sqr!(r, a, k)
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



#---
@inline function identity!(c::TaylorN{T}, a::TaylorN{T}, k::Int) where
        {T<:NumberNotSeries}
    @inbounds for l in eachindex(c[k])
        identity!(c[k], a[k], l)
    end
    return nothing
end

@inline function one!(c::TaylorN{T}, a::TaylorN{T}, k::Int) where {T<:Number}
    zero!(c, k)
    (k == 0) && (@inbounds c[0][1] = one(constant_term(a)))
    return nothing
end

@inline function exp!(c::TaylorN{T}, a::TaylorN{T}, k::Int) where {T<:Number}
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

@inline function expm1!(c::TaylorN{T}, a::TaylorN{T}, k::Int) where {T<:Number}
    if k == 0
        @inbounds c[0] = expm1(constant_term(a))
        return nothing
    end
    zero!(c, k)
    c0 = c[0]+one(c[0])
    @inbounds mul_scalar!(c[k], k, a[k], c0)
    @inbounds for i = 1:k-1
        mul_scalar!(c[k], k-i, a[k-i], c[i])
    end
    @inbounds div!(c[k], c[k], k)
    return nothing
end

@inline function log!(c::TaylorN{T}, a::TaylorN{T}, k::Int) where {T<:Number}
    if k == 0
        @inbounds c[0] = log(constant_term(a))
        return nothing
    elseif k == 1
        @inbounds c[1] = a[1] / constant_term(a)
        return nothing
    end
    zero!(c, k)
    @inbounds for i = 1:k-1
        mul_scalar!(c[k], k-i, a[i], c[k-i])
    end
    @inbounds c[k] = (a[k] - c[k]/k) / constant_term(a)
    return nothing
end

@inline function log1p!(c::TaylorN{T}, a::TaylorN{T}, k::Int) where {T<:Number}
    if k == 0
        a0 = constant_term(a)
        @inbounds c[0] = log1p(a0)
        return nothing
    elseif k == 1
        a0 = constant_term(a)
        a0p1 = a0+one(a0)
        @inbounds c[1] = a[1] / a0p1
        return nothing
    end
    a0 = constant_term(a)
    a0p1 = a0+one(a0)
    zero!(c, k)
    @inbounds for i = 1:k-1
        mul_scalar!(c[k], k-i, a[i], c[k-i])
    end
    @inbounds c[k] = (a[k] - c[k]/k) / a0p1
    return nothing
end

@inline function sincos!(s::TaylorN{T}, c::TaylorN{T}, a::TaylorN{T}, k::Int) where
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

@inline function sincospi!(s::TaylorN{T}, c::TaylorN{T}, a::TaylorN{T}, k::Int) where
        {T<:Number}
    if k == 0
        a0 = constant_term(a)
        @inbounds s[0][1], c[0][1] = sincospi( a0 )
        return nothing
    end
    mul!(a, pi, a, k)
    sincos!(s, c, a, k)
    return nothing
end

@inline function tan!(c::TaylorN{T}, a::TaylorN{T}, c2::TaylorN{T}, k::Int) where
        {T<:Number}
    if k == 0
        @inbounds aux = tan( constant_term(a) )
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
    sqr!(c2, c, k)
    return nothing
end

@inline function asin!(c::TaylorN{T}, a::TaylorN{T}, r::TaylorN{T}, k::Int) where
        {T<:Number}
    if k == 0
        a0 = constant_term(a)
        @inbounds c[0][1] = asin( a0 )
        @inbounds r[0][1] = sqrt( 1 - a0^2 )
        return nothing
    end
    zero!(c, k)
    @inbounds for i in 1:k-1
        mul_scalar!(c[k], k-i, r[i], c[k-i])
    end
    # Compute k-th coefficient of auxiliary term s=1-a^2
    zero!(r, k) # r[k] <- 0
    sqr!(r, a, k) # r[k] <- (a^2)[k]
    subst!(r, r, k) # r[k] <- -r[k]
    sqrt!(r, r, k) # r[k] <- (sqrt(r))[k]
    for l in eachindex(c[k])
        @inbounds c[k][l] = (a[k][l] - c[k][l]/k) / constant_term(r)
    end
    return nothing
end

@inline function acos!(c::TaylorN{T}, a::TaylorN{T}, r::TaylorN{T}, k::Int) where
        {T<:Number}
    if k == 0
        a0 = constant_term(a)
        @inbounds c[0][1] = acos( a0 )
        @inbounds r[0][1] = sqrt( 1 - a0^2 )
        return nothing
    end
    zero!(c, k)
    @inbounds for i in 1:k-1
        mul_scalar!(c[k], k-i, r[i], c[k-i])
    end
    # Compute k-th coefficient of auxiliary term s=1-a^2
    zero!(r, k) # r[k] <- 0
    sqr!(r, a, k) # r[k] <- (a^2)[k]
    subst!(r, r, k) # r[k] <- -r[k]
    sqrt!(r, r, k) # r[k] <- (sqrt(r))[k]
    for l in eachindex(c[k])
        @inbounds c[k][l] = -(a[k][l] + c[k][l]/k) / constant_term(r)
    end
    return nothing
end

@inline function atan!(c::TaylorN{T}, a::TaylorN{T}, r::TaylorN{T}, k::Int) where
        {T<:Number}
    if k == 0
        a0 = constant_term(a)
        @inbounds c[0] = atan( a0 )
        @inbounds r[0] = 1 + a0^2
        return nothing
    end
    zero!(c, k)
    @inbounds for i in 1:k-1
        mul_scalar!(c[k], k-i, r[i], c[k-i])
    end
    @inbounds sqr!(r, a, k)
    for l in eachindex(c[k])
        @inbounds c[k][l] = (a[k][l] - c[k][l]/k) / constant_term(r)
    end
    return nothing
end

@inline function sinhcosh!(s::TaylorN{T}, c::TaylorN{T}, a::TaylorN{T}, k::Int) where
        {T<:Number}
    if k == 0
        @inbounds s[0] = sinh( constant_term(a) )
        @inbounds c[0] = cosh( constant_term(a) )
        return nothing
    end
    x = a[1]
    zero!(s, k)
    zero!(c, k)
    @inbounds for i = 1:k
        mul_scalar!(s[k], i, a[i], c[k-i])
        mul_scalar!(c[k], i, a[i], s[k-i])
    end
    @inbounds div!(s[k], s[k], k)
    @inbounds div!(c[k], c[k], k)
    return nothing
end

@inline function tanh!(c::TaylorN{T}, a::TaylorN{T}, c2::TaylorN{T}, k::Int) where
        {T<:Number}
    if k == 0
        @inbounds aux = tanh( constant_term(a) )
        @inbounds c[0] = aux
        @inbounds c2[0] = aux^2
        return nothing
    end
    zero!(c, k)
    @inbounds for i = 0:k-1
        mul_scalar!(c[k], k-i, a[k-i], c2[i])
    end
    @inbounds for l in eachindex(c[k])
        c[k][l] = a[k][l] - c[k][l]/k
    end
    sqr!(c2, c, k)
    return nothing
end

@inline function asinh!(c::TaylorN{T}, a::TaylorN{T}, r::TaylorN{T}, k::Int) where
        {T<:Number}
    if k == 0
        a0 = constant_term(a)
        @inbounds c[0] = asinh( a0 )
        @inbounds r[0] = sqrt( a0^2 + 1 )
        return nothing
    end
    zero!(c, k)
    @inbounds for i in 1:k-1
        mul_scalar!(c[k], k-i, r[i], c[k-i])
    end
    sqrt!(r, a^2+1, k)
    @inbounds c[k] = (a[k] - c[k]/k) / constant_term(r)
    return nothing
end

@inline function acosh!(c::TaylorN{T}, a::TaylorN{T}, r::TaylorN{T}, k::Int) where
        {T<:Number}
    if k == 0
        a0 = constant_term(a)
        @inbounds c[0] = acosh( a0 )
        @inbounds r[0] = sqrt( a0^2 - 1 )
        return nothing
    end
    zero!(c, k)
    @inbounds for i in 1:k-1
        mul_scalar!(c[k], k-i, r[i], c[k-i])
    end
    sqrt!(r, a^2-1, k)
    @inbounds c[k] = (a[k] - c[k]/k) / constant_term(r)
    return nothing
end

@inline function atanh!(c::TaylorN{T}, a::TaylorN{T}, r::TaylorN{T}, k::Int) where
        {T<:Number}
    if k == 0
        a0 = constant_term(a)
        @inbounds c[0] = atanh( a0 )
        @inbounds r[0] = 1 - a0^2
        return nothing
    end
    zero!(c, k)
    @inbounds for i in 1:k-1
        mul_scalar!(c[k], k-i, r[i], c[k-i])
    end
    @inbounds sqr!(r, a, k)
    @inbounds c[k] = (a[k] + c[k]/k) / constant_term(r)
    return nothing
end

#---

for T in (:Taylor1, :TaylorN)
    @eval begin
        @inline function identity!(c::$T{T}, a::$T{T}, k::Int) where {T<:Number}
            @inbounds for l in eachindex(c[k])
                identity!(c[k], a[k], l)
            end
            return nothing
        end

        @inline function zero!(c::$T{T}, a::$T{T}, k::Int) where {T<:Number}
            zero!(c, k)
            return nothing
        end

        @inline function abs!(c::$T{T}, a::$T{T}, k::Int) where {T<:Number}
            z = zero(constant_term(a))
            if constant_term(constant_term(a)) > constant_term(z)
                return add!(c, a, k)
            elseif constant_term(constant_term(a)) < constant_term(z)
                return subst!(c, a, k)
            else
                throw(DomainError(a,
                    """The 0th order coefficient must be non-zero
                    (abs(x) is not differentiable at x=0)."""))
            end
            return nothing
        end

        @inline abs2!(c::$T{T}, a::$T{T}, k::Int) where {T<:Number} = sqr!(c, a, k)

    end
end


# Mutating functions for Taylor1{TaylorN{T}}
@inline function exp!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
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

@inline function expm1!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        k::Int) where {T<:NumberNotSeries}
    if k == 0
        @inbounds for ordQ in eachindex(a[0])
            # zero!(res[0], a[0], ordQ)
            expm1!(res[0], a[0], ordQ)
        end
        return nothing
    end
    # The recursion formula
    tmp = TaylorN( zero(a[k][0][1]), a[0].order)
    zero!(res[k])
    # i=0 term of sum
    @inbounds for ordQ in eachindex(a[0])
        one!(tmp, a[0], ordQ)
        add!(tmp, res[0], tmp, ordQ)
        tmp[ordQ] = k * tmp[ordQ]
        # zero!(res[k], a[0], ordQ)
        mul!(res[k], a[k], tmp, ordQ)
    end
    for i = 1:k-1
        @inbounds for ordQ in eachindex(a[0])
            tmp[ordQ] = (k-i) * res[i][ordQ]
            mul!(res[k], tmp, a[k-i], ordQ)
        end
    end
    div!(res, res, k, k)
    return nothing
end

@inline function log!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        k::Int) where {T<:NumberNotSeries}
    if k == 0
        @inbounds for ordQ in eachindex(a[0])
            # zero!(res[k], a[0], ordQ)
            log!(res[k], a[0], ordQ)
        end
        return nothing
    elseif k == 1
        @inbounds for ordQ in eachindex(a[0])
            zero!(res[k][ordQ])
            div!(res[k], a[1], a[0], ordQ)
        end
        return nothing
    end
    # The recursion formula
    tmp = TaylorN( zero(a[k][0][1]), a[0].order)
    zero!(res[k])
    for i = 1:k-1
        @inbounds for ordQ in eachindex(a[0])
            tmp[ordQ] = (k-i) * res[k-i][ordQ]
            mul!(res[k], tmp, a[i], ordQ)
        end
    end
    div!(res, res, k, k)
    @inbounds for ordQ in eachindex(a[0])
        subst!(tmp, a[k], res[k], ordQ)
        zero!(res[k][ordQ])
        div!(res[k], tmp, a[0], ordQ)
    end
    return nothing
end

@inline function log1p!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        k::Int) where {T<:NumberNotSeries}
    if k == 0
        @inbounds for ordQ in eachindex(a[0])
            # zero!(res[k], a[0], ordQ)
            log1p!(res[k], a[0], ordQ)
        end
        return nothing
    end
    tmp1 = TaylorN( zero(a[k][0][1]), a[0].order)
    zero!(res[k])
    @inbounds for ordQ in eachindex(a[0])
        # zero!(res[k], a[0], ordQ)
        one!(tmp1, a[0], ordQ)
        add!(tmp1, tmp1, a[0], ordQ)
    end
    if k == 1
        @inbounds for ordQ in eachindex(a[0])
            div!(res[k], a[1], tmp1, ordQ)
        end
        return nothing
    end
    # The recursion formula
    tmp = TaylorN( zero(a[k][0][1]), a[0].order)
    for i = 1:k-1
        @inbounds for ordQ in eachindex(a[0])
            tmp[ordQ] = (k-i) * res[k-i][ordQ]
            mul!(res[k], tmp, a[i], ordQ)
        end
    end
    div!(res, res, k, k)
    @inbounds for ordQ in eachindex(a[0])
        subst!(tmp, a[k], res[k], ordQ)
        zero!(res[k][ordQ])
        div!(res[k], tmp, tmp1, ordQ)
    end
    return nothing
end

@inline function sincos!(s::Taylor1{TaylorN{T}}, c::Taylor1{TaylorN{T}},
        a::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
    if k == 0
        @inbounds for ordQ in eachindex(a[0])
            sincos!(s[0], c[0], a[0], ordQ)
        end
        return nothing
    end
    # The recursion formula
    # x = TaylorN( a[1][0][1], a[0].order )
    zero!(s[k])
    zero!(c[k])
    @inbounds for i = 1:k
        for ordQ in eachindex(a[0])
            # x[ordQ].coeffs .= i .* a[i][ordQ].coeffs
            mul_scalar!(s[k], i, a[i], c[k-i], ordQ)
            mul_scalar!(c[k], i, a[i], s[k-i], ordQ)
        end
    end
    div!(s, s, k, k)
    subst!(c, c, k)
    div!(c, c, k, k)
    return nothing
end

@inline function sincospi!(s::Taylor1{TaylorN{T}}, c::Taylor1{TaylorN{T}},
        a::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
    if k == 0
        @inbounds for ordQ in eachindex(a[0])
            sincospi!(s[0], c[0], a[0], ordQ)
        end
        return nothing
    end
    # aa = pi * a
    aa = Taylor1(zero(a[0]), a.order)
    @inbounds for ordT in eachindex(a)
        mul!(aa, pi, a, ordT)
    end
    sincos!(s, c, aa, k)
    return nothing
end

@inline function tan!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
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

@inline function asin!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
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
    # Compute k-th coefficient of auxiliary term s=1-a^2
    @inbounds zero!(r, k) # r[k] <- 0
    @inbounds sqr!(r, a, k) # r[k] <- (a^2)[k]
    @inbounds subst!(r, r, k) # r[k] <- -r[k]
    @inbounds sqrt!(r, r, k) # r[k] <- (sqrt(r))[k]
    return nothing
end

@inline function acos!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
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
    # Compute k-th coefficient of auxiliary term s=1-a^2
    @inbounds zero!(r, k) # r[k] <- 0
    @inbounds sqr!(r, a, k) # r[k] <- (a^2)[k]
    @inbounds subst!(r, r, k) # r[k] <- -r[k]
    @inbounds sqrt!(r, r, k) # r[k] <- (sqrt(r))[k]
    return nothing
end

@inline function atan!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        r::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
    if k == 0
        res[0] = atan( a[0] )
        # zero!(r, a, 0)
        sqr!(r, a, 0)
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
    tmp = TaylorN( zero(a[0][0][1]), a[0].order )
    @inbounds for ordQ in eachindex(a[0])
        # zero!(tmp, res[k], ordQ)
        tmp[ordQ] = - res[k][ordQ] / k
        add!(tmp, a[k], tmp, ordQ)
        zero!(res[k][ordQ])
        div!(res[k], tmp, r[0], ordQ)
    end
    zero!(r[k])
    sqr!(r, a, k)
    return nothing
end

@inline function sinhcosh!(s::Taylor1{TaylorN{T}}, c::Taylor1{TaylorN{T}},
        a::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
    if k == 0
        @inbounds for ordQ in eachindex(a[0])
            sinhcosh!(s[0], c[0], a[0], ordQ)
        end
        return nothing
    end
    # The recursion formula
    zero!(s[k])
    zero!(c[k])
    @inbounds for i = 1:k
        for ordQ in eachindex(a[0])
            mul_scalar!(s[k], i, a[i], c[k-i], ordQ)
            mul_scalar!(c[k], i, a[i], s[k-i], ordQ)
        end
    end
    div!(s, s, k, k)
    div!(c, c, k, k)
    return nothing
end

@inline function tanh!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        res2::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
    if k == 0
        @inbounds res[0] = tanh( a[0] )
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
    tmp = TaylorN( zero(a[0][0][1]), a[0].order)
    @inbounds for ordQ in eachindex(a[0])
        # zero!(tmp, res[k], ordQ)
        tmp[ordQ] = res[k][ordQ] / k
        subst!(res[k], a[k], tmp, ordQ)
    end
    zero!(res2[k])
    sqr!(res2, res, k)
    return nothing
end

@inline function asinh!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        r::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
    if k == 0
        @inbounds res[0] = asinh( a[0] )
        # r[0] = sqrt(1+a[0]^2)
        tmp = TaylorN( zero(a[0][0][1]), a[0].order)
        r[0] = square(a[0])
        for ordQ in eachindex(a[0])
            one!(tmp, a[0], ordQ)
            add!(tmp, tmp, r[0], ordQ)
            # zero!(r[0], tmp, ordQ)
            sqrt!(r[0], tmp, ordQ)
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
    tmp = TaylorN( zero(a[0][0][1]), a[0].order)
    @inbounds for ordQ in eachindex(a[0])
        subst!(tmp, a[k], res[k], ordQ)
        zero!(res[k][ordQ])
        div!(res[k], tmp, r[0], ordQ)
    end
    # Compute auxiliary term s=1+a^2
    s = Taylor1(zero(a[0]), a.order)
    for i = 0:k
        sqr!(s, a, i)
        if i == 0
            s[0] = one(s[0]) + s[0]
            add!(s, one(s), s, 0)
        end
    end
    # Update aux term r = sqrt(s) = sqrt(1+a^2)
    sqrt!(r, s, k)
    return nothing
end

@inline function acosh!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        r::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
    if k == 0
        @inbounds res[0] = acosh( a[0] )
        # r[0] = sqrt(a[0]^2-1)
        tmp = TaylorN( zero(a[0][0][1]), a[0].order)
        r[0] = square(a[0])
        for ordQ in eachindex(a[0])
            one!(tmp, a[0], ordQ)
            subst!(tmp, r[0], tmp, ordQ)
            # zero!(r[0], tmp, ordQ)
            sqrt!(r[0], tmp, ordQ)
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
    tmp = TaylorN( zero(a[0][0][1]), a[0].order)
    @inbounds for ordQ in eachindex(a[0])
        subst!(tmp, a[k], res[k], ordQ)
        zero!(res[k][ordQ])
        div!(res[k], tmp, r[0], ordQ)
    end
    # Compute auxiliary term s=a^2-1
    s = Taylor1(zero(a[0]), a.order)
    for i = 0:k
        sqr!(s, a, i)
        if i == 0
            s[0] = one(s[0]) + s[0]
            subst!(s, s, one(s), 0)
        end
    end
    # Update aux term r = sqrt(s) = sqrt(a^2-1)
    sqrt!(r, s, k)
    return nothing
end

@inline function atanh!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        r::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
    if k == 0
        res[0] = atanh( a[0] )
        # zero!(r, a, 0)
        sqr!(r, a, 0)
        subst!(r, one(a[0][0][1]), r, 0)
        return nothing
    end
    # The recursion formula
    tmp = TaylorN( zero(a[0][0][1]), a[0].order )
    zero!(res[k])
    for i in 1:k-1
        @inbounds for ordQ in eachindex(a[0])
            mul_scalar!(res[k], k-i, res[k-i], r[i], ordQ)
        end
    end
    @inbounds for ordQ in eachindex(a[0])
        # zero!(tmp, res[k], ordQ)
        tmp[ordQ] = res[k][ordQ] / k
        add!(tmp, a[k], tmp, ordQ)
        zero!(res[k][ordQ])
        div!(res[k], tmp, r[0], ordQ)
    end
    zero!(r[k])
    sqr!(r, a, k)
    return nothing
end




@doc doc"""
    inverse(f)

Return the Taylor expansion of ``f^{-1}(t)``, of order `N = f.order`,
for `f::Taylor1` polynomial, assuming the first coefficient of `f` is zero.
Otherwise, a `DomainError` is thrown.

The algorithm implements Lagrange inversion at ``t=0`` if ``f(0)=0``:
```math
\begin{equation*}
f^{-1}(t) = \sum_{n=1}^{N} \frac{t^n}{n!} \left.
    \frac{{\rm d}^{n-1}}{{\rm d} z^{n-1}}\left(\frac{z}{f(z)}\right)^n
    \right\vert_{z=0}.
\end{equation*}
```

""" inverse

function inverse(f::Taylor1{T}) where {T<:Number}
    if !iszero(f[0])
        throw(DomainError(f,
        """
        Evaluation of Taylor1 series at 0 is non-zero; revert
        a Taylor1 series with constant coefficient 0 and re-expand about f(0).
        """))
    end
    z = Taylor1(T, f.order)
    zdivf = z/f
    zdivfpown = zdivf
    res = Taylor1(zero(TS.numtype(zdivf)), f.order)

    @inbounds for ord in 1:f.order
        res[ord] = zdivfpown[ord-1]/ord
        zdivfpown *= zdivf
    end
    return res
end


@doc doc"""
    inverse_map(f)

Return the Taylor expansion of ``f^{-1}(t)``, of order `N = f.order`,
for `Taylor1` or `TaylorN` polynomials, assuming the first coefficient of `f` is zero.
Otherwise, a `DomainError` is thrown.

This method is based in the algorithm by M. Berz, Modern map methods in
Particle Beam Physics, Academic Press (1999), Sect 2.3.1.
See [`inverse`](@ref) (for `f::Taylor1`).

"""
function inverse_map(p::Taylor1)
    if !iszero(constant_term(p))
        throw(DomainError(p,
        """
        Evaluation of Taylor1 series at 0 is non-zero; revert
        a Taylor1 series with constant coefficient 0 and re-expand about f(0).
        """))
    end
    inv_m_pol = inv(linear_polynomial(p)[1])
    n_pol = inv_m_pol * nonlinear_polynomial(p)
    scaled_ident = inv_m_pol * Taylor1(p.order)
    res = scaled_ident
    aux1 = zero(res)
    aux2 = zero(res)
    for ord in 1:p.order
        _horner!(aux2, n_pol, res, aux1)
        subst!(res, scaled_ident, aux2, ord)
    end
    return res
end

function inverse_map(p::Vector{TaylorN{T}}) where {T<:NumberNotSeries}
    if !iszero(constant_term(p))
        throw(DomainError(p,
        """
        Evaluation of Taylor1 series at 0 is non-zero; revert
        a Taylor1 series with constant coefficient 0 and re-expand about f(0).
        """))
    end
    @assert length(p) == get_numvars()
    inv_m_pol = inv(jacobian(p))
    n_pol = inv_m_pol * nonlinear_polynomial(p)
    scaled_ident = inv_m_pol * TaylorN.(1:get_numvars(), order=get_order(p[1]))
    res = deepcopy(scaled_ident)
    aux = zero.(res)
    auxvec = [zero(res[1]) for val in eachindex(res[1])]
    valscache = [zero(val) for val in res]
    aaux = zero(res[1])
    for ord in 1:get_order(p[1])
        t_res = (res...,)
        for i = 1:get_numvars()
            zero!.(auxvec)
            _evaluate!(auxvec, n_pol[i], t_res, valscache, aaux)
            aux[i] = sum( auxvec )
            subst!(res[i], scaled_ident[i], aux[i], ord)
        end
    end
    return res
end


# Documentation for the recursion relations
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
