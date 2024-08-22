# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

# Arithmetic operations: +, -, *, /

## Equality ##
for T in (:Taylor1, :TaylorN)
    @eval begin
        ==(a::$T{T}, b::$T{S}) where {T<:Number,S<:Number} = ==(promote(a,b)...)

        function ==(a::$T{T}, b::$T{T}) where {T<:Number}
            if a.order != b.order
                a, b = fixorder(a, b)
            end
            return a.coeffs == b.coeffs
        end
    end
end

function ==(a::Taylor1{TaylorN{T}}, b::TaylorN{Taylor1{S}}) where {T, S}
    R = promote_type(T, S)
    return a == convert(Taylor1{TaylorN{R}}, b)
end
==(b::TaylorN{Taylor1{S}}, a::Taylor1{TaylorN{T}}) where {T, S} = a == b

function ==(a::HomogeneousPolynomial, b::HomogeneousPolynomial)
    a.order == b.order && return a.coeffs == b.coeffs
    return iszero(a.coeffs) && iszero(b.coeffs)
end


## Total ordering ##
for T in (:Taylor1, :TaylorN)
    @eval begin
        @inline function isless(a::$T{<:Number}, b::Real)
            a0 = constant_term(a)
            a0 != b && return isless(a0, b)
            nz = findfirst(a-b)
            if nz == -1
                return isless(zero(a0), zero(b))
            else
                return isless(a[nz], zero(b))
            end
        end
        @inline function isless(b::Real, a::$T{<:Number})
            a0 = constant_term(a)
            a0 != b && return isless(b, a0)
            nz = findfirst(b-a)
            if nz == -1
                return isless(zero(b), zero(a0))
            else
                return isless(zero(b), a[nz])
            end
        end
        #
        @inline isless(a::$T{T}, b::$T{S}) where {T<:Number, S<:Number} =
            isless(promote(a,b)...)
        @inline isless(a::$T{T}, b::$T{T}) where {T<:Number} =
            isless(a - b, zero(constant_term(a)))
    end
end

@inline function isless(a::HomogeneousPolynomial{<:Number}, b::Real)
    orda = get_order(a)
    if orda == 0
        return isless(a[1], b)
    else
        !iszero(b) && return isless(zero(a[1]), b)
        nz = max(findfirst(a), 1)
        return isless(a[nz], b)
    end
end
@inline function isless(b::Real, a::HomogeneousPolynomial{<:Number})
    orda = get_order(a)
    if orda == 0
        return isless(b, a[1])
    else
        !iszero(b) && return isless(b, zero(a[1]))
        nz = max(findfirst(a),1)
        return isless(b, a[nz])
    end
end
#
@inline isless(a::HomogeneousPolynomial{T}, b::HomogeneousPolynomial{S}) where
    {T<:Number, S<:Number} = isless(promote(a,b)...)
@inline function isless(a::HomogeneousPolynomial{T},
        b::HomogeneousPolynomial{T}) where {T<:Number}
    orda = get_order(a)
    ordb = get_order(b)
    if orda == ordb
        return isless(a-b, zero(a[1]))
    elseif orda < ordb
        return isless(a, zero(a[1]))
    else
        return isless(-b, zero(a[1]))
    end
end

# Mixtures
@inline isless(a::Taylor1{TaylorN{T}}, b::Taylor1{TaylorN{T}}) where
    {T<:NumberNotSeries} = isless(a - b, zero(T))
@inline function isless(a::HomogeneousPolynomial{Taylor1{T}},
        b::HomogeneousPolynomial{Taylor1{T}}) where {T<:NumberNotSeries}
    orda = get_order(a)
    ordb = get_order(b)
    if orda == ordb
        return isless(a-b, zero(T))
    elseif orda < ordb
        return isless(a, zero(T))
    else
        return isless(-b, zero(T))
    end
end
@inline isless(a::TaylorN{Taylor1{T}}, b::TaylorN{Taylor1{T}}) where
    {T<:NumberNotSeries} = isless(a - b, zero(T))

#= TODO: Nested Taylor1s; needs careful thinking; iss #326. The following works:
@inline isless(a::Taylor1{Taylor1{T}}, b::Taylor1{Taylor1{T}}) where {T<:Number} = isless(a - b, zero(T))
# Is the following correct?
# ti = Taylor1(3)
# to = Taylor1([zero(ti), one(ti)], 9)
# tito = ti * to
# ti > to > 0 # ok
# to^2 < toti < ti^2 # ok
# ti > ti^2 > to # is this ok?
=#

@doc doc"""
    isless(a::Taylor1{<:Real}, b::Real)
    isless(a::TaylorN{<:Real}, b::Real)

Compute `isless` by comparing the `constant_term(a)` and `b`. If they are equal,
returns `a[nz] < 0`, with `nz` the first
non-zero coefficient after the constant term. This defines a total order.

For many variables, the ordering includes a lexicographical convention in order to be
total. We have opted for the simplest one, where the *larger* variable appears *before*
when the `TaylorN` variables are defined (e.g., through [`set_variables`](@ref)).

Refs:
- M. Berz, AIP Conference Proceedings 177, 275 (1988); https://doi.org/10.1063/1.37800
- M. Berz, "Automatic Differentiation as Nonarchimedean Analysis", Computer Arithmetic and
    Enclosure Methods, (1992), Elsevier, 439-450.

---

    isless(a::Taylor1{<:Real}, b::Taylor1{<:Real})
    isless(a::TaylorN{<:Real}, b::Taylor1{<:Real})

Returns `isless(a - b, zero(b))`.
""" isless


## zero and one ##
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval iszero(a::$T) = iszero(a.coeffs)
end

for T in (:Taylor1, :TaylorN)
    @eval zero(a::$T) = $T(zero.(a.coeffs))
    @eval function one(a::$T)
        b = zero(a)
        b[0] = one(b[0])
        return b
    end
end

zero(a::HomogeneousPolynomial{T}) where {T<:Number} =
    HomogeneousPolynomial(zero.(a.coeffs), a.order)

function zeros(a::HomogeneousPolynomial{T}, order::Int) where {T<:Number}
    order == 0 && return [HomogeneousPolynomial([zero(a[1])], 0)]
    v = Array{HomogeneousPolynomial{T}}(undef, order+1)
    @simd for ord in eachindex(v)
        @inbounds v[ord] = HomogeneousPolynomial(zero(a[1]), ord-1)
    end
    return v
end

zeros(::Type{HomogeneousPolynomial{T}}, order::Int) where {T<:Number} =
    zeros( HomogeneousPolynomial([zero(T)], 0), order)

function one(a::HomogeneousPolynomial{T}) where {T<:Number}
    v = one.(a.coeffs)
    return HomogeneousPolynomial(v, a.order)
end

function ones(a::HomogeneousPolynomial{T}, order::Int) where {T<:Number}
    order == 0 && return [HomogeneousPolynomial([one(a[1])], 0)]
    v = Array{HomogeneousPolynomial{T}}(undef, order+1)
    @simd for ord in eachindex(v)
        @inbounds num_coeffs = size_table[ord]
        @inbounds v[ord] = HomogeneousPolynomial(ones(T, num_coeffs), ord-1)
    end
    return v
end

ones(::Type{HomogeneousPolynomial{T}}, order::Int) where {T<:Number} =
    ones( HomogeneousPolynomial([one(T)], 0), order)



## Addition and subtraction ##
for (f, fc) in ((:+, :(add!)), (:-, :(subst!)))

    for T in (:Taylor1, :TaylorN)
        @eval begin
            ($f)(a::$T{T}, b::$T{S}) where {T<:Number, S<:Number} =
                ($f)(promote(a, b)...)

            function ($f)(a::$T{T}, b::$T{T}) where {T<:Number}
                if a.order != b.order
                    a, b = fixorder(a, b)
                end
                c = $T( zero(constant_term(a)), a.order)
                for k in eachindex(a)
                    ($fc)(c, a, b, k)
                end
                return c
            end

            function ($f)(a::$T)
                c = $T( zero(constant_term(a)), a.order)
                for k in eachindex(a)
                    ($fc)(c, a, k)
                end
                return c
            end

            ($f)(a::$T{T}, b::S) where {T<:Number, S<:Number} = ($f)(promote(a, b)...)

            function ($f)(a::$T{T}, b::T) where {T<:Number}
                coeffs = copy(a.coeffs)
                @inbounds coeffs[1] = $f(a[0], b)
                return $T(coeffs, a.order)
           end

            # ($f)(b::S, a::$T{T}) where {T<:Number,S<:Number} = $f(promote(b, a)...)

            function ($f)(b::T, a::$T{T}) where {T<:Number}
                coeffs = similar(a.coeffs)
                @__dot__ coeffs = ($f)(a.coeffs)
                @inbounds coeffs[1] = $f(b, a[0])
                return $T(coeffs, a.order)
            end

            ## add! and subst! ##
            function ($fc)(v::$T{T}, a::$T{T}, k::Int) where {T<:Number}
                if $T == Taylor1
                    @inbounds v[k] = ($f)(a[k])
                else
                    @inbounds for l in eachindex(v[k])
                        v[k][l] = ($f)(a[k][l])
                    end
                end
                return nothing
            end

            function ($fc)(v::$T{T}, a::T, k::Int) where {T<:Number}
                @inbounds v[k] = k==0 ? ($f)(a) : zero(a)
                return nothing
            end

            if $T == Taylor1
                function ($fc)(v::$T, a::$T, b::$T, k::Int)
                    @inbounds v[k] = ($f)(a[k], b[k])
                    return nothing
                end
            else
                function ($fc)(v::$T, a::$T, b::$T, k::Int)
                    @inbounds for i in eachindex(v[k])
                        v[k][i] = ($f)(a[k][i], b[k][i])
                    end
                    return nothing
                end
            end

            function ($fc)(v::$T, a::$T, b::Number, k::Int)
                @inbounds v[k] = k==0 ?
                    ($f)(constant_term(a), b) : ($f)(a[k], zero(b))
                return nothing
            end

            function ($fc)(v::$T, a::Number, b::$T, k::Int)
                @inbounds v[k] = k==0 ?
                    ($f)(a, constant_term(b)) : ($f)(zero(a), b[k])
                return nothing
            end
        end

    end

    @eval ($f)(a::T, b::S) where {T<:Taylor1, S<:TaylorN} = ($f)(promote(a, b)...)
    @eval ($f)(a::T, b::S) where {T<:TaylorN, S<:Taylor1} = ($f)(promote(a, b)...)

    @eval begin
        ($f)(a::HomogeneousPolynomial{T}, b::HomogeneousPolynomial{S}) where
            {T<:NumberNotSeriesN,S<:NumberNotSeriesN} = ($f)(promote(a,b)...)

        function ($f)(a::HomogeneousPolynomial{T}, b::HomogeneousPolynomial{T}) where
                {T<:NumberNotSeriesN}
            @assert a.order == b.order
            v = similar(a.coeffs)
            @__dot__ v = ($f)(a.coeffs, b.coeffs)
            return HomogeneousPolynomial(v, a.order)
        end

        # NOTE add! and subst! for HomogeneousPolynomial's act as += or -=
        function ($fc)(res::HomogeneousPolynomial{T}, a::HomogeneousPolynomial{T},
                b::HomogeneousPolynomial{T}, k::Int) where {T<:NumberNotSeriesN}
            res[k] += ($f)(a[k], b[k])
            return nothing
        end

        function ($f)(a::HomogeneousPolynomial)
            v = similar(a.coeffs)
            @__dot__ v = ($f)(a.coeffs)
            return HomogeneousPolynomial(v, a.order)
        end

        function ($f)(a::TaylorN{Taylor1{T}}, b::Taylor1{S}) where
                {T<:NumberNotSeries, S<:NumberNotSeries}
            @inbounds aux = $f(a[0][1], b)
            R = TS.numtype(aux)
            coeffs = Array{HomogeneousPolynomial{Taylor1{R}}}(undef, a.order+1)
            coeffs .= a.coeffs
            @inbounds coeffs[1] = aux
            return TaylorN(coeffs, a.order)
        end

        function ($f)(b::Taylor1{S}, a::TaylorN{Taylor1{T}}) where
                {T<:NumberNotSeries, S<:NumberNotSeries}
            @inbounds aux = $f(b, a[0][1])
            R = TS.numtype(aux)
            coeffs = Array{HomogeneousPolynomial{Taylor1{R}}}(undef, a.order+1)
            @__dot__ coeffs = $f(a.coeffs)
            @inbounds coeffs[1] = aux
            return TaylorN(coeffs, a.order)
        end

        function ($f)(a::Taylor1{TaylorN{T}}, b::TaylorN{S}) where
                {T<:NumberNotSeries,S<:NumberNotSeries}
            @inbounds aux = $f(a[0], b)
            c = Taylor1( zero(aux), a.order)
            for k in eachindex(a)
                ($fc)(c, a, b, k)
            end
            return c
    end

        function ($f)(b::TaylorN{S}, a::Taylor1{TaylorN{T}}) where
                {T<:NumberNotSeries,S<:NumberNotSeries}
            @inbounds aux = $f(b, a[0])
            c = Taylor1( zero(aux), a.order)
            for k in eachindex(a)
                ($fc)(c, a, b, k)
            end
            return c
        end

    end
end

for (f, fc) in ((:+, :(add!)), (:-, :(subst!)))
    @eval begin
        function ($f)(a::Taylor1{TaylorN{T}}, b::Taylor1{TaylorN{T}}) where {T<:NumberNotSeries}
            if a.order != b.order || any(get_order.(a.coeffs) .!= get_order.(b.coeffs))
                a, b = fixorder(a, b)
            end
            c = zero(a)
            for k in eachindex(a)
                ($fc)(c, a, b, k)
            end
            return c
        end
        function ($fc)(v::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}}, b::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
            @inbounds for i in eachindex(v[k])
                for j in eachindex(v[k][i])
                    v[k][i][j] = ($f)(a[k][i][j], b[k][i][j])
                end
            end
            return nothing
        end
        function ($fc)(v::Taylor1{TaylorN{T}}, a::NumberNotSeries, b::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
            @inbounds for i in eachindex(v[k])
                for j in eachindex(v[k][i])
                    v[k][i][j] = ($f)(k==0 && i==0 && j==1 ? a : zero(a), b[k][i][j])
                end
            end
            return nothing
        end
        function ($fc)(v::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
            @inbounds for l in eachindex(v[k])
                for m in eachindex(v[k][l])
                    v[k][l][m] = ($f)(a[k][l][m])
                end
            end
            return nothing
        end
    end
end

for T in (:Taylor1, :TaylorN)
    @eval begin
        function sum!(v::$T{S}, a::AbstractArray{$T{S}}) where {S <: Number}
            for i in eachindex(a)
                for k in eachindex(v)
                    add!(v, v, a[i], k)
                end
            end
            return nothing
        end
    end
end

function sum!(v::TaylorN{S}, a::AbstractArray{HomogeneousPolynomial{S}}) where {S <: Number}
    for i in eachindex(a)
        for k in eachindex(v)
            add!(v, v, a[i], k)
        end
    end
    return nothing
end

## Multiplication ##
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)

    @eval begin
        function *(a::T, b::$T{S}) where {T<:NumberNotSeries, S<:NumberNotSeries}
            @inbounds aux = a * b.coeffs[1]
            v = Array{typeof(aux)}(undef, length(b.coeffs))
            @__dot__ v = a * b.coeffs
            return $T(v, b.order)
        end

        *(b::$T{S}, a::T) where {T<:NumberNotSeries, S<:NumberNotSeries} = a * b

        function *(a::T, b::$T{T}) where {T<:Number}
            v = Array{T}(undef, length(b.coeffs))
            @__dot__ v = a * b.coeffs
            return $T(v, b.order)
        end

        *(b::$T{T}, a::T) where {T<:Number} = a * b
    end
end

for T in (:HomogeneousPolynomial, :TaylorN)

    @eval begin
        function *(a::Taylor1{T}, b::$T{Taylor1{S}}) where
                {T<:NumberNotSeries, S<:NumberNotSeries}
            @inbounds aux = a * b.coeffs[1]
            R = typeof(aux)
            coeffs = Array{R}(undef, length(b.coeffs))
            @__dot__ coeffs = a * b.coeffs
            return $T(coeffs, b.order)
        end

        *(b::$T{Taylor1{R}}, a::Taylor1{T}) where
            {T<:NumberNotSeries, R<:NumberNotSeries} = a * b

        function *(a::$T{T}, b::Taylor1{$T{S}}) where {T<:NumberNotSeries, S<:NumberNotSeries}
            @inbounds aux = a * b[0]
            R = typeof(aux)
            coeffs = Array{R}(undef, length(b.coeffs))
            @__dot__ coeffs = a * b.coeffs
            return Taylor1(coeffs, b.order)
        end

        *(b::Taylor1{$T{S}}, a::$T{T}) where {T<:NumberNotSeries, S<:NumberNotSeries} = a * b
    end
end

for (T, W) in ((:Taylor1, :Number), (:TaylorN, :NumberNotSeriesN))
    @eval function *(a::$T{T}, b::$T{T}) where {T<:$W}
        if a.order != b.order
            a, b = fixorder(a, b)
        end
        c = $T(zero(constant_term(a)), a.order)
        for ord in eachindex(c)
            mul!(c, a, b, ord) # updates c[ord]
        end
        return c
    end
end


*(a::HomogeneousPolynomial{T}, b::HomogeneousPolynomial{S}) where
    {T<:NumberNotSeriesN,S<:NumberNotSeriesN} = *(promote(a,b)...)

function *(a::HomogeneousPolynomial{T}, b::HomogeneousPolynomial{T}) where
        {T<:NumberNotSeriesN}

    order = a.order + b.order

    # NOTE: the following returns order 0, but could be get_order(), or get_order(a)
    order > get_order() && return HomogeneousPolynomial(zero(a[1]), get_order(a))

    res = HomogeneousPolynomial(zero(a[1]), order)
    mul!(res, a, b)
    return res
end

function *(a::Taylor1{TaylorN{T}}, b::Taylor1{TaylorN{S}}) where
        {T<:NumberNotSeries, S<:NumberNotSeries}
    R = promote_type(T,S)
    return *(convert(Taylor1{TaylorN{R}}, a), convert(Taylor1{TaylorN{R}}, b))
end

function *(a::Taylor1{TaylorN{T}}, b::Taylor1{TaylorN{T}}) where {T<:NumberNotSeries}
    if (a.order != b.order) || any(get_order.(a.coeffs) .!= get_order.(b.coeffs))
        a, b = fixorder(a, b)
    end
    res = Taylor1(zero(a[0]), a.order)
    for ordT in eachindex(a)
        mul!(res, a, b, ordT)
    end
    return res
end


# Internal multiplication functions
for T in (:Taylor1, :TaylorN)
    # NOTE: For $T = TaylorN, `mul!` *accumulates* the result of a * b in c[k]
    @eval @inline function mul!(c::$T{T}, a::$T{T}, b::$T{T}, k::Int) where {T<:Number}
        if $T == Taylor1
            @inbounds c[k] = a[0] * b[k]
            @inbounds for i = 1:k
                c[k] += a[i] * b[k-i]
            end
        else
            @inbounds mul!(c[k], a[0], b[k])
            @inbounds for i = 1:k
                mul!(c[k], a[i], b[k-i])
            end
        end
        return nothing
    end

    @eval @inline function mul_scalar!(c::$T{T}, scalar::NumberNotSeries, a::$T{T}, b::$T{T}, k::Int) where {T<:Number}
        if $T == Taylor1
            @inbounds c[k] = scalar * a[0] * b[k]
            @inbounds for i = 1:k
                c[k] += scalar * a[i] * b[k-i]
            end
        else
            @inbounds mul_scalar!(c[k], scalar, a[0], b[k])
            @inbounds for i = 1:k
                mul_scalar!(c[k], scalar, a[i], b[k-i])
            end
        end
        return nothing
    end

    @eval begin
        if $T == Taylor1
            @inline function mul!(v::$T, a::$T, b::NumberNotSeries, k::Int)
                @inbounds v[k] = a[k] * b
                return nothing
            end
            @inline function mul!(v::$T, a::NumberNotSeries, b::$T, k::Int)
                @inbounds v[k] = a * b[k]
                return nothing
            end
            @inline function muladd!(v::$T, a::$T, b::NumberNotSeries, k::Int)
                @inbounds v[k] += a[k] * b
                return nothing
            end
            @inline function muladd!(v::$T, a::NumberNotSeries, b::$T, k::Int)
                @inbounds v[k] += a * b[k]
                return nothing
            end
        else
            @inline function mul!(v::$T, a::$T, b::NumberNotSeries, k::Int)
                @inbounds for i in eachindex(v[k])
                    v[k][i] = a[k][i] * b
                end
                return nothing
            end
            @inline function mul!(v::$T, a::NumberNotSeries, b::$T, k::Int)
                @inbounds for i in eachindex(v[k])
                    v[k][i] = a * b[k][i]
                end
                return nothing
            end
            @inline function muladd!(v::$T, a::$T, b::NumberNotSeries, k::Int)
                @inbounds for i in eachindex(v[k])
                    v[k][i] += a[k][i] * b
                end
                return nothing
            end
            @inline function muladd!(v::$T, a::NumberNotSeries, b::$T, k::Int)
                @inbounds for i in eachindex(v[k])
                    v[k][i] += a * b[k][i]
                end
                return nothing
            end
        end
    end

    @eval @inline function mul!(v::$T, a::$T, b::NumberNotSeries)
        for k in eachindex(v)
            mul!(v, a, b, k)
        end
        return nothing
    end
    @eval @inline function mul!(v::$T, a::NumberNotSeries, b::$T)
        for k in eachindex(v)
            mul!(v, a, b, k)
        end
        return nothing
    end
end

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
function mul!(a::Taylor1{T}, b::Taylor1{T}) where {T<:Number}
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
function mul!(a::Taylor1{TaylorN{T}}, b::Taylor1{TaylorN{T}}) where {T<:NumberNotSeries}
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

function mul!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}}, b::Taylor1{TaylorN{T}},
        ordT::Int) where {T<:NumberNotSeries}
    # Sanity
    zero!(res, ordT)
    for k in 0:ordT
        @inbounds for ordQ in eachindex(a[ordT])
            mul!(res[ordT], a[k], b[ordT-k], ordQ)
        end
    end
    return nothing
end

@inline function mul!(res::Taylor1{TaylorN{T}}, a::NumberNotSeries,
    b::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
for l in eachindex(b[k])
    for m in eachindex(b[k][l])
        res[k][l][m] = a*b[k][l][m]
    end
end
return nothing
end

mul!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
b::NumberNotSeries, k::Int) where {T<:NumberNotSeries} = mul!(res, b, a, k)

# in-place product (assumes equal order among TaylorNs)
# NOTE: the result of the product is *accumulated* in c[k]
function mul!(c::TaylorN, a::TaylorN, b::TaylorN)
    for k in eachindex(c)
        mul!(c, a, b, k)
    end
end

function mul_scalar!(c::TaylorN, scalar::NumberNotSeries, a::TaylorN, b::TaylorN)
    for k in eachindex(c)
        mul_scalar!(c, scalar, a, b, k)
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


"""
    mul!(c, a, b) --> nothing

Accumulates in `c` the result of `a*b` with minimum allocation. Arguments
c, a and b are `HomogeneousPolynomial`.

"""
@inline function mul!(c::HomogeneousPolynomial, a::HomogeneousPolynomial,
        b::HomogeneousPolynomial)

    (iszero(b) || iszero(a)) && return nothing

    @inbounds num_coeffs_a = size_table[a.order+1]
    @inbounds num_coeffs_b = size_table[b.order+1]

    @inbounds posTb = pos_table[c.order+1]

    @inbounds indTa = index_table[a.order+1]
    @inbounds indTb = index_table[b.order+1]

    @inbounds for na in 1:num_coeffs_a
        ca = a[na]
        # iszero(ca) && continue
        inda = indTa[na]

        @inbounds for nb in 1:num_coeffs_b
            cb = b[nb]
            # iszero(cb) && continue
            indb = indTb[nb]

            pos = posTb[inda + indb]
            c[pos] += ca * cb
        end
    end

    return nothing
end


"""
    mul_scalar!(c, scalar, a, b) --> nothing

Accumulates in `c` the result of `scalar*a*b` with minimum allocation. Arguments
c, a and b are `HomogeneousPolynomial`; `scalar` is a NumberNotSeries.

"""
@inline function mul_scalar!(c::HomogeneousPolynomial, scalar::NumberNotSeries, a::HomogeneousPolynomial,
        b::HomogeneousPolynomial)

    (iszero(b) || iszero(a)) && return nothing

    @inbounds num_coeffs_a = size_table[a.order+1]
    @inbounds num_coeffs_b = size_table[b.order+1]

    @inbounds posTb = pos_table[c.order+1]

    @inbounds indTa = index_table[a.order+1]
    @inbounds indTb = index_table[b.order+1]

    @inbounds for na in 1:num_coeffs_a
        ca = a[na]
        # iszero(ca) && continue
        inda = indTa[na]

        @inbounds for nb in 1:num_coeffs_b
            cb = b[nb]
            # iszero(cb) && continue
            indb = indTb[nb]

            pos = posTb[inda + indb]
            c[pos] += scalar * ca * cb
        end
    end

    return nothing
end



## Division ##
function /(a::Taylor1{Rational{T}}, b::S) where {T<:Integer, S<:NumberNotSeries}
    R = typeof( a[0] // b)
    v = Array{R}(undef, a.order+1)
    @__dot__ v = a.coeffs // b
    return Taylor1(v, a.order)
end

for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval function /(a::$T{T}, b::S) where {T<:NumberNotSeries, S<:NumberNotSeries}
        @inbounds aux = a.coeffs[1] / b
        v = Array{typeof(aux)}(undef, length(a.coeffs))
        @__dot__ v = a.coeffs / b
        return $T(v, a.order)
    end

    @eval function /(a::$T{T}, b::T) where {T<:Number}
        @inbounds aux = a.coeffs[1] / b
        # v = Array{typeof(aux)}(undef, length(a.coeffs))
        # @__dot__ v = a.coeffs / b
        # return $T(v, a.order)
        c = $T( zero(aux), a.order )
        for ord in eachindex(c)
            div!(c, a, b, ord) # updates c[ord]
        end
        return c
    end
end

for T in (:HomogeneousPolynomial, :TaylorN)
    @eval function /(b::$T{Taylor1{S}}, a::Taylor1{T}) where
            {T<:NumberNotSeries, S<:NumberNotSeries}
        @inbounds aux = b.coeffs[1] / a
        R = typeof(aux)
        coeffs = Array{R}(undef, length(b.coeffs))
        @__dot__ coeffs = b.coeffs / a
        return $T(coeffs, b.order)
    end

    @eval function /(b::$T{Taylor1{T}}, a::S) where {T<:NumberNotSeries, S<:NumberNotSeries}
        @inbounds aux = b.coeffs[1] / a
        R = typeof(aux)
        coeffs = Array{R}(undef, length(b.coeffs))
        @__dot__ coeffs = b.coeffs / a
        return $T(coeffs, b.order)
    end

    @eval function /(b::Taylor1{$T{S}}, a::$T{T}) where
            {T<:NumberNotSeries, S<:NumberNotSeries}
        @inbounds aux = b[0] / a
        # R = typeof(aux)
        # coeffs = Array{R}(undef, length(b.coeffs))
        # @__dot__ coeffs = b.coeffs / a
        # return Taylor1(coeffs, b.order)
        v = Taylor1(zero(aux), b.order)
        @inbounds for k in eachindex(b)
            v[k] = b[k] / a
        end
        return v
    end
end

/(a::Taylor1{T}, b::Taylor1{S}) where {T<:Number, S<:Number} = /(promote(a,b)...)

function /(a::Taylor1{T}, b::Taylor1{T}) where {T<:Number}
    iszero(a) && !iszero(b) && return zero(a)
    if a.order != b.order
        a, b = fixorder(a, b)
    end

    # order and coefficient of first factorized term
    ordfact, cdivfact = divfactorization(a, b)

    c = Taylor1(cdivfact, a.order-ordfact)
    for ord in eachindex(c)
        div!(c, a, b, ord) # updates c[ord]
    end

    return c
end

/(a::TaylorN{T}, b::TaylorN{S}) where
    {T<:NumberNotSeriesN, S<:NumberNotSeriesN} = /(promote(a,b)...)

function /(a::TaylorN{T}, b::TaylorN{T}) where {T<:NumberNotSeriesN}
    @assert !iszero(constant_term(b))

    if a.order != b.order
        a, b = fixorder(a, b)
    end

    # first coefficient
    @inbounds cdivfact = a[0] / constant_term(b)
    c = TaylorN(cdivfact, a.order)
    for ord in eachindex(c)
        div!(c, a, b, ord) # updates c[ord]
    end

    return c
end

function /(a::Taylor1{TaylorN{T}}, b::Taylor1{TaylorN{T}}) where {T<:NumberNotSeries}
    iszero(a) && !iszero(b) && return zero(a)
    if (a.order != b.order) || any(get_order.(a.coeffs) .!= get_order.(b.coeffs))
        a, b = fixorder(a, b)
    end

    # order and coefficient of first factorized term
    ordfact, cdivfact = divfactorization(a, b)

    res = Taylor1(cdivfact, a.order-ordfact)
    for ordT in eachindex(res)
        div!(res, a, b, ordT)
    end
    return res
end

function /(a::S, b::Taylor1{TaylorN{T}}) where {S<:NumberNotSeries, T<:NumberNotSeries}
    R = promote_type(TaylorN{S}, TaylorN{T})
    res = convert(Taylor1{R}, zero(b))
    iszero(a) && !iszero(b) && return res

    for ordT in eachindex(res)
        div!(res, a, b, ordT)
    end
    return res
end

function /(a::TaylorN{T}, b::Taylor1{TaylorN{T}}) where {T<:NumberNotSeries}
    res = zero(b)
    iszero(a) && !iszero(b) && return res

    aa = Taylor1(a, b.order)
    for ordT in eachindex(res)
        div!(res, aa, b, ordT)
    end
    return res
end

@inline function divfactorization(a1::Taylor1, b1::Taylor1)
    # order of first factorized term; a1 and b1 assumed to be of the same order
    a1nz = findfirst(a1)
    b1nz = findfirst(b1)
    a1nz = a1nz ≥ 0 ? a1nz : a1.order
    b1nz = b1nz ≥ 0 ? b1nz : a1.order
    ordfact = min(a1nz, b1nz)
    cdivfact = a1[ordfact] / b1[ordfact]

    # Is the polynomial factorizable?
    iszero(b1[ordfact]) && throw( ArgumentError(
        """Division does not define a Taylor1 polynomial;
        order k=$(ordfact) => coeff[$(ordfact)]=$(cdivfact).""") )

    return ordfact, cdivfact
end


## TODO: Implement factorization (divfactorization) for TaylorN polynomials


# Homogeneous coefficient for the division
@doc doc"""
    div!(c, a, b, k::Int)

Compute the `k-th` expansion coefficient `c[k]` of `c = a / b`,
where all `c`, `a` and `b` are either `Taylor1` or `TaylorN`.

The coefficients are given by

```math
c_k =  \frac{1}{b_0} \big(a_k - \sum_{j=0}^{k-1} c_j b_{k-j}\big).
```

For `Taylor1` polynomials, a similar formula is implemented which
exploits `k_0`, the order of the first non-zero coefficient of `a`.
""" div!

@inline function div!(c::Taylor1, a::Taylor1, b::Taylor1, k::Int)

    # order and coefficient of first factorized term
    ordfact, cdivfact = divfactorization(a, b)
    if k == 0
        @inbounds c[0] = cdivfact
        return nothing
    end

    imin = max(0, k+ordfact-b.order)
    @inbounds c[k] = c[imin] * b[k+ordfact-imin]
    @inbounds for i = imin+1:k-1
        c[k] += c[i] * b[k+ordfact-i]
    end
    if k+ordfact ≤ b.order
        @inbounds c[k] = (a[k+ordfact]-c[k]) / b[ordfact]
    else
        @inbounds c[k] = - c[k] / b[ordfact]
    end
    return nothing
end

@inline function div!(v::Taylor1, a::Taylor1, b::NumberNotSeries, k::Int)
    @inbounds v[k] = a[k] / b
    return nothing
end

@inline function div!(c::Taylor1{T}, a::NumberNotSeries,
        b::Taylor1{T}, k::Int) where {T<:Number}
    zero!(c, k)
    iszero(a) && !iszero(b) && return nothing
    # order and coefficient of first factorized term
    # In this case, since a[k]=0 for k>0, we can simplify to:
    # ordfact, cdivfact = 0, a/b[0]
    if k == 0
        @inbounds c[0] = a/b[0]
        return nothing
    end

    @inbounds c[k] = c[0] * b[k]
    @inbounds for i = 1:k-1
        c[k] += c[i] * b[k-i]
    end
    @inbounds c[k] = -c[k]/b[0]
    return nothing
end

@inline function div!(c::Taylor1{TaylorN{T}}, a::NumberNotSeries,
        b::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
    zero!(c, k)
    iszero(a) && !iszero(b) && return nothing
    # order and coefficient of first factorized term
    # In this case, since a[k]=0 for k>0, we can simplify to:
    # ordfact, cdivfact = 0, a/b[0]
    if k == 0
        @inbounds div!(c[0], a, b[0])
        return nothing
    end

    @inbounds mul!(c[k], c[0], b[k])
    @inbounds for i = 1:k-1
        # c[k] += c[i] * b[k-i]
        mul!(c[k], c[i], b[k-i])
    end
    # @inbounds c[k] = -c[k]/b[0]
    @inbounds div_scalar!(c[k], -1, b[0])
    return nothing
end

# TODO: avoid allocations when T isa Taylor1
@inline function div!(v::HomogeneousPolynomial{T}, a::HomogeneousPolynomial{T}, b::NumberNotSeriesN) where {T <: Number}
    @inbounds for k in eachindex(v)
        v[k] = a[k] / b
    end
    return nothing
end

# NOTE: Due to the use of `zero!`, this `div!` method does *not* accumulate the result of a / b in c[k] (k > 0)
@inline function div!(c::TaylorN, a::TaylorN, b::TaylorN, k::Int)
    if k==0
        @inbounds c[0][1] = constant_term(a) / constant_term(b)
        return nothing
    end

    zero!(c, k)

    @inbounds for i = 0:k-1
        mul!(c[k], c[i], b[k-i])
    end
    @inbounds for i in eachindex(c[k])
        c[k][i] = (a[k][i] - c[k][i]) / constant_term(b)
    end
    return nothing
end

# In-place division and assignment: c[k] = (c/a)[k]
# NOTE: Here `div!` *accumulates* the result of (c/a)[k] in c[k] (k > 0)
#
# Recursion algorithm:
#
# k = 0: c[0] <- c[0]/a[0]
# k = 1: c[1] <- c[1] - c[0]*a[1]
#        c[1] <- c[1]/a[0]
# k = 2: c[2] <- c[2] - c[0]*a[2] - c[1]*a[1]
#        c[2] <- c[2]/a[0]
# etc.
@inline function div!(c::TaylorN, a::TaylorN, k::Int)
    if k==0
        @inbounds c[0][1] = constant_term(c) / constant_term(a)
        return nothing
    end

    @inbounds for i = 0:k-1
        mul_scalar!(c[k], -1, c[i], a[k-i])
    end
    @inbounds for i in eachindex(c[k])
        c[k][i] = c[k][i] / constant_term(a)
    end
    return nothing
end

# In-place division and assignment: c[k] <- scalar * (c/a)[k]
# NOTE: Here `div!` *accumulates* the result of scalar * (c/a)[k] in c[k] (k > 0)
#
# Recursion algorithm:
#
# k = 0: c[0] <- scalar*c[0]/a[0]
# k = 1: c[1] <- scalar*c[1] - c[0]*a[1]
#        c[1] <- c[1]/a[0]
# k = 2: c[2] <- scalar*c[2] - c[0]*a[2] - c[1]*a[1]
#        c[2] <- c[2]/a[0]
# etc.
@inline function div_scalar!(c::TaylorN, scalar::NumberNotSeries, a::TaylorN, k::Int)
    if k==0
        @inbounds c[0][1] = scalar*constant_term(c) / constant_term(a)
        return nothing
    end

    @inbounds mul!(c, scalar, c, k)
    @inbounds for i = 0:k-1
        mul_scalar!(c[k], -1, c[i], a[k-i])
    end
    @inbounds for i in eachindex(c[k])
        c[k][i] = c[k][i] / constant_term(a)
    end
    return nothing
end

# NOTE: Here `div!` *accumulates* the result of a[k] / b[k] in c[k] (k > 0)
@inline function div!(c::TaylorN, a::NumberNotSeries, b::TaylorN, k::Int)
    if k==0
        @inbounds c[0][1] = a / constant_term(b)
        return nothing
    end

    @inbounds for i = 0:k-1
        mul!(c[k], c[i], b[k-i])
    end
    @inbounds for i in eachindex(c[k])
        c[k][i] = ( -c[k][i] ) / constant_term(b)
    end
    return nothing
end

# in-place division c <- c/a (assumes equal order among TaylorNs)
function div!(c::TaylorN, a::TaylorN)
    @inbounds for k in eachindex(c)
        div!(c, a, k)
    end
    return nothing
end

# in-place division c <- scalar*c/a (assumes equal order among TaylorNs)
function div_scalar!(c::TaylorN, scalar::NumberNotSeries, a::TaylorN)
    @inbounds for k in eachindex(c)
        div_scalar!(c, scalar, a, k)
    end
    return nothing
end

# c[k] <- (a/b)[k]
function div!(c::TaylorN, a::TaylorN, b::TaylorN)
    @inbounds for k in eachindex(c)
        div!(c, a, b, k)
    end
    return nothing
end

# c[k] <- (a/b)[k], where a is a scalar
function div!(c::TaylorN, a::NumberNotSeries, b::TaylorN)
    @inbounds for k in eachindex(c)
        div!(c, a, b, k)
    end
    return nothing
end

# c[k] <- a[k]/b, where b is a scalar
function div!(c::TaylorN, a::TaylorN, b::NumberNotSeries)
    @inbounds for k in eachindex(c)
        div!(c[k], a[k], b)
    end
    return nothing
end

# NOTE: Here `div!` *accumulates* the result of a / b in res[k] (k > 0)
@inline function div!(c::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        b::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeriesN}

    # order and coefficient of first factorized term
    # ordfact, cdivfact = divfactorization(a, b)
    anz = findfirst(a)
    bnz = findfirst(b)
    anz = anz ≥ 0 ? anz : a.order
    bnz = bnz ≥ 0 ? bnz : a.order
    ordfact = min(anz, bnz)

    # Is the polynomial factorizable?
    iszero(b[ordfact]) && throw( ArgumentError(
        """Division does not define a Taylor1 polynomial;
        order k=$(ordfact) => coeff[$(ordfact)]=$(cdivfact).""") )

    zero!(c, k)

    if k == 0
        # @inbounds c[0] = a[ordfact]/b[ordfact]
        @inbounds div!(c[0], a[ordfact], b[ordfact])
        return nothing
    end

    imin = max(0, k+ordfact-b.order)
    @inbounds mul!(c[k], c[imin], b[k+ordfact-imin])
    @inbounds for i = imin+1:k-1
        mul!(c[k], c[i], b[k+ordfact-i])
    end
        if k+ordfact ≤ b.order
        # @inbounds c[k] = (a[k+ordfact]-c[k]) / b[ordfact]
        @inbounds for l in eachindex(c[k])
            subst!(c[k], a[k+ordfact], c[k], l)
        end
        @inbounds div!(c[k], b[ordfact])
    else
        # @inbounds c[k] = (-c[k]) / b[ordfact]
        @inbounds div_scalar!(c[k], -1, b[ordfact])
    end
    return nothing
end

@inline function div!(res::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
        b::NumberNotSeries, k::Int) where {T<:NumberNotSeries}
    for l in eachindex(a[k])
        for m in eachindex(a[k][l])
            res[k][l][m] = a[k][l][m]/b
        end
    end
    return nothing
end



"""
    mul!(Y, A, B)

Multiply A*B and save the result in Y.
"""
function mul!(y::Vector{Taylor1{T}},
        a::Union{Matrix{T},SparseMatrixCSC{T}},
        b::Vector{Taylor1{T}}) where {T<:Number}

    n, k = size(a)
    @assert (length(y)== n && length(b)== k)

    # determine the maximal order of b
    # order = maximum([b1.order for b1 in b])
    order = maximum(get_order.(b))

    # Use matrices of coefficients (of proper size) and mul!
    # B = zeros(T, k, order+1)
    B = Array{T}(undef, k, order+1)
    B = zero.(B)
    for i = 1:k
        @inbounds ord = b[i].order
        @inbounds for j = 1:ord+1
            B[i,j] = b[i][j-1]
        end
    end
    Y = Array{T}(undef, n, order+1)
    mul!(Y, a, B)
    @inbounds for i = 1:n
        # y[i] = Taylor1( collect(Y[i,:]), order)
        y[i] = Taylor1( Y[i,:], order)
    end

    return y
end


# Adapted from (Julia v1.2) stdlib/v1.2/LinearAlgebra/src/dense.jl#721-734,
# licensed under MIT "Expat".
# Specialize a method of `inv` for Matrix{Taylor1{T}}. Simply, avoid pivoting,
# since the polynomial field is not an ordered one.
# function Base.inv(A::StridedMatrix{Taylor1{T}}) where T
#     checksquare(A)
#     S = Taylor1{typeof((one(T)*zero(T) + one(T)*zero(T))/one(T))}
#     AA = convert(AbstractArray{S}, A)
#     if istriu(AA)
#         Ai = triu!(parent(inv(UpperTriangular(AA))))
#     elseif istril(AA)
#         Ai = tril!(parent(inv(LowerTriangular(AA))))
#     else
#         # Do not use pivoting !!
#         Ai = inv!(lu(AA, Val(false)))
#         Ai = convert(typeof(parent(Ai)), Ai)
#     end
#     return Ai
# end

# see https://github.com/JuliaLang/julia/pull/40623
const LU_RowMaximum = VERSION >= v"1.7.0-DEV.1188" ? RowMaximum() : Val(true)
const LU_NoPivot = VERSION >= v"1.7.0-DEV.1188" ? NoPivot() : Val(false)

# Adapted from (Julia v1.2) stdlib/v1.2/LinearAlgebra/src/lu.jl#240-253
# and (Julia v1.4.0-dev) stdlib/LinearAlgebra/v1.4/src/lu.jl#270-274,
# licensed under MIT "Expat".
# Specialize a method of `lu` for Matrix{Taylor1{T}}, which avoids pivoting,
# since the polynomial field is not an ordered one.
# We can't assume an ordered field so we first try without pivoting
function lu(A::AbstractMatrix{Taylor1{T}}; check::Bool = true) where {T<:Number}
    S = Taylor1{lutype(T)}
    F = lu!(copy_oftype(A, S), LU_NoPivot; check = false)
    if issuccess(F)
        return F
    else
        return lu!(copy_oftype(A, S), LU_RowMaximum; check = check)
    end
end
