# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

# Arithmetic operations: +, -, *, /
## Addition and subtraction ##
for (f, fc) in ((:+, :(add!)), (:-, :(subst!)))

    for T in (:Taylor1, :TaylorN)
        @eval begin
            ($f)(a::$T{T}, b::$T{S}) where {T<:Number, S<:Number} =
                ($f)(promote(a, b)...)

            function ($f)(a::$T{T}, b::$T{T}) where {T<:Number}
                _check_same_space(a, b)
                if order(a) != order(b)
                    a, b = fixorder(a, b)
                end
                c = zero(a)
                for k in eachindex(a)
                    ($fc)(c, a, b, k)
                end
                return c
            end

            function ($f)(a::$T)
                c = zero(a)
                for k in eachindex(a)
                    ($fc)(c, a, k)
                end
                return c
            end

            function ($f)(a::$T{T}, b::T) where {T<:Number}
                c = $T(a.coeffs[:])
                c[0] = $f(a[0], b)
                return c
            end

            function ($f)(b::T, a::$T{T}) where {T<:Number}
                c = $T($f(a.coeffs[:]))
                c[0] = $f(b, a[0])
                return c
            end

            ## add! and subst! ##
            function ($fc)(v::$T{T}, a::T, k::Int) where {T<:Number}
                @inbounds v[k] = k==0 ? ($f)(a) : zero(a)
                return nothing
            end
        end

        if T == :Taylor1
            @eval begin
                function ($f)(a::$T{T}, b::S) where {T<:Number, S<:NumberNotSeries}
                    return ($f)(promote(a, b)...)
                end

                function ($f)(b::S, a::$T{T}) where {T<:Number, S<:NumberNotSeries}
                    return ($f)(promote(b, a)...)
                end

                function ($fc)(v::$T{T}, a::$T{T}, k::Int) where {T<:Number}
                    @inbounds v[k] = ($f)(a[k])
                    return nothing
                end

                function ($fc)(v::$T, a::$T, b::$T, k::Int)
                    @inbounds v[k] = ($f)(a[k], b[k])
                    return nothing
                end

                function ($fc)(v::$T, a::$T, b::Number, k::Int)
                    if k == 0
                        v[k] = ($f)(a[k], b)
                    else
                        v[k] = ($f)(a[k], zero(b))
                    end
                    return nothing
                end

                function ($fc)(v::$T, a::Number, b::$T, k::Int)
                    if k == 0
                        v[k] = ($f)(a, b[k])
                    else
                        v[k] = ($f)(zero(a), b[k])
                    end
                    return nothing
                end

                # Nested Taylor1s
                function ($fc)(v::$T{$T{T}}, a::$T{$T{T}}, k::Int) where
                        {T<:NumberNotSeriesN}
                    @inbounds for i in eachindex(v[k])
                        v[k][i] = ($f)(a[k][i])
                    end
                    return nothing
                end

                function ($fc)(v::$T{$T{T}}, a::$T{$T{T}}, b::$T{$T{T}}, k::Int) where
                        {T<:NumberNotSeriesN}
                    @inbounds for i in eachindex(v[k])
                        ($fc)(v[k], a[k], b[k], i)
                    end
                    return nothing
                end

                function ($fc)(v::$T{$T{T}}, a::$T{$T{T}}, b::$T{T}, k::Int) where
                        {T<:NumberNotSeriesN}
                    @inbounds for i in eachindex(v[k])
                        ($fc)(v[k], a[k], b, i)
                    end
                    return nothing
                end

                function ($fc)(v::$T{$T{T}}, a::$T{T}, b::$T{$T{T}}, k::Int) where
                        {T<:NumberNotSeriesN}
                    @inbounds for i in eachindex(v[k])
                        ($fc)(v[k], a, b[k], i)
                    end
                    return nothing
                end

                function ($fc)(v::$T{$T{T}}, a::$T{$T{T}}, b::T, k::Int) where
                        {T<:NumberNotSeriesN}
                    bb = k == 0 ? b : zero(b)
                    @inbounds for i in eachindex(v[k])
                        ($fc)(v[k], a[k], bb, i)
                    end
                    return nothing
                end

                function ($fc)(v::$T{$T{T}}, a::T, b::$T{$T{T}}, k::Int) where
                        {T<:NumberNotSeriesN}
                    aa = k == 0 ? a : zero(a)
                    @inbounds for i in eachindex(v[k])
                        ($fc)(v[k], aa, b[k], i)
                    end
                    return nothing
                end

            end
        else
            @eval begin
                function ($f)(a::$T{T}, b::S) where {T<:Number, S<:NumberNotSeries}
                    R = promote_type(T, S)
                    aa = convert(TaylorN{R}, a)
                    bb = convert(R, b)
                    c = TaylorN(aa.coeffs[:], order(aa))
                    c[0][1] = ($f)(aa[0][1], bb)
                    return c
                end

                function ($f)(b::S, a::$T{T}) where {T<:Number, S<:NumberNotSeries}
                    R = promote_type(T, S)
                    aa = convert(TaylorN{R}, a)
                    bb = convert(R, b)
                    c = TaylorN(($f)(aa.coeffs[:]), order(aa))
                    c[0][1] = ($f)(bb, aa[0][1])
                    return c
                end

                function ($fc)(v::$T{T}, a::$T{T}, k::Int) where {T<:Number}
                    @inbounds for l in eachindex(v[k])
                        v[k][l] = ($f)(a[k][l])
                    end
                    return nothing
                end

                function ($fc)(v::$T, a::$T, b::$T, k::Int)
                    _check_same_space(v, a, b)
                    @inbounds for i in eachindex(v[k])
                        v[k][i] = ($f)(a[k][i], b[k][i])
                    end
                    return nothing
                end

                function ($fc)(v::$T, a::$T, b::Number, k::Int)
                    bb = k == 0 ? b : zero(b)
                    for i in eachindex(v[k])
                        v[k][i] = ($f)(a[k][i], bb)
                    end
                    return nothing
                end

                function ($fc)(v::$T, a::Number, b::$T, k::Int)
                    aa = k == 0 ? a : zero(a)
                    for i in eachindex(v[k])
                        v[k][i] = ($f)(aa, b[k][i])
                    end
                    return nothing
                end

            end
        end
    end

    @eval ($f)(a::T, b::S) where {T<:Taylor1, S<:TaylorN} = ($f)(promote(a, b)...)
    @eval ($f)(a::T, b::S) where {T<:TaylorN, S<:Taylor1} = ($f)(promote(a, b)...)

    @eval begin
        function ($f)(a::HomogeneousPolynomial{T}, b::HomogeneousPolynomial{S}) where
                {T<:NumberNotSeriesN, S<:NumberNotSeriesN}
            _check_same_space(a, b)
            @assert order(a) == order(b)
            v = ($f)(a.coeffs, b.coeffs)
            return HomogeneousPolynomial(a.space, v, order(a))
        end

        function ($f)(a::HomogeneousPolynomial{T}, b::HomogeneousPolynomial{T}) where
                {T<:NumberNotSeriesN}
            _check_same_space(a, b)
            @assert order(a) == order(b)
            v = ($f)(a.coeffs, b.coeffs)
            return HomogeneousPolynomial(a.space, v, order(a))
        end

        # NOTE add! and subst! for HomogeneousPolynomial's act as += or -=
        function ($fc)(res::HomogeneousPolynomial{T}, a::HomogeneousPolynomial{T},
                b::HomogeneousPolynomial{T}, k::Int) where {T<:NumberNotSeriesN}
            _check_same_space(res, a, b)
            res[k] += ($f)(a[k], b[k])
            return nothing
        end

        ($f)(a::HomogeneousPolynomial) =
            HomogeneousPolynomial(a.space, ($f)(a.coeffs), order(a))

        function ($f)(a::TaylorN{Taylor1{T}}, b::S) where
                {T<:NumberNotSeries, S<:NumberNotSeries}
            @inbounds aux = $f(a[0][1], b)
            R = TS.numtype(aux)
            coeffs = FixedSizeVectorDefault{HomogeneousPolynomial{Taylor1{R}}}(
                    undef, order(a)+1)
            coeffs .= a.coeffs
            @inbounds coeffs[1] = aux
            return TaylorN(coeffs, order(a))
        end

        function ($f)(b::S, a::TaylorN{Taylor1{T}}) where
                {T<:NumberNotSeries, S<:NumberNotSeries}
            @inbounds aux = $f(b, a[0][1])
            R = TS.numtype(aux)
            coeffs = FixedSizeVectorDefault{HomogeneousPolynomial{Taylor1{R}}}(
                    undef, order(a)+1)
            coeffs .= $f.(a.coeffs)
            @inbounds coeffs[1] = aux
            return TaylorN(coeffs, order(a))
        end

        function ($f)(a::TaylorN{Taylor1{T}}, b::Taylor1{S}) where
                {T<:NumberNotSeries, S<:NumberNotSeries}
            @inbounds aux = $f(a[0][1], b)
            R = TS.numtype(aux)
            coeffs = FixedSizeVectorDefault{HomogeneousPolynomial{Taylor1{R}}}(undef, order(a)+1)
            coeffs .= a.coeffs
            @inbounds coeffs[1] = aux
            return TaylorN(coeffs, order(a))
        end

        function ($f)(b::Taylor1{S}, a::TaylorN{Taylor1{T}}) where
                {T<:NumberNotSeries, S<:NumberNotSeries}
            @inbounds aux = $f(b, a[0][1])
            R = TS.numtype(aux)
            coeffs = FixedSizeVectorDefault{HomogeneousPolynomial{Taylor1{R}}}(undef, order(a)+1)
            coeffs .= $f.(a.coeffs)
            @inbounds coeffs[1] = aux
            return TaylorN(coeffs, order(a))
        end

        function ($f)(a::Taylor1{TaylorN{T}}, b::S) where
                {T<:NumberNotSeries, S<:NumberNotSeries}
            @inbounds aux = ($f)(a[0][0][1], b)
            c = Taylor1( TaylorN(a[0].space, zero(aux), order(a[0])), order(a))
            for k in eachindex(a)
                ($fc)(c, a, b, k)
            end
            return c
        end

        function ($f)(b::S, a::Taylor1{TaylorN{T}}) where
                {T<:NumberNotSeries, S<:NumberNotSeries}
            @inbounds aux = ($f)(b, a[0][0][1])
            c = Taylor1( TaylorN(a[0].space, zero(aux), order(a[0])), order(a))
            for k in eachindex(a)
                ($fc)(c, b, a, k)
            end
            return c
        end

        function ($f)(a::Taylor1{TaylorN{T}}, b::TaylorN{S}) where
                {T<:NumberNotSeries, S<:NumberNotSeries}
            _check_same_space(a[0], b)
            @inbounds aux = $f(a[0], b)
            c = Taylor1( zero(aux), order(a))
            for k in eachindex(a)
                ($fc)(c, a, b, k)
            end
            return c
        end

        function ($f)(b::TaylorN{S}, a::Taylor1{TaylorN{T}}) where
                {T<:NumberNotSeries,S<:NumberNotSeries}
            _check_same_space(b, a[0])
            @inbounds aux = $f(b, a[0])
            c = Taylor1( zero(aux), order(a))
            for k in eachindex(a)
                ($fc)(c, a, b, k)
            end
            return c
        end

    end
end

function add!(v::Taylor1{T}, a::Taylor1{T}, b::Taylor1{T}) where
        {T<:NumberNotSeries}
    v_coeffs = v.coeffs
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    @inbounds for i in eachindex(v_coeffs)
        v_coeffs[i] = a_coeffs[i] + b_coeffs[i]
    end
    return nothing
end

function subst!(v::Taylor1{T}, a::Taylor1{T}, b::Taylor1{T}) where
        {T<:NumberNotSeries}
    v_coeffs = v.coeffs
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    @inbounds for i in eachindex(v_coeffs)
        v_coeffs[i] = a_coeffs[i] - b_coeffs[i]
    end
    return nothing
end

function add!(v::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}},
        b::Taylor1{Taylor1{T}}) where {T<:NumberNotSeries}
    v_coeffs = v.coeffs
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    @inbounds for i in eachindex(v_coeffs)
        add!(v_coeffs[i], a_coeffs[i], b_coeffs[i])
    end
    return nothing
end

function subst!(v::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}},
        b::Taylor1{Taylor1{T}}) where {T<:NumberNotSeries}
    v_coeffs = v.coeffs
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    @inbounds for i in eachindex(v_coeffs)
        subst!(v_coeffs[i], a_coeffs[i], b_coeffs[i])
    end
    return nothing
end

function add!(v::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}},
        b::Taylor1{Taylor1{T}}, k::Int) where {T<:NumberNotSeries}
    @inbounds add!(v.coeffs[k+1], a.coeffs[k+1], b.coeffs[k+1])
    return nothing
end

function subst!(v::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}},
        b::Taylor1{Taylor1{T}}, k::Int) where {T<:NumberNotSeries}
    @inbounds subst!(v.coeffs[k+1], a.coeffs[k+1], b.coeffs[k+1])
    return nothing
end

function +(a::Taylor1{T}, b::Taylor1{T}) where {T<:NumberNotSeries}
    if order(a) != order(b)
        a, b = fixorder(a, b)
    end
    c = zero(a)
    add!(c, a, b)
    return c
end

function -(a::Taylor1{T}, b::Taylor1{T}) where {T<:NumberNotSeries}
    if order(a) != order(b)
        a, b = fixorder(a, b)
    end
    c = zero(a)
    subst!(c, a, b)
    return c
end

function +(a::Taylor1{Taylor1{T}}, b::Taylor1{Taylor1{T}}) where
        {T<:NumberNotSeries}
    if order(a) != order(b)
        a, b = fixorder(a, b)
    end
    c = zero(a)
    add!(c, a, b)
    return c
end

function -(a::Taylor1{Taylor1{T}}, b::Taylor1{Taylor1{T}}) where
        {T<:NumberNotSeries}
    if order(a) != order(b)
        a, b = fixorder(a, b)
    end
    c = zero(a)
    subst!(c, a, b)
    return c
end

for (f, fc) in ((:+, :(add!)), (:-, :(subst!)))
    @eval begin
        function ($fc)(v::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
                b::Taylor1{TaylorN{T}}) where {T<:NumberNotSeries}
            v_coeffs = v.coeffs
            a_coeffs = a.coeffs
            b_coeffs = b.coeffs
            _check_same_space(v_coeffs[1], a_coeffs[1], b_coeffs[1])
            @inbounds for i in eachindex(v_coeffs)
                v_hps = v_coeffs[i].coeffs
                a_hps = a_coeffs[i].coeffs
                b_hps = b_coeffs[i].coeffs
                for j in eachindex(v_hps)
                    v_hp = v_hps[j].coeffs
                    a_hp = a_hps[j].coeffs
                    b_hp = b_hps[j].coeffs
                    for k in eachindex(v_hp)
                        v_hp[k] = ($f)(a_hp[k], b_hp[k])
                    end
                end
            end
            return nothing
        end
        function ($f)(a::Taylor1{TaylorN{T}}, b::Taylor1{TaylorN{T}}) where
                {T<:NumberNotSeries}
            _check_same_space(a[0], b[0])
            if order(a) != order(b) ||
                    any(order.(a.coeffs) .!= order.(b.coeffs))
                a, b = fixorder(a, b)
            end
            c = zero(a)
            ($fc)(c, a, b)
            return c
        end
        function ($fc)(v::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
                b::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
            v_k = v.coeffs[k+1]
            a_k = a.coeffs[k+1]
            b_k = b.coeffs[k+1]
            v_hps = v_k.coeffs
            a_hps = a_k.coeffs
            b_hps = b_k.coeffs
            @inbounds for i in eachindex(v_hps)
                v_hp = v_hps[i].coeffs
                a_hp = a_hps[i].coeffs
                b_hp = b_hps[i].coeffs
                for j in eachindex(v_hp)
                    v_hp[j] = ($f)(a_hp[j], b_hp[j])
                end
            end
            return nothing
        end
        function ($fc)(v::Taylor1{TaylorN{T}}, a::NumberNotSeries,
                b::Taylor1{TaylorN{T}}, k::Int) where {T<:NumberNotSeries}
            v_k = v.coeffs[k+1]
            b_k = b.coeffs[k+1]
            v_hps = v_k.coeffs
            b_hps = b_k.coeffs
            @inbounds for i in eachindex(v_hps)
                aaa = ifelse(k == 0 && i == 1, a, zero(a))
                v_hp = v_hps[i].coeffs
                b_hp = b_hps[i].coeffs
                for j in eachindex(v_hp)
                    v_hp[j] = ($f)(aaa, b_hp[j])
                end
            end
            return nothing
        end
        function ($fc)(v::Taylor1{TaylorN{T}}, b::Taylor1{TaylorN{T}},
                a::NumberNotSeries, k::Int) where {T<:NumberNotSeries}
            v_k = v.coeffs[k+1]
            b_k = b.coeffs[k+1]
            v_hps = v_k.coeffs
            b_hps = b_k.coeffs
            @inbounds for i in eachindex(v_hps)
                aaa = ifelse(k == 0 && i == 1, a, zero(a))
                v_hp = v_hps[i].coeffs
                b_hp = b_hps[i].coeffs
                for j in eachindex(v_hp)
                    v_hp[j] = ($f)(b_hp[j], aaa)
                end
            end
            return nothing
        end
        function ($fc)(v::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
                k::Int) where {T<:NumberNotSeries}
            v_k = v.coeffs[k+1]
            a_k = a.coeffs[k+1]
            v_hps = v_k.coeffs
            a_hps = a_k.coeffs
            @inbounds for l in eachindex(v_hps)
                v_hp = v_hps[l].coeffs
                a_hp = a_hps[l].coeffs
                for m in eachindex(v_hp)
                    v_hp[m] = ($f)(a_hp[m])
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



## Division ##
function /(a::Taylor1{Rational{T}}, b::S) where {T<:Integer, S<:NumberNotSeries}
    R = typeof( a[0] // b)
    v = FixedSizeVectorDefault{R}(undef, order(a)+1)
    v .= a.coeffs .// b
    return Taylor1(v, order(a))
end

function /(a::Taylor1{T}, b::S) where {T<:NumberNotSeries, S<:NumberNotSeries}
    @inbounds aux = a.coeffs[1] / b
    v = FixedSizeVectorDefault{typeof(aux)}(undef, length(a.coeffs))
    v .= a.coeffs ./ b
    return Taylor1(v, order(a))
end

for T in (:HomogeneousPolynomial, :TaylorN)
    @eval function /(a::$T{T}, b::S) where {T<:NumberNotSeries, S<:NumberNotSeries}
        @inbounds aux = a.coeffs[1] / b
        v = FixedSizeVectorDefault{typeof(aux)}(undef, length(a.coeffs))
        v .= a.coeffs ./ b
        return $T(a.space, v, order(a))
    end

    # @eval function /(a::$T{T}, b::T) where {T<:Number}
    #     @inbounds aux = a.coeffs[1] / b
    #     c = $T( zero(aux), order(a) )
    #     for ord in eachindex(c)
    #         div!(c, a, b, ord) # updates c[ord]
    #     end
    #     return c
    # end
end

for T in (:HomogeneousPolynomial, :TaylorN)
    @eval function /(b::$T{Taylor1{S}}, a::Taylor1{T}) where
            {T<:NumberNotSeries, S<:NumberNotSeries}
        @inbounds aux = b.coeffs[1] / a
        R = typeof(aux)
        coeffs = FixedSizeVectorDefault{R}(undef, length(b.coeffs))
        coeffs .= b.coeffs ./ a
        return $T(b.space, coeffs, order(b))
    end

    @eval function /(b::$T{Taylor1{T}}, a::S) where {T<:NumberNotSeries, S<:NumberNotSeries}
        @inbounds aux = b.coeffs[1] / a
        R = typeof(aux)
        coeffs = FixedSizeVectorDefault{R}(undef, length(b.coeffs))
        coeffs .= b.coeffs ./ a
        return $T(b.space, coeffs, order(b))
    end

    @eval function /(b::Taylor1{$T{S}}, a::$T{T}) where
            {T<:NumberNotSeries, S<:NumberNotSeries}
        @inbounds aux = b[0] / a
        v = Taylor1(zero(aux), order(b))
        @inbounds for k in eachindex(b)
            v[k] = b[k] / a
        end
        return v
    end
end

/(a::Taylor1{T}, b::Taylor1{S}) where {T<:Number, S<:Number} = /(promote(a,b)...)

function /(a::Taylor1{T}, b::Taylor1{T}) where {T<:Number}
    iszero(a) && !iszero(b) && return zero(a)
    if order(a) != order(b)
        a, b = fixorder(a, b)
    end

    # order and coefficient of first factorized term
    ordfact, cdivfact = divfactorization(a, b)
    R = typeof(cdivfact)
    if R == T
        aa = a
        bb = b
    else
        aa = convert(Taylor1{R}, a)
        bb = convert(Taylor1{R}, b)
    end
    c = Taylor1(cdivfact, order(a)-ordfact)
    for ord in eachindex(c)
        div!(c, aa, bb, ord) # updates c[ord]
    end

    return c
end

/(a::TaylorN{T}, b::TaylorN{S}) where
    {T<:NumberNotSeriesN, S<:NumberNotSeriesN} = /(promote(a,b)...)

function /(a::TaylorN{T}, b::TaylorN{T}) where {T<:NumberNotSeriesN}
    _check_same_space(a, b)
    @assert !iszero(constant_term(b))

    if order(a) != order(b)
        a, b = fixorder(a, b)
    end

    # first coefficient
    @inbounds cdivfact = a[0] / constant_term(b)
    c = TaylorN(cdivfact, order(a))
    for ord in eachindex(c)
        div!(c, a, b, ord) # updates c[ord]
    end

    return c
end

function /(a::S, b::TaylorN{T}) where {S<:NumberNotSeriesN, T<:NumberNotSeriesN}
    @assert !iszero(constant_term(b))
    R = typeof(a / constant_term(b))
    bb = convert(TaylorN{R}, b)
    res = TaylorN(b.space, zero(R), order(b))
    iszero(a) && !iszero(b) && return res
    aa = convert(R, a)
    for ord in eachindex(res)
        div!(res, aa, bb, ord)
    end
    return res
end

function /(a::Taylor1{TaylorN{T}}, b::Taylor1{TaylorN{T}}) where {T<:NumberNotSeries}
    _check_same_space(a[0], b[0])
    iszero(a) && !iszero(b) && return zero(a)
    if (order(a) != order(b)) || any(order.(a.coeffs) .!= order.(b.coeffs))
        a, b = fixorder(a, b)
    end
    # order and coefficient of first factorized term
    ordfact, cdivfact = divfactorization(a, b)
    R = numtype(cdivfact)
    if R == T
        aa = a
        bb = b
    else
        aa = convert(Taylor1{TaylorN{R}}, a)
        bb = convert(Taylor1{TaylorN{R}}, b)
    end
    res = Taylor1(cdivfact, order(a)-ordfact)
    for ordT in eachindex(res)
        div!(res, aa, bb, ordT)
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

    aa = Taylor1(a, order(b))
    for ordT in eachindex(res)
        div!(res, aa, b, ordT)
    end
    return res
end

# @inline
function divfactorization(a1::Taylor1, b1::Taylor1)
    # order of first factorized term; a1 and b1 assumed to be of the same order
    a1nz = findfirst(a1)
    b1nz = findfirst(b1)
    a1nz = a1nz ≥ 0 ? a1nz : order(a1)
    b1nz = b1nz ≥ 0 ? b1nz : order(a1)
    ordfact = min(a1nz, b1nz)
    cdivfact = a1[ordfact] / b1[ordfact]

    # Is the polynomial factorizable?
    TS._isthinzero(b1[ordfact]) && throw( ArgumentError(
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

# @inline
function div!(c::Taylor1{T}, a::Taylor1{T}, b::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    kk = k+1
    @inbounds c_coeffs[kk] = zero(c_coeffs[kk])
    iszero(a) && !iszero(b) && return nothing
    # order and coefficient of first factorized term
    ordfact, cdivfact = divfactorization(a, b)
    if k == 0
        @inbounds c_coeffs[1] = cdivfact
        return nothing
    end
    b_order = order(b)
    imin = max(0, k+ordfact-b_order)
    @inbounds acc = c_coeffs[imin+1] * b_coeffs[k+ordfact-imin+1]
    @inbounds for i = imin+1:k-1
        acc += c_coeffs[i+1] * b_coeffs[k+ordfact-i+1]
    end
    if k+ordfact ≤ b_order
        @inbounds acc = a_coeffs[k+ordfact+1] - acc
    else
        acc = -acc
    end
    @inbounds c_coeffs[kk] = acc / b_coeffs[ordfact+1]
    return nothing
end

# @inline
function div!(v::Taylor1{T}, a::Taylor1{T}, b::NumberNotSeries,
        k::Int) where {T<:Number}
    @inbounds v.coeffs[k+1] = a.coeffs[k+1] / b
    return nothing
end

function div!(v::Taylor1{T}, a::Taylor1{S}, b::NumberNotSeries) where
        {T<:NumberNotSeries, S<:NumberNotSeries}
    v_coeffs = v.coeffs
    a_coeffs = a.coeffs
    @inbounds for i in eachindex(v_coeffs)
        v_coeffs[i] = a_coeffs[i] / b
    end
    return nothing
end

function div!(v::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}},
        b::NumberNotSeries) where {T<:NumberNotSeries}
    v_coeffs = v.coeffs
    a_coeffs = a.coeffs
    @inbounds for i in eachindex(v_coeffs)
        div!(v_coeffs[i], a_coeffs[i], b)
    end
    return nothing
end

# @inline
function div!(c::Taylor1{T}, a::NumberNotSeries, b::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    c_coeffs = c.coeffs
    b_coeffs = b.coeffs
    kk = k+1
    @inbounds c_coeffs[kk] = zero(c_coeffs[kk])
    iszero(a) && !iszero(b) && return nothing
    # order and coefficient of first factorized term
    # In this case, since a[k]=0 for k>0, we can simplify to:
    # ordfact, cdivfact = 0, a/b[0]
    if k == 0
        @inbounds c_coeffs[1] = a / b_coeffs[1]
        return nothing
    end
    @inbounds acc = c_coeffs[1] * b_coeffs[kk]
    @inbounds for i = 1:k-1
        acc += c_coeffs[i+1] * b_coeffs[k-i+1]
    end
    @inbounds c_coeffs[kk] = -acc / b_coeffs[1]
    return nothing
end

function div!(c::Taylor1, a::NumberNotSeries, b::Taylor1)
    @inbounds for k in eachindex(c)
        div!(c, a, b, k)
    end
    return nothing
end

@inline function div!(c::Taylor1{Taylor1{T}}, a::NumberNotSeries,
        b::Taylor1{Taylor1{T}}, k::Int) where {T<:NumberNotSeriesN}
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
        muladd!(c[k], c[i], b[k-i])
    end
    # @inbounds c[k] = -c[k]/b[0]
    @inbounds div_scalar!(c[k], -1, b[0])
    return nothing
end

#
# @inline
function div!(c::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}},
        b::Taylor1{Taylor1{T}}, k::Int) where {T<:NumberNotSeriesN}
    zero!(c, k)
    iszero(a) && !iszero(b) && return nothing
    # order and coefficient of first factorized term
    ordfact, aux = divfactorization(a, b)
    if k == 0
        identity!(c, aux, k)
        return nothing
    end
    b_order = order(b)
    imin = max(0, k+ordfact-b_order)
    # c[k] = c[imin] * b[k+ordfact-imin]
    mul!(c[k], c[imin], b[k+ordfact-imin])
    for i = imin+1:k-1
        # c[k] += c[i] * b[k+ordfact-i]
        for ord in eachindex(minlength(c[k], c[i], b[k+ordfact-i]))
            muladd!(c[k], c[i], b[k+ordfact-i], ord)
        end
    end
    zero!(aux)
    if k+ordfact ≤ b_order
        # @inbounds aux <- a[k+ordfact] - c[k]
        for ord in eachindex(minlength(aux, a[k+ordfact]))
            subst!(aux, a[k+ordfact], c[k], ord)
        end
    else
        # @inbounds aux <- - c[k]
        for ord in eachindex(minlength(aux, a[k+ordfact]))
            subst!(aux, c[k], ord)
        end
    end
    # c[k] <- aux / b[ordfact]
    for ord in eachindex(c[k])
        div!(c[k], aux, b[ordfact], ord)
    end
    return nothing
end

# @inline function div!(v::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}},
#         b::NumberNotSeries, k::Int) where {T<:NumberNotSeriesN}
#     # @inbounds v[k] = a[k] / b
#     for ord in eachindex(v)
#         div!(v, a, b, ord)
#     end
#     return nothing
# end

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
    _check_same_space(v, a)
    @inbounds for k in eachindex(v)
        v[k] = a[k] / b
    end
    return nothing
end

# NOTE: Due to the use of `zero!`, this `div!` method does *not* accumulate the result of a / b in c[k] (k > 0)
@inline function div!(c::TaylorN, a::TaylorN, b::TaylorN, k::Int)
    _check_same_space(c, a, b)
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
    _check_same_space(c, a)
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
    _check_same_space(c, a)
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

@inline function div_scalar!(c::Taylor1{T}, scalar::NumberNotSeries,
        a::Taylor1{T}, k::Int) where {T <: NumberNotSeries}
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    if k==0
        @inbounds c_coeffs[1] = scalar*c_coeffs[1] / a_coeffs[1]
        return nothing
    end

    kk = k+1
    @inbounds aux = scalar * c_coeffs[kk]
    @inbounds acc = zero(c_coeffs[kk])
    @inbounds for i = 0:k-1
        acc -= c_coeffs[i+1] * a_coeffs[k-i+1]
    end
    @inbounds c_coeffs[kk] = (acc + aux) / a_coeffs[1]
    return nothing
end

# NOTE: Here `div!` *accumulates* the result of a[k] / b[k] in c[k] (k > 0)
@inline function div!(c::TaylorN, a::NumberNotSeries, b::TaylorN, k::Int)
    _check_same_space(c, b)
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

# c[k] <- a[k]/b, where b is a scalar
@inline function div!(c::TaylorN{T}, a::TaylorN{T}, b::NumberNotSeries,
        k::Int) where {T<:Number}
    _check_same_space(c, a)
    @inbounds for i in eachindex(c[k])
        c[k][i] = a[k][i] / b
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

# in-place division c <- scalar*c/a (assumes equal order among TaylorNs)
function div_scalar!(c::Taylor1, scalar::NumberNotSeries, a::Taylor1)
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
function div!(c::TaylorN{T}, a::TaylorN{T}, b::NumberNotSeries) where {T<:Number}
    @inbounds for k in eachindex(c)
        div!(c, a, b, k)
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
    anz = anz ≥ 0 ? anz : order(a)
    bnz = bnz ≥ 0 ? bnz : order(a)
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
    b_order = order(b)
    imin = max(0, k+ordfact-b_order)
    @inbounds mul!(c[k], c[imin], b[k+ordfact-imin])
    @inbounds for i = imin+1:k-1
        mul!(c[k], c[i], b[k+ordfact-i])
    end
        if k+ordfact ≤ b_order
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
    res_k = res.coeffs[k+1]
    a_k = a.coeffs[k+1]
    res_hps = res_k.coeffs
    a_hps = a_k.coeffs
    @inbounds for l in eachindex(res_hps)
        res_hp = res_hps[l].coeffs
        a_hp = a_hps[l].coeffs
        for m in eachindex(res_hp)
            res_hp[m] = a_hp[m]/b
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
    order = maximum(TS.order.(b))

    # Use matrices of coefficients (of proper size) and mul!
    # B = zeros(T, k, order+1)
    B = Array{T}(undef, k, order+1)
    B = zero.(B)
    for i = 1:k
        @inbounds ord = TS.order(b[i])
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
const LU_RowMaximum = RowMaximum()
const LU_NoPivot = NoPivot()

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


# Fast allocation-free matrix multiplication
for T in (:Taylor1, :TaylorN)
    @eval function matmul!(C::Matrix{$T{T}},
                           A::Matrix{$T{T}}, B::Matrix{$T{T}}) where {T}
        mc, nc = size(C)
        ma, na = size(A)
        mb, nb = size(B)
        @assert (na == mb && mc == ma && nc == nb)
        for j in axes(C,2)
            for i in axes(C,1)
                TS.zero!(C[i,j])
                for k in 1:na
                    TS.muladd!(C[i,j], A[i,k], B[k,j])
                end
            end
        end
        return nothing
    end
end
