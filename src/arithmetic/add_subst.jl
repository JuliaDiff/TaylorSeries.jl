# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

## Addition and subtraction ##
for (f, fc) in ((:+, :(add!)), (:-, :(subst!)))

    for T in (:Taylor1, :TaylorN)
        @eval begin
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

            ## add! and subst! ##
            function ($fc)(v::$T{T}, a::T, k::Int) where {T<:Number}
                @inbounds v[k] = k==0 ? ($f)(a) : zero(a)
                return nothing
            end
        end

        if T == :Taylor1
            @eval begin
                function ($f)(a::$T{T}, b::S) where {T<:Number, S<:NumberNotSeries}
                    c = $T(($f)(a.coeffs[1], b), order(a))
                    for k in eachindex(a)
                        ($fc)(c, a, b, k)
                    end
                    return c
                end

                function ($f)(b::S, a::$T{T}) where {T<:Number, S<:NumberNotSeries}
                    c = $T(($f)(b, a.coeffs[1]), order(a))
                    for k in eachindex(a)
                        ($fc)(c, b, a, k)
                    end
                    return c
                end

                function ($f)(a::$T{T}, b::$T{S}) where {T<:NumberNotSeries, S<:NumberNotSeries}
                    if order(a) != order(b)
                        a, b = fixorder(a, b)
                    end
                    z = zero(a.coeffs[1] + b.coeffs[1])
                    c = $T(z, order(a))
                    for k in eachindex(a)
                        ($fc)(c, a, b, k)
                    end
                    return c
                end

                function ($fc)(v::$T{T}, a::$T{T}, k::Int) where {T<:Number}
                    @inbounds v.coeffs[k+1] = ($f)(a.coeffs[k+1])
                    return nothing
                end

                @inline function ($fc)(v::$T, a::$T, b::$T, k::Int)
                    v_coeffs = v.coeffs
                    a_coeffs = a.coeffs
                    b_coeffs = b.coeffs
                    kk = k + 1
                    @inbounds v_coeffs[kk] = ($f)(a_coeffs[kk], b_coeffs[kk])
                    return nothing
                end

                function ($fc)(v::$T, a::$T, b::Number, k::Int)
                    bb = k==0 ? b : zero(b)
                    v.coeffs[k+1] = ($f)(a.coeffs[k+1], bb)
                    return nothing
                end

                function ($fc)(v::$T, a::Number, b::$T, k::Int)
                    aa = k==0 ? a : zero(a)
                    v.coeffs[k+1] = ($f)(aa, b.coeffs[k+1])
                    return nothing
                end

                # Nested Taylor1s
                function ($fc)(v::$T{$T{T}}, a::$T{$T{T}}, k::Int) where
                        {T<:NumberNotSeriesN}
                    v_coeffs = v.coeffs
                    a_coeffs = a.coeffs
                    kk = k + 1
                    @inbounds for i in eachindex(v_coeffs[kk])
                        ($fc)(v_coeffs[kk], a_coeffs[kk], i)
                    end
                    return nothing
                end

                function ($fc)(v::$T{$T{T}}, a::$T{$T{T}}, b::$T{$T{T}}, k::Int) where
                        {T<:NumberNotSeriesN}
                    v_coeffs = v.coeffs
                    a_coeffs = a.coeffs
                    b_coeffs = b.coeffs
                    kk = k + 1
                    @inbounds for i in eachindex(v_coeffs[kk])
                        ($fc)(v_coeffs[kk], a_coeffs[kk], b_coeffs[kk], i)
                    end
                    return nothing
                end

                function ($fc)(v::$T{$T{T}}, a::$T{$T{T}}, b::$T{T}, k::Int) where
                        {T<:NumberNotSeriesN}
                    v_coeffs = v.coeffs
                    a_coeffs = a.coeffs
                    kk = k + 1
                    @inbounds for i in eachindex(v_coeffs[kk])
                        ($fc)(v_coeffs[kk], a_coeffs[kk], b, i)
                    end
                    return nothing
                end

                function ($fc)(v::$T{$T{T}}, a::$T{T}, b::$T{$T{T}}, k::Int) where
                        {T<:NumberNotSeriesN}
                    v_coeffs = v.coeffs
                    b_coeffs = b.coeffs
                    kk = k + 1
                    @inbounds for i in eachindex(v_coeffs[kk])
                        ($fc)(v_coeffs[kk], a, b_coeffs[kk], i)
                    end
                    return nothing
                end

                function ($fc)(v::$T{$T{T}}, a::$T{$T{T}}, b::T, k::Int) where
                        {T<:NumberNotSeriesN}
                    bb = k == 0 ? b : zero(b)
                    v_coeffs = v.coeffs
                    a_coeffs = a.coeffs
                    kk = k + 1
                    @inbounds for i in eachindex(v_coeffs[kk])
                        ($fc)(v_coeffs[kk], a_coeffs[kk], bb, i)
                    end
                    return nothing
                end

                function ($fc)(v::$T{$T{T}}, a::T, b::$T{$T{T}}, k::Int) where
                        {T<:NumberNotSeriesN}
                    aa = k == 0 ? a : zero(a)
                    v_coeffs = v.coeffs
                    b_coeffs = b.coeffs
                    kk = k + 1
                    @inbounds for i in eachindex(v_coeffs[kk])
                        ($fc)(v_coeffs[kk], aa, b_coeffs[kk], i)
                    end
                    return nothing
                end

                function ($fc)(v::$T{T}, a::$T{T}, b::$T{T}) where
                        {T<:NumberNotSeries}
                    for k in eachindex(v)
                        ($fc)(v, a, b, k)
                    end
                    return nothing
                end

                function ($fc)(v::$T{$T{T}}, a::$T{$T{T}},
                        b::$T{$T{T}}) where {T<:NumberNotSeriesN}
                    v_coeffs = v.coeffs
                    a_coeffs = a.coeffs
                    b_coeffs = b.coeffs
                    @inbounds for i in eachindex(v_coeffs)
                        ($fc)(v_coeffs[i], a_coeffs[i], b_coeffs[i])
                    end
                    return nothing
                end

                function ($fc)(v::$T{$T{T}}, a::$T{$T{T}},
                        b::$T{$T{T}}, k::Int) where {T<:NumberNotSeries}
                    v_coeffs = v.coeffs
                    a_coeffs = a.coeffs
                    b_coeffs = b.coeffs
                    kk = k + 1
                    @inbounds ($fc)(v_coeffs[kk], a_coeffs[kk], b_coeffs[kk])
                    return nothing
                end

                function ($f)(a::$T{T}, b::$T{T}) where {T<:NumberNotSeries}
                    if order(a) != order(b)
                        a, b = fixorder(a, b)
                    end
                    c = zero(a)
                    ($fc)(c, a, b)
                    return c
                end

                function ($f)(a::$T{$T{T}}, b::$T{$T{T}}) where
                        {T<:NumberNotSeries}
                    if order(a) != order(b)
                        a, b = fixorder(a, b)
                    end
                    c = zero(a)
                    ($fc)(c, a, b)
                    return c
                end

            end
        else # TaylorN
            @eval begin
                function ($f)(a::$T{T}, b::S) where {T<:Number, S<:NumberNotSeries}
                    c = $T(space(a), ($f)(constant_term(a), b), order(a))
                    for k in eachindex(c)
                        ($fc)(c, a, b, k)
                    end
                    return c
                end

                function ($f)(b::S, a::$T{T}) where {T<:Number, S<:NumberNotSeries}
                    c = $T(space(a), ($f)(b, constant_term(a)), order(a))
                    for k in eachindex(c)
                        ($fc)(c, b, a, k)
                    end
                    return c
                end

                function ($f)(a::$T{T}, b::$T{S}) where {T<:Number, S<:Number}
                    _check_same_space(a, b)
                    z = zero(a.coeffs[1] + b.coeffs[1])
                    c = $T(z, order(a))
                    for k in eachindex(a)
                        ($fc)(c, a, b, k)
                    end
                    return c
                end

                function ($fc)(v::$T{T}, a::$T{T}, k::Int) where {T<:Number}
                    v_coeffs = v.coeffs[k+1].coeffs
                    a_coeffs = a.coeffs[k+1].coeffs
                    @inbounds for l in eachindex(v_coeffs)
                        v_coeffs[l] = ($f)(a_coeffs[l])
                    end
                    return nothing
                end

                function ($fc)(v::$T, a::$T, b::$T, k::Int)
                    kk = k + 1
                    v_coeffs = v.coeffs[kk].coeffs
                    a_coeffs = a.coeffs[kk].coeffs
                    b_coeffs = b.coeffs[kk].coeffs
                    @inbounds for i in eachindex(v_coeffs)
                        v_coeffs[i] = ($f)(a_coeffs[i], b_coeffs[i])
                    end
                    return nothing
                end

                function ($fc)(v::$T, a::$T, b::Number, k::Int)
                    v_coeffs = v.coeffs[k+1].coeffs
                    a_coeffs = a.coeffs[k+1].coeffs
                    copyto!(v_coeffs, a_coeffs)
                    constant_term!(v, ($f)(constant_term(a), b))
                    return nothing
                end

                function ($fc)(v::$T, a::Number, b::$T, k::Int)
                    v_coeffs = v.coeffs[k+1].coeffs
                    b_coeffs = b.coeffs[k+1].coeffs
                    copyto!(v_coeffs, ($f)(b_coeffs))
                    constant_term!(v, ($f)(a, constant_term(b)))
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
            res.coeffs[k] += ($f)(a.coeffs[k], b.coeffs[k])
            return nothing
        end

        ($f)(a::HomogeneousPolynomial) =
            HomogeneousPolynomial(a.space, ($f).(a.coeffs), order(a))

        function ($f)(a::TaylorN{Taylor1{T}}, b::S) where
                {T<:NumberNotSeries, S<:NumberNotSeries}
            @inbounds aux = $f(a.coeffs[1].coeffs[1], b)
            R = TS.numtype(aux)
            coeffs = FixedSizeVectorDefault{HomogeneousPolynomial{Taylor1{R}}}(
                    undef, order(a)+1)
            coeffs .= a.coeffs
            c = TaylorN(space(a), coeffs, order(a))
            constant_term!(c, aux)
            return c
        end

        function ($f)(b::S, a::TaylorN{Taylor1{T}}) where
                {T<:NumberNotSeries, S<:NumberNotSeries}
            @inbounds aux = $f(b, a.coeffs[1].coeffs[1])
            R = TS.numtype(aux)
            coeffs = FixedSizeVectorDefault{HomogeneousPolynomial{Taylor1{R}}}(
                    undef, order(a)+1)
            coeffs .= ($f)(a.coeffs)
            c = TaylorN(space(a), coeffs, order(a))
            constant_term!(c, aux)
            return c
        end

        function ($f)(a::TaylorN{Taylor1{T}}, b::Taylor1{S}) where
                {T<:NumberNotSeries, S<:NumberNotSeries}
            @inbounds aux = $f(a.coeffs[1].coeffs[1], b)
            R = TS.numtype(aux)
            coeffs = FixedSizeVectorDefault{HomogeneousPolynomial{Taylor1{R}}}(
                    undef, order(a)+1)
            coeffs .= a.coeffs
            c = TaylorN(space(a), coeffs, order(a))
            constant_term!(c, aux)
            return c
        end

        function ($f)(b::Taylor1{S}, a::TaylorN{Taylor1{T}}) where
                {T<:NumberNotSeries, S<:NumberNotSeries}
            @inbounds aux = $f(b, a.coeffs[1].coeffs[1])
            R = TS.numtype(aux)
            coeffs = FixedSizeVectorDefault{HomogeneousPolynomial{Taylor1{R}}}(
                    undef, order(a)+1)
            coeffs .= ($f)(a.coeffs)
            c = TaylorN(space(a), coeffs, order(a))
            constant_term!(c, aux)
            return c
        end

        function ($f)(a::Taylor1{TaylorN{T}}, b::S) where
                {T<:NumberNotSeries, S<:NumberNotSeries}
            @inbounds aux = ($f)(a.coeffs[1].coeffs[1].coeffs[1], b)
            c = Taylor1(
                TaylorN(space(a.coeffs[1]), zero(aux), order(a.coeffs[1])), order(a))
            for k in eachindex(a)
                ($fc)(c, a, b, k)
            end
            return c
        end

        function ($f)(b::S, a::Taylor1{TaylorN{T}}) where
                {T<:NumberNotSeries, S<:NumberNotSeries}
            @inbounds aux = ($f)(b, a.coeffs[1].coeffs[1].coeffs[1])
            c = Taylor1(
                TaylorN(space(a.coeffs[1]), zero(aux), order(a.coeffs[1])), order(a))
            for k in eachindex(a)
                ($fc)(c, b, a, k)
            end
            return c
        end

        function ($f)(a::Taylor1{TaylorN{T}}, b::TaylorN{S}) where
                {T<:NumberNotSeries, S<:NumberNotSeries}
            _check_same_space(a.coeffs[1], b)
            @inbounds aux = $f(a.coeffs[1], b)
            c = Taylor1( zero(aux), order(a))
            for k in eachindex(a)
                ($fc)(c, a, b, k)
            end
            return c
        end

        function ($f)(b::TaylorN{S}, a::Taylor1{TaylorN{T}}) where
                {T<:NumberNotSeries,S<:NumberNotSeries}
            _check_same_space(b, a.coeffs[1])
            @inbounds aux = $f(b, a.coeffs[1])
            c = Taylor1( zero(aux), order(a))
            for k in eachindex(a)
                ($fc)(c, a, b, k)
            end
            return c
        end

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
            kk = k+1
            v_hps = v.coeffs[kk].coeffs
            a_hps = a.coeffs[kk].coeffs
            b_hps = b.coeffs[kk].coeffs
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
            v_hps = v.coeffs[k+1].coeffs
            b_hps = b.coeffs[k+1].coeffs
            za = zero(a)
            @inbounds for i in eachindex(v_hps)
                aaa = ifelse(k == 0 && i == 1, a, za)
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
            v_hps = v.coeffs[k+1].coeffs
            b_hps = b.coeffs[k+1].coeffs
            za = zero(a)
            @inbounds for i in eachindex(v_hps)
                aaa = ifelse(k == 0 && i == 1, a, za)
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
            v_hps = v.coeffs[k+1].coeffs
            a_hps = a.coeffs[k+1].coeffs
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
