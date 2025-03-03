module TaylorSeriesIAExt

using TaylorSeries

import Base: ^, log, asin, acos, acosh, atanh, power_by_squaring

import TaylorSeries: evaluate, _evaluate, normalize_taylor, square

isdefined(Base, :get_extension) ? (using IntervalArithmetic) : (using ..IntervalArithmetic)

for T in (:Taylor1, :TaylorN)
    @eval begin
        @eval ^(a::$T{Interval{T}}, n::Integer) where {T<:Real} = TS.power_by_squaring(a, n)

        @eval ^(a::$T{Interval{T}}, r::Rational) where {T<:Real} = a^float(r)

        @eval function ^(a::$T{Interval{T}}, r::S) where {T<:Real, S<:Real}
            isinteger(r) && r ≥ 0 && return TS.power_by_squaring(a, Integer(r))
            a0 = constant_term(a) ∩ Interval(zero(T), T(Inf))
            @assert !isempty(a0)
            aux = one(a0^r)
            a[0] = aux * a0
            r == 0.5 && return sqrt(a)
            if $T == TaylorN
                if iszero(a0)
                    throw(DomainError(a,
                    """The 0-th order TaylorN coefficient must be non-zero
                    in order to expand `^` around 0."""))
                end
            end
            return TS._pow(a, r)
        end

        # _pow
        function TS._pow(a::$T{Interval{S}}, n::Integer) where {S<:Real}
            n < 0 && return TS._pow(a, float(n))
            return TS.power_by_squaring(a, n)
        end

        function TS._pow(a::$T{Interval{T}}, r::S) where {T<:Real, S<:Real}
            isinteger(r) && r ≥ 0 && return TS.power_by_squaring(a, Integer(r))
            a0 = constant_term(a) ∩ Interval(zero(T), T(Inf))
            @assert !isempty(a0)
            aux = one(a0^r)
            a[0] = aux * a0
            r == 0.5 && return sqrt(a)
            if $T == Taylor1
                l0 = findfirst(a)
                # Index of first non-zero coefficient of the result; must be integer
                !isinteger(r*l0) && throw(DomainError(a,
                    """The 0-th order Taylor1 coefficient must be non-zero
                    to raise the Taylor1 polynomial to a non-integer exponent."""))
                lnull = trunc(Int, r*l0 )
                (lnull > a.order) && return $T( zero(aux), a.order)
                c_order = l0 == 0 ? a.order : min(a.order, trunc(Int, r*a.order))
            else
                if iszero(a0)
                    throw(DomainError(a,
                    """The 0-th order TaylorN coefficient must be non-zero
                    in order to expand `^` around 0."""))
                end
                c_order = a.order
            end
            #
            c = $T(zero(aux), c_order)
            aux0 = zero(c)
            for k in eachindex(c)
                TS.pow!(c, a, aux0, r, k)
            end
            return c
        end

        function TS.pow!(res::$T{Interval{T}}, a::$T{Interval{T}}, aux::$T{Interval{T}},
                r::S, k::Int) where {T<:Real, S<:Integer}
            (r == 0) && return TS.one!(res, a, k)
            (r == 1) && return TS.identity!(res, a, k)
            (r == 2) && return TS.sqr!(res, a, k)
            TS.power_by_squaring!(res, a, aux, r)
            return nothing
        end
    end
end


for T in (:Taylor1, :TaylorN)
    @eval function log(a::$T{Interval{S}}) where {S<:Real}
        iszero(constant_term(a)) && throw(DomainError(a,
            """The 0-th order coefficient must be non-zero in order to expand `log` around 0."""))
        a0 = constant_term(a) ∩ Interval(S(0.0), S(Inf))
        order = a.order
        aux = log(a0)
        aa = one(aux) * a
        aa[0] = one(aux) * a0
        c = $T( aux, order )
        for k in eachindex(a)
            TS.log!(c, aa, k)
        end
        return c
    end

    @eval function asin(a::$T{Interval{S}}) where {S<:Real}
        a0 = constant_term(a) ∩ Interval(-one(S), one(S))
        a0^2 == one(a0) && throw(DomainError(a,
            """Series expansion of asin(x) diverges at x = ±1."""))

        order = a.order
        aux = asin(a0)
        aa = one(aux) * a
        aa[0] = one(aux) * a0
        c = $T( aux, order )
        r = $T( sqrt(1 - a0^2), order )
        for k in eachindex(a)
            TS.asin!(c, aa, r, k)
        end
        return c
    end

    @eval function acos(a::$T{Interval{S}}) where {S<:Real}
        a0 = constant_term(a) ∩ Interval(-one(S), one(S))
        a0^2 == one(a0) && throw(DomainError(a,
        """Series expansion of asin(x) diverges at x = ±1."""))

        order = a.order
        aux = acos(a0)
        aa = one(aux) * a
        aa[0] = one(aux) * a0
        c = $T( aux, order )
        r = $T( sqrt(1 - a0^2), order )
        for k in eachindex(a)
            TS.acos!(c, aa, r, k)
        end
        return c
    end

    @eval function acosh(a::$T{Interval{S}}) where {S<:Real}
        a0 = constant_term(a) ∩ Interval(one(S), S(Inf))
        a0^2 == one(a0) && throw(DomainError(a,
            """Series expansion of acosh(x) diverges at x = ±1."""))

        order = a.order
        aux = acosh(a0)
        aa = one(aux) * a
        aa[0] = one(aux) * a0
        c = $T( aux, order )
        r = $T( sqrt(a0^2 - 1), order )
        for k in eachindex(a)
            TS.acosh!(c, aa, r, k)
        end
        return c
    end

    @eval function atanh(a::$T{Interval{S}}) where {S<:Real}
        order = a.order
        a0 = constant_term(a) ∩ Interval(-one(S), one(S))
        aux = atanh(a0)
        aa = one(aux) * a
        aa[0] = one(aux) * a0
        c = $T( aux, order)
        r = $T(one(aux) - a0^2, order)
        iszero(constant_term(r)) && throw(DomainError(a,
            """Series expansion of atanh(x) diverges at x = ±1."""))

        for k in eachindex(a)
            TS.atanh!(c, aa, r, k)
        end
        return c
    end
end


function evaluate(a::Taylor1, dx::Interval{S}) where {S<:Real}
    order = a.order
    uno = one(dx)
    dx2 = dx^2
    if iseven(order)
        kend = order-2
        @inbounds sum_even = a[end]*uno
        @inbounds sum_odd = a[end-1]*zero(dx)
    else
        kend = order-3
        @inbounds sum_odd = a[end]*uno
        @inbounds sum_even = a[end-1]*uno
    end
    @inbounds for k in kend:-2:0
        sum_odd = sum_odd*dx2 + a[k+1]
        sum_even = sum_even*dx2 + a[k]
    end
    return sum_even + sum_odd*dx
end

function evaluate(a::TaylorN, dx::IntervalBox{N,T}) where {T<:Real,N}
    @assert N == get_numvars()
    a_length = length(a)
    suma = zero(constant_term(a)) + Interval{T}(0, 0)
    @inbounds for homPol in reverse(eachindex(a))
        suma += evaluate(a[homPol], dx)
    end

    return suma
end

function evaluate(a::HomogeneousPolynomial, dx::IntervalBox{N,T}) where {T<:Real,N}
    @assert N == get_numvars()
    dx == IntervalBox(-1..1, Val(N)) && return _evaluate(a, dx, Val(true))
    dx == IntervalBox( 0..1, Val(N)) && return _evaluate(a, dx, Val(false))

    return evaluate(a, dx...)
end

function _evaluate(a::HomogeneousPolynomial, dx::IntervalBox{N,T}, ::Val{true} ) where {T<:Real,N}
    a.order == 0 && return a[1] + Interval{T}(0, 0)

    ct = TS.coeff_table[a.order+1]
    @inbounds suma = a[1]*Interval{T}(0,0)

    Ieven = Interval{T}(0,1)
    for (i, a_coeff) in enumerate(a.coeffs)
        iszero(a_coeff) && continue
        if isodd(sum(ct[i]))
            tmp = dx[1]
        else
            tmp = Ieven
            for n in eachindex(ct[i])
                iseven(ct[i][n]) && continue
                tmp *= dx[1]
            end
        end
        suma += a_coeff * tmp
    end
    return suma
end

function _evaluate(a::HomogeneousPolynomial, dx::IntervalBox{N,T}, ::Val{false} ) where {T<:Real,N}
    a.order == 0 && return a[1] + Interval{T}(0, 0)

    @inbounds suma = zero(a[1])*dx[1]
    @inbounds for homPol in a.coeffs
        suma += homPol*dx[1]
    end
    return suma
end


"""
    normalize_taylor(a::Taylor1, I::Interval, symI::Bool=true)

Normalizes `a::Taylor1` such that the interval `I` is mapped
by an affine transformation to the interval `-1..1` (`symI=true`)
or to `0..1` (`symI=false`).
"""
normalize_taylor(a::Taylor1, I::Interval{T}, symI::Bool=true) where {T<:Real} =
    _normalize(a, I, Val(symI))

"""
    normalize_taylor(a::TaylorN, I::IntervalBox, symI::Bool=true)

Normalize `a::TaylorN` such that the intervals in `I::IntervalBox`
are mapped by an affine transformation to the intervals `-1..1`
(`symI=true`) or to `0..1` (`symI=false`).
"""
normalize_taylor(a::TaylorN, I::IntervalBox{N,T}, symI::Bool=true) where {T<:Real,N} =
    _normalize(a, I, Val(symI))

#  I -> -1..1
function _normalize(a::Taylor1, I::Interval{T}, ::Val{true}) where {T<:Real}
    order = get_order(a)
    t = Taylor1(T, order)
    tnew = mid(I) + t*radius(I)
    return a(tnew)
end

#  I -> 0..1
function _normalize(a::Taylor1, I::Interval{T}, ::Val{false}) where {T<:Real}
    order = get_order(a)
    t = Taylor1(T, order)
    tnew = inf(I) + t*diam(I)
    return a(tnew)
end

#  I -> IntervalBox(-1..1, Val(N))
function _normalize(a::TaylorN, I::IntervalBox{N,T}, ::Val{true}) where {T<:Real,N}
    order = get_order(a)
    x = Vector{typeof(a)}(undef, N)
    for ind in eachindex(x)
        x[ind] = mid(I[ind]) + TaylorN(ind, order=order)*radius(I[ind])
    end
    return a(x)
end

#  I -> IntervalBox(0..1, Val(N))
function _normalize(a::TaylorN, I::IntervalBox{N,T}, ::Val{false}) where {T<:Real,N}
    order = get_order(a)
    x = Vector{typeof(a)}(undef, N)
    for ind in eachindex(x)
        x[ind] = inf(I[ind]) + TaylorN(ind, order=order)*diam(I[ind])
    end
    return a(x)
end

end