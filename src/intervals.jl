using .IntervalArithmetic

# Method used for Taylor1{Interval{T}}^n
function ^(a::Taylor1{T}, n::Integer) where {T<:Interval}
    n == 0 && return one(a)
    n == 1 && return copy(a)
    n == 2 && return square(a)
    n < 0 && return a^float(n)
    return power_by_squaring(a, n)
end

function ^(a::Taylor1{Interval{T}}, r::S) where {T<:Real, S<:Real}
    a0 = constant_term(a) ∩ Interval(zero(T), T(Inf))
    aux = one(a0)^r

    iszero(r) && return Taylor1(aux, a.order)
    aa = one(aux) * a
    aa[0] = one(aux) * a0
    r == 1 && return aa
    r == 2 && return square(aa)
    r == 1/2 && return sqrt(aa)

    l0 = findfirst(a)
    lnull = trunc(Int, r*l0 )
    if (a.order-lnull < 0) || (lnull > a.order)
        return Taylor1( zero(aux), a.order)
    end
    c_order = l0 == 0 ? a.order : min(a.order, trunc(Int,r*a.order))
    c = Taylor1(zero(aux), c_order)
    for k = 0:c_order
        pow!(c, aa, r, k)
    end

    return c
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
            log!(c, aa, k)
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
            asin!(c, aa, r, k)
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
            acos!(c, aa, r, k)
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
            acosh!(c, aa, r, k)
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
            atanh!(c, aa, r, k)
        end
        return c
    end
end

function evaluate(a::Taylor1, dx::Interval)
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
    @inbounds for homPol in length(a):-1:1
        suma += evaluate(a.coeffs[homPol], dx)
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

    ct = coeff_table[a.order+1]
    @inbounds suma = a[1]*Interval{T}(0,0)

    Ieven = Interval{T}(0,1)
    for (i,a_coeff) in enumerate(a.coeffs)
        iszero(a_coeff) && continue
        if isodd(sum(ct[i]))
            suma += sum(a_coeff) * dx[1]
            continue
        end
        @inbounds tmp = iseven(ct[i][1]) ? Ieven : dx[1]
        for n in 2:N
            @inbounds vv = iseven(ct[i][n]) ? Ieven : dx[1]
            tmp *= vv
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
