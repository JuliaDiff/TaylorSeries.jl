module TaylorSeriesIAExt

using TaylorSeries

import Base: ^, sqrt, log, asin, acos, acosh, atanh, iszero, ==

import TaylorSeries: evaluate, _evaluate, normalize_taylor

isdefined(Base, :get_extension) ? (using IntervalArithmetic) : (using ..IntervalArithmetic)

const NumTypes = IntervalArithmetic.NumTypes

"""
    intersect_interval_nonstd(x, y)

Returns the intersection of the intervals `x` and `y`, considered as (extended)
sets of real numbers. That is, the set that contains the points common in `x`
and `y`.

This function is similar to [`intersect_interval`](@ref), but allows for higher
than `trv` decorations depending on `issubset_interval(x, y)`: if
`issubset_interval(x,y) == true` the decoration is that of the intersection
of the bare intervals, otherwise it is `trv` (Section 11.7.1).
"""
_intersect_domain_nonstd(x::BareInterval, y::BareInterval) = intersect_interval(x, y)

function _intersect_domain_nonstd(x::Interval{T}, y::Interval{T}) where {T<:NumTypes}
    d = ifelse(issubset_interval(x, y), decoration(x), trv)
    x = intersect_interval(x, y)
    return IntervalArithmetic._unsafe_interval(bareinterval(x), d, isguaranteed(x))
end

_intersect_domain_nonstd(x::Interval, y::Interval) =
    _intersect_domain_nonstd(promote(x, y)...)

# Some functions require special interval functions (isequal_interval, isthinzero)
for I in (:Interval, :ComplexI)
    @eval begin
        TS._isthinzero(x::$I{T}) where {T<:Real} = isthinzero(x)

        function ==(a::Taylor1{$I{T}}, b::Taylor1{$I{S}}) where {T<:NumTypes, S<:NumTypes}
            if a.order != b.order
                a, b = fixorder(a, b)
            end
            return all(isequal_interval.(a.coeffs, b.coeffs))
        end

        function ==(a::HomogeneousPolynomial{$I{T}},
                b::HomogeneousPolynomial{$I{S}}) where {T<:NumTypes, S<:NumTypes}
            a.order == b.order && return all(isequal_interval.(a.coeffs, b.coeffs))
            return all(TS._isthinzero.(a.coeffs)) && all(TS._isthinzero.(b.coeffs))
        end

        iszero(a::Taylor1{$I{T}}) where {T<:NumTypes} = all(TS._isthinzero.(a.coeffs))

        iszero(a::HomogeneousPolynomial{$I{T}}) where {T<:NumTypes} =
            all(TS._isthinzero.(a.coeffs))
    end

# Methods related to power, sqr, sqrt, ...
    for T in (:Taylor1, :TaylorN)
        @eval begin
            function ^(a::$T{$I{T}}, n::S) where {T<:NumTypes, S<:Integer}
                n == 0 && return one(a)
                n == 1 && return copy(a)
                n == 2 && return TS.square(a)
                n < 0 && return a^float(n)
                return Base.power_by_squaring(a, n)
            end

            ^(a::$T{$I{T}}, r::Rational) where {T<:NumTypes} = a^float(r)
        end
    end
end

function ^(a::Taylor1{Interval{T}}, r::S) where {T<:NumTypes, S<:Real}
    a0 = intersect_interval(constant_term(a), interval(zero(T), T(Inf)))
    aux = one(a0)^r
    iszero(r) && return Taylor1(aux, a.order)
    # aa = one(aux) * a
    aa = convert(Taylor1{typeof(aux)}, a)
    aa[0] = one(aux) * a0
    r == 1 && return aa
    r == 2 && return TS.square(aa)
    r == 1/2 && return sqrt(aa)
    l0 = findfirst(a)
    lnull = trunc(Int, r*l0 )
    if (a.order-lnull < 0) || (lnull > a.order)
        return Taylor1( zero(aux), a.order)
    end
    c_order = l0 == 0 ? a.order : min(a.order, trunc(Int,r*a.order))
    c = Taylor1(zero(aux), c_order)
    for k in eachindex(c)
        TS.pow!(c, aa, c, r, k)
    end
    return c
end

function ^(a::TaylorN{Interval{T}}, r::S) where {T<:NumTypes, S<:Real}
    isinteger(r) && return a^Int(r)
    a0 = intersect_interval(constant_term(a), interval(zero(T), T(Inf)))
    a0r = a0^r
    aux = one(a0r)
    iszero(r) && return TaylorN(aux, a.order)
    aa = aux * a
    aa[0] = aux * a0
    r == 1 && return aa
    r == 2 && return TS.square(aa)
    r == 1/2 && return sqrt(aa)
    TS._isthinzero(a0) && throw(DomainError(a,
        """The 0-th order TaylorN coefficient must be non-zero
        in order to expand `^` around 0."""))
    c = TaylorN( a0r, a.order)
    for ord in 1:a.order
        TS.pow!(c, aa, c, r, ord)
    end
    return c
end


# sqr!
for T = (:Taylor1, :TaylorN)
    @eval begin
        @inline function TS.sqr!(c::$T{Interval{T}}, a::$T{Interval{T}},
                k::Int) where {T<:NumTypes}
            if k == 0
                TS.sqr_orderzero!(c, a)
                return nothing
            end
            # Sanity
            TS.zero!(c, k)
            # Recursion formula
            kodd = k%2
            kend = div(k - 2 + kodd, 2)
            if $T == Taylor1
                @inbounds for i = 0:kend
                    c[k] += a[i] * a[k-i]
                end
            else
                @inbounds for i = 0:kend
                    TS.mul!(c[k], a[i], a[k-i])
                end
            end
            @inbounds TS.mul!(c, interval(2), c, k)
            kodd == 1 && return nothing
            if $T == Taylor1
                @inbounds c[k] += a[k >> 1]^2
            else
                TS.accsqr!(c[k], a[k >> 1])
            end
            return nothing
        end

        @inline function TS.sqr!(c::$T{Interval{T}}, k::Int) where {T<:NumTypes}
            if k == 0
                TS.sqr_orderzero!(c, c)
                return nothing
            end
            # Recursion formula
            kodd = k%2
            kend = (k - 2 + kodd) >> 1
            if $T == Taylor1
                (kend ≥ 0) && ( @inbounds c[k] = c[0] * c[k] )
                @inbounds for i = 1:kend
                    c[k] += c[i] * c[k-i]
                end
                @inbounds c[k] = interval(2) * c[k]
                (kodd == 0) && ( @inbounds c[k] += c[k >> 1]^2 )
            else
                (kend ≥ 0) && ( @inbounds TS.mul!(c, c[0][1], c, k) )
                @inbounds for i = 1:kend
                    TS.mul!(c[k], c[i], c[k-i])
                end
                @inbounds TS.mul!(c, interval(2), c, k)
                if (kodd == 0)
                    TS.accsqr!(c[k], c[k >> 1])
                end
            end

            return nothing
        end
    end
end
@inline function TS.accsqr!(c::HomogeneousPolynomial{Interval{T}},
        a::HomogeneousPolynomial{Interval{T}}) where {T<:NumTypes}
    iszero(a) && return nothing
    @inbounds num_coeffs_a = TS.size_table[a.order+1]
    @inbounds posTb = TS.pos_table[c.order+1]
    @inbounds idxTb = TS.index_table[a.order+1]
    @inbounds for na = 1:num_coeffs_a
        ca = a[na]
        TS._isthinzero(ca) && continue
        inda = idxTb[na]
        pos = posTb[2*inda]
        c[pos] += ca^2
        @inbounds for nb = na+1:num_coeffs_a
            cb = a[nb]
            TS._isthinzero(cb) && continue
            indb = idxTb[nb]
            pos = posTb[inda+indb]
            c[pos] += interval(2) * ca * cb
        end
    end
    return nothing
end


function sqrt(a::Taylor1{Interval{T}}) where {T<:NumTypes}
    domain = interval(zero(T), typemax(T))
    a0 = _intersect_domain_nonstd(constant_term(a), domain)
    aux = sqrt(a0)
    isempty_interval(aux) && throw(DomainError(a,
        """The 0-th order coefficient must have a positive part
        in order to expand `sqrt`."""))
    # First non-zero coefficient
    aa = convert(Taylor1{typeof(aux)}, a)
    aa[0] = one(aux)*a0
    l0nz = findfirst(aa)
    order = a.order
    if l0nz < 0
        return Taylor1(zero(aux), order)
    elseif l0nz%2 == 1 # l0nz must be pair
        throw(DomainError(aa,
        """First non-vanishing Taylor1 coefficient must correspond
        to a **even power** in order to expand `sqrt`."""))
    end
    # The last l0nz coefficients are set to zero.
    lnull = l0nz >> 1 # integer division by 2
    c_order = l0nz == 0 ? order : order >> 1
    c = Taylor1( zero(aux), c_order )
    @inbounds c[lnull] = aux
    for k = lnull+1:c_order
        TS.sqrt!(c, aa, k, lnull)
    end
    return c
end

function sqrt(a::TaylorN{Interval{T}}) where {T<:NumTypes}
    domain = interval(zero(T), typemax(T))
    a0 = _intersect_domain_nonstd(constant_term(a), domain)
    aux = sqrt(a0)
    (isempty_interval(aux) || TS._isthinzero(a0)) && throw(DomainError(a,
        """The 0-th order coefficient must have a positive part
        in order to expand `sqrt`."""))
    # First non-zero coefficient
    aa = convert(TaylorN{typeof(aux)}, a)
    aa[0] = one(aux)*a0
    order = a.order
    c = TaylorN( zero(aux), order)
    for k in 1:order
        TS.sqrt!(c, aa, k)
    end
    return c
end

@inline function TS.sqrt!(c::Taylor1{Interval{T}}, a::Taylor1{Interval{T}},
        k::Int, k0::Int=0) where {T<:AbstractFloat}
    if k == k0
        @inbounds c[k] = sqrt(a[2*k0])
        return nothing
    end
    kodd = (k - k0)%2
    kend = div(k - k0 - 2 + kodd, 2)
    imax = min(k0+kend, a.order)
    imin = max(k0+1, k+k0-a.order)
    imin ≤ imax && ( @inbounds c[k] = c[imin] * c[k+k0-imin] )
    @inbounds for i = imin+1:imax
        c[k] += c[i] * c[k+k0-i]
    end
    if k+k0 ≤ a.order
        @inbounds aux = a[k+k0] - interval(2) * c[k]
    else
        @inbounds aux = - interval(2) * c[k]
    end
    if kodd == 0
        @inbounds aux = aux - c[kend+k0+1]^2
    end
    @inbounds c[k] = aux / (interval(2) * c[k0])
    return nothing
end
# Faltan métodos de sqrt! para TaylorN y Taylor1{TaylorN}

# several math functions
for T in (:Taylor1, :TaylorN)
    @eval begin
        function log(a::$T{Interval{T}}) where {T<:NumTypes}
            domain = interval(zero(T), typemax(T))
            a0 = _intersect_domain_nonstd(constant_term(a), domain)
            aux = log(a0)
            isempty_interval(aux) && throw(DomainError(a,
                """The 0-th order coefficient must be positive in order to expand `log`."""))
            aa = convert($T{typeof(aux)}, a)
            aa[0] = one(aux) * a0
            order = a.order
            c = $T( aux, order )
            for k in eachindex(a)
                TS.log!(c, aa, k)
            end
            return c
        end

        function asin(a::$T{Interval{T}}) where {T<:NumTypes}
            domain = interval(-one(T), one(T))
            a0 = _intersect_domain_nonstd(constant_term(a), domain)
            aux = asin(a0)
            isempty_interval(aux) && throw(DomainError(a,
                """The 0-th order coefficient must have a non-empty intersection with $domain."""))
            a0sqr = a0^2
            uno = one(aux)
            isequal_interval(a0sqr, uno) && throw(DomainError(a,
                    """Series expansion of asin(x) diverges at x = ±1."""))
            order = a.order
            aa = convert($T{typeof(aux)}, a)
            aa[0] = uno * a0
            c = $T( aux, order )
            r = $T( sqrt(uno - a0sqr), order )
            for k in eachindex(a)
                TS.asin!(c, aa, r, k)
            end
            return c
        end

        function acos(a::$T{Interval{T}}) where {T<:NumTypes}
            domain = interval(-one(T), one(T))
            a0 = _intersect_domain_nonstd(constant_term(a), domain)
            aux = acos(a0)
            isempty_interval(aux) && throw(DomainError(a,
                """The 0-th order coefficient must have a non-empty intersection with $domain."""))
            a0sqr = a0^2
            uno = one(a0)
            isequal_interval(a0sqr, uno) && throw(DomainError(a,
                    """Series expansion of acos(x) diverges at x = ±1."""))
            order = a.order
            aa = convert($T{typeof(aux)}, a)
            aa[0] = uno * a0
            c = $T( aux, order )
            r = $T( sqrt(uno - a0sqr), order )
            for k in eachindex(a)
                TS.acos!(c, aa, r, k)
            end
            return c
        end

        function acosh(a::$T{Interval{T}}) where {T<:NumTypes}
            domain = interval(one(T), typemax(T))
            a0 = _intersect_domain_nonstd(constant_term(a), domain)
            aux = acosh(a0)
            isempty_interval(aux) && throw(DomainError(a,
                """The 0-th order coefficient must have a non-empty intersection with $domain."""))
            a0sqr = a0^2
            uno = one(a0)
            isequal_interval(a0sqr, uno) && throw(DomainError(a,
                """Series expansion of acosh(x) diverges at x = ±1."""))
            order = a.order
            aa = convert($T{typeof(aux)}, a)
            aa[0] = uno * a0
            c = $T( aux, order )
            r = $T( sqrt(a0sqr - uno), order )
            for k in eachindex(a)
                TS.acosh!(c, aa, r, k)
            end
            return c
        end

        function atanh(a::$T{Interval{T}}) where {T<:NumTypes}
            domain = interval(-one(T), one(T))
            a0 = _intersect_domain_nonstd(constant_term(a), domain)
            aux = atanh(a0)
            isempty_interval(aux) && throw(DomainError(a,
                """The 0-th order coefficient must have a non-empty intersection with $domain."""))
            order = a.order
            uno = one(a0)
            aa = convert($T{typeof(aux)}, a)
            aa[0] = uno * a0
            c = $T( aux, order)
            r = $T(uno - a0^2, order)
            TS._isthinzero(constant_term(r)) && throw(DomainError(a,
                """Series expansion of atanh(x) diverges at x = ±1."""))
            for k in eachindex(a)
                TS.atanh!(c, aa, r, k)
            end
            return c
        end

        # Some internal functions
        @inline function TS.exp!(c::$T{Interval{T}}, a::$T{Interval{T}}, k::Int) where {T<:AbstractFloat}
            if k == 0
                @inbounds c[0] = exp(constant_term(a))
                return nothing
            end
            if $T == Taylor1
                @inbounds c[k] = interval(k) * a[k] * c[0]
            else
                @inbounds TS.mul!(c[k], interval(k) * a[k], c[0])
            end
            @inbounds for i = 1:k-1
                if $T == Taylor1
                    c[k] += interval(k-i) * a[k-i] * c[i]
                else
                    TS.mul!(c[k], interval(k-i) * a[k-i], c[i])
                end
            end
            @inbounds c[k] = c[k] / interval(k)
            return nothing
        end

        @inline function TS.expm1!(c::$T{Interval{T}}, a::$T{Interval{T}}, k::Int) where {T<:AbstractFloat}
            if k == 0
                @inbounds c[0] = expm1(constant_term(a))
                return nothing
            end
            c0 = c[0]+one(c[0])
            if $T == Taylor1
                @inbounds c[k] = interval(k) * a[k] * c0
            else
                @inbounds TS.mul!(c[k], interval(k) * a[k], c0)
            end
            @inbounds for i = 1:k-1
                if $T == Taylor1
                    c[k] += interval(k-i) * a[k-i] * c[i]
                else
                    TS.mul!(c[k], interval(k-i) * a[k-i], c[i])
                end
            end
            @inbounds c[k] = c[k] / interval(k)
            return nothing
        end

        @inline function TS.log!(c::$T{Interval{T}}, a::$T{Interval{T}}, k::Int) where {T<:AbstractFloat}
            if k == 0
                @inbounds c[0] = log(constant_term(a))
                return nothing
            elseif k == 1
                @inbounds c[1] = a[1] / constant_term(a)
                return nothing
            end
            if $T == Taylor1
                @inbounds c[k] = interval(k-1) * a[1] * c[k-1]
            else
                @inbounds TS.mul!(c[k], interval(k-1)*a[1], c[k-1])
            end
            @inbounds for i = 2:k-1
                if $T == Taylor1
                    c[k] += interval(k-i) * a[i] * c[k-i]
                else
                    TS.mul!(c[k], interval(k-i)*a[i], c[k-i])
                end
            end
            @inbounds c[k] = (a[k] - c[k]/interval(k)) / constant_term(a)
            return nothing
        end

        @inline function TS.log1p!(c::$T{Interval{T}}, a::$T{Interval{T}}, k::Int) where {T<:AbstractFloat}
            a0 = constant_term(a)
            a0p1 = a0+one(a0)
            if k == 0
                @inbounds c[0] = log1p(a0)
                return nothing
            elseif k == 1
                @inbounds c[1] = a[1] / a0p1
                return nothing
            end
            if $T == Taylor1
                @inbounds c[k] = interval(k-1) * a[1] * c[k-1]
            else
                @inbounds TS.mul!(c[k], interval(k-1)*a[1], c[k-1])
            end
            @inbounds for i = 2:k-1
                if $T == Taylor1
                    c[k] += interval(k-i) * a[i] * c[k-i]
                else
                    TS.mul!(c[k], interval(k-i)*a[i], c[k-i])
                end
            end
            @inbounds c[k] = (a[k] - c[k]/interval(k)) / a0p1
            return nothing
        end

        @inline function TS.sincos!(s::$T{Interval{T}}, c::$T{Interval{T}}, a::$T{Interval{T}},
                k::Int) where {T<:AbstractFloat}
            if k == 0
                a0 = constant_term(a)
                @inbounds s[0], c[0] = sincos( a0 )
                return nothing
            end
            x = a[1]
            if $T == Taylor1
                @inbounds s[k] = x * c[k-1]
                @inbounds c[k] = -x * s[k-1]
            else
                TS.mul!(s[k], x, c[k-1])
                TS.mul!(c[k], -x, s[k-1])
            end
            @inbounds for i = 2:k
                x = interval(i) * a[i]
                if $T == Taylor1
                    s[k] += x * c[k-i]
                    c[k] -= x * s[k-i]
                else
                    TS.mul!(s[k], x, c[k-i])
                    TS.mul!(c[k], -x, s[k-i])
                end
            end
            @inbounds s[k] = s[k] / interval(k)
            @inbounds c[k] = c[k] / interval(k)
            return nothing
        end

        @inline function TS.tan!(c::$T{Interval{T}}, a::$T{Interval{T}}, c2::$T{Interval{T}},
                k::Int) where {T<:AbstractFloat}
            if k == 0
                @inbounds aux = tan( constant_term(a) )
                @inbounds c[0] = aux
                @inbounds c2[0] = aux^2
                return nothing
            end
            if $T == Taylor1
                @inbounds c[k] = interval(k) * a[k] * c2[0]
            else
                @inbounds TS.mul!(c[k], interval(k) * a[k], c2[0])
            end
            @inbounds for i = 1:k-1
                if $T == Taylor1
                    c[k] += interval(k-i) * a[k-i] * c2[i]
                else
                    TS.mul!(c[k], interval(k-i) * a[k-i], c2[i])
                end
            end
            @inbounds c[k] = a[k] + c[k]/interval(k)
            TS.sqr!(c2, c, k)
            return nothing
        end

        @inline function TS.asin!(c::$T{Interval{T}}, a::$T{Interval{T}}, r::$T{Interval{T}},
                k::Int) where {T<:AbstractFloat}
            if k == 0
                a0 = constant_term(a)
                @inbounds c[0] = asin( a0 )
                @inbounds r[0] = sqrt( one(a0) - a0^2 )
                return nothing
            end
            if $T == Taylor1
                @inbounds c[k] = interval(k-1) * r[1] * c[k-1]
            else
                @inbounds TS.mul!(c[k], interval(k-1) * r[1], c[k-1])
            end
            @inbounds for i in 2:k-1
                if $T == Taylor1
                    c[k] += interval(k-i) * r[i] * c[k-i]
                else
                    TS.mul!(c[k], interval(k-i) * r[i], c[k-i])
                end
            end
            TS.sqrt!(r, one(a[0])-a^2, k)
            @inbounds c[k] = (a[k] - c[k]/interval(k)) / constant_term(r)
            return nothing
        end

        @inline function TS.acos!(c::$T{Interval{T}}, a::$T{Interval{T}}, r::$T{Interval{T}},
                k::Int) where {T<:AbstractFloat}
            if k == 0
                a0 = constant_term(a)
                @inbounds c[0] = acos( a0 )
                @inbounds r[0] = sqrt( one(a0) - a0^2 )
                return nothing
            end
            if $T == Taylor1
                @inbounds c[k] = interval(k-1) * r[1] * c[k-1]
            else
                @inbounds TS.mul!(c[k], interval(k-1) * r[1], c[k-1])
            end
            @inbounds for i in 2:k-1
                if $T == Taylor1
                    c[k] += interval(k-i) * r[i] * c[k-i]
                else
                    TS.mul!(c[k], interval(k-i) * r[i], c[k-i])
                end
            end
            TS.sqrt!(r, one(a[0])-a^2, k)
            @inbounds c[k] = -(a[k] + c[k]/interval(k)) / constant_term(r)
            return nothing
        end

        @inline function TS.atan!(c::$T{Interval{T}}, a::$T{Interval{T}}, r::$T{Interval{T}},
                k::Int) where {T<:AbstractFloat}
            if k == 0
                a0 = constant_term(a)
                @inbounds c[0] = atan( a0 )
                @inbounds r[0] = one(a0) + a0^2
                return nothing
            end
            if $T == Taylor1
                @inbounds c[k] = interval(k-1) * r[1] * c[k-1]
            else
                @inbounds TS.mul!(c[k], interval(k-1) * r[1], c[k-1])
            end
            @inbounds for i in 2:k-1
                if $T == Taylor1
                    c[k] += interval(k-i) * r[i] * c[k-i]
                else
                    TS.mul!(c[k], interval(k-i) * r[i], c[k-i])
                end
            end
            @inbounds TS.sqr!(r, a, k)
            @inbounds c[k] = (a[k] - c[k]/interval(k)) / constant_term(r)
            return nothing
        end

        @inline function TS.sinhcosh!(s::$T{Interval{T}}, c::$T{Interval{T}}, a::$T{Interval{T}},
                k::Int) where {T<:AbstractFloat}
            if k == 0
                @inbounds s[0] = sinh( constant_term(a) )
                @inbounds c[0] = cosh( constant_term(a) )
                return nothing
            end
            x = a[1]
            if $T == Taylor1
                @inbounds s[k] = x * c[k-1]
                @inbounds c[k] = x * s[k-1]
            else
                @inbounds TS.mul!(s[k], x, c[k-1])
                @inbounds TS.mul!(c[k], x, s[k-1])
            end
            @inbounds for i = 2:k
                x = interval(i) * a[i]
                if $T == Taylor1
                    s[k] += x * c[k-i]
                    c[k] += x * s[k-i]
                else
                    TS.mul!(s[k], x, c[k-i])
                    TS.mul!(c[k], x, s[k-i])
                end
            end
            s[k] = s[k] / interval(k)
            c[k] = c[k] / interval(k)
            return nothing
        end

        @inline function TS.tanh!(c::$T{Interval{T}}, a::$T{Interval{T}}, c2::$T{Interval{T}},
                k::Int) where {T<:AbstractFloat}
            if k == 0
                @inbounds aux = tanh( constant_term(a) )
                @inbounds c[0] = aux
                @inbounds c2[0] = aux^2
                return nothing
            end
            if $T == Taylor1
                @inbounds c[k] = k * a[k] * c2[0]
            else
                @inbounds TS.mul!(c[k], k * a[k], c2[0])
            end
            @inbounds for i = 1:k-1
                if $T == Taylor1
                    c[k] += (k-i) * a[k-i] * c2[i]
                else
                    TS.mul!(c[k], (k-i) * a[k-i], c2[i])
                end
            end
            @inbounds c[k] = a[k] - c[k]/k
            TS.sqr!(c2, c, k)
            return nothing
        end

        @inline function TS.asinh!(c::$T{Interval{T}}, a::$T{Interval{T}}, r::$T{Interval{T}},
                k::Int) where {T<:AbstractFloat}
            if k == 0
                a0 = constant_term(a)
                @inbounds c[0] = asinh( a0 )
                @inbounds r[0] = sqrt( a0^2 + one(a0) )
                return nothing
            end
            if $T == Taylor1
                @inbounds c[k] = interval(k-1) * r[1] * c[k-1]
            else
                @inbounds TS.mul!(c[k], interval(k-1) * r[1], c[k-1])
            end
            @inbounds for i in 2:k-1
                if $T == Taylor1
                    c[k] += interval(k-i) * r[i] * c[k-i]
                else
                    TS.mul!(c[k], interval(k-i) * r[i], c[k-i])
                end
            end
            TS.sqrt!(r, a^2+one(a[0]), k)
            @inbounds c[k] = (a[k] - c[k]/interval(k)) / constant_term(r)
            return nothing
        end

        @inline function TS.acosh!(c::$T{Interval{T}}, a::$T{Interval{T}}, r::$T{Interval{T}},
                k::Int) where {T<:AbstractFloat}
            if k == 0
                a0 = constant_term(a)
                @inbounds c[0] = acosh( a0 )
                @inbounds r[0] = sqrt( a0^2 - one(a0) )
                return nothing
            end
            if $T == Taylor1
                @inbounds c[k] = interval(k-1) * r[1] * c[k-1]
            else
                @inbounds TS.mul!(c[k], interval(k-1) * r[1], c[k-1])
            end
            @inbounds for i in 2:k-1
                if $T == Taylor1
                    c[k] += interval(k-i) * r[i] * c[k-i]
                else
                    TS.mul!(c[k], interval(k-i) * r[i], c[k-i])
                end
            end
            TS.sqrt!(r, a^2-one(a[0]), k)
            @inbounds c[k] = (a[k] - c[k]/interval(k)) / constant_term(r)
            return nothing
        end

        @inline function TS.atanh!(c::$T{Interval{T}}, a::$T{Interval{T}}, r::$T{Interval{T}},
                k::Int) where {T<:AbstractFloat}
            if k == 0
                a0 = constant_term(a)
                @inbounds c[0] = atanh( a0 )
                @inbounds r[0] = one(a0) - a0^2
                return nothing
            end
            if $T == Taylor1
                @inbounds c[k] = interval(k-1) * r[1] * c[k-1]
            else
                @inbounds TS.mul!(c[k], interval(k-1) * r[1], c[k-1])
            end
            @inbounds for i in 2:k-1
                if $T == Taylor1
                    c[k] += interval(k-i) * r[i] * c[k-i]
                else
                    TS.mul!(c[k], interval(k-i) * r[i], c[k-i])
                end
            end
            @inbounds TS.sqr!(r, a, k)
            @inbounds c[k] = (a[k] + c[k]/interval(k)) / constant_term(r)
            return nothing
        end
    end
end


function evaluate(a::Taylor1{T}, dx::Interval{S}) where {T<:Real, S<:NumTypes}
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

function evaluate(a::TaylorN{T}, dx::AbstractVector{Interval{S}}) where
        {T<:Real, S<:NumTypes}
    @assert length(dx) == get_numvars()
    suma = zero(constant_term(a)) + interval(zero(T))
    @inbounds for homPol in reverse(eachindex(a))
        suma += evaluate(a.coeffs[homPol+1], dx)
    end
    return suma
end

function evaluate(a::Taylor1{TaylorN{T}}, dx::Interval{S}) where {T<:Real, S<:Real}
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


function evaluate(a::HomogeneousPolynomial{T},
        dx::AbstractVector{Interval{S}}) where {T<:Real, S<:NumTypes}
    @assert length(dx) == get_numvars()
    all(isequal_interval.(dx, interval(-one(T), one(T)))) &&
        return TS._evaluate(a, dx, Val(true))
    all(isequal_interval.(dx, interval(zero(T), one(T)))) &&
        return TS._evaluate(a, dx, Val(false))
    return TS._evaluate(a, dx)
end

evaluate(a::TaylorN{Taylor1{T}}, vals::AbstractVector{Interval{S}}) where
    {T<:Real, S<:NumTypes} = TS._evaluate(a, (vals...,), Val(false))

# _evaluate
function TS._evaluate(a::HomogeneousPolynomial{T},
        dx::AbstractVector{Interval{S}}) where {T<:Real, S<:NumTypes}
    a.order == 0 && return a[1] + interval(zero(T))
    ct = TS.coeff_table[a.order+1]
    @inbounds suma = a[1]*interval(zero(T))
    for (i, a_coeff) in enumerate(a.coeffs)
        TS._isthinzero(a_coeff) && continue
        @inbounds tmp = prod(dx .^ ct[i])
        suma += a_coeff * tmp
    end
    return suma
end

function TS._evaluate(a::HomogeneousPolynomial{T}, dx::AbstractVector{Interval{S}},
        ::Val{true} ) where {T<:Real, S<:NumTypes}
    a.order == 0 && return a[1] + interval(zero(T))
    ct = TS.coeff_table[a.order+1]
    @inbounds suma = a[1]*interval(zero(T))
    Ieven = interval(zero(T), one(T))
    for (i, a_coeff) in enumerate(a.coeffs)
        TS._isthinzero(a_coeff) && continue
        if isodd(sum(ct[i]))
            suma += sum(a_coeff) * dx[1]
            continue
        end
        @inbounds tmp = iseven(ct[i][1]) ? Ieven : dx[1]
        for n in 2:length(dx)
            @inbounds vv = iseven(ct[i][n]) ? Ieven : dx[1]
            tmp *= vv
        end
        suma += a_coeff * tmp
    end
    return suma
end

function TS._evaluate(a::HomogeneousPolynomial{T}, dx::AbstractVector{Interval{S}},
        ::Val{false} ) where {T<:Real, S<:NumTypes}
    a.order == 0 && return a[1] + interval(zero(T))
    @inbounds suma = zero(a[1])*dx[1]
    @inbounds for homPol in a.coeffs
        suma += homPol*dx[1]
    end
    return suma
end

function TS._evaluate(a::TaylorN{T}, vals::NTuple{N,TaylorN{Interval{S}}}) where
        {N, T<:Real, S<:NumTypes}
    @assert get_numvars() == N
    R = promote_type(TS.numtype(a), typeof(vals[1]))
    a_length = length(a)
    suma = zeros(R, a_length)
    @inbounds for homPol in 1:a_length
        suma[homPol] = TS._evaluate(a.coeffs[homPol], vals)
    end
    return suma
end

function TS._evaluate(a::HomogeneousPolynomial{T},
        vals::NTuple{N,TaylorN{Interval{S}}}) where {N, T<:Real, S<:NumTypes}
    ct = TS.coeff_table[a.order+1]
    suma = zero(a[1])*vals[1]
    for (i, a_coeff) in enumerate(a.coeffs)
        TS._isthinzero(a_coeff) && continue
        @inbounds tmp = prod( vals .^ ct[i] )
        suma += a_coeff * tmp
    end
    return suma
end


"""
    normalize_taylor(a::Taylor1, I::Interval, symI::Bool=true)

Normalizes `a::Taylor1` such that the interval `I` is mapped
by an affine transformation to the interval `-1..1` (`symI=true`)
or to `0..1` (`symI=false`).
"""
normalize_taylor(a::Taylor1, I::Interval{T}, symI::Bool=true) where {T<:NumTypes} =
    _normalize(a, I, Val(symI))

"""
    normalize_taylor(a::TaylorN, I::AbstractVector{Interval{T}}, symI::Bool=true)

Normalize `a::TaylorN` such that the intervals in `I::AbstractVector{Interval{T}}`
are mapped by an affine transformation to the intervals `-1..1`
(`symI=true`) or to `0..1` (`symI=false`).
"""
normalize_taylor(a::TaylorN, I::AbstractVector{Interval{T}},
    symI::Bool=true) where {T<:NumTypes} = _normalize(a, I, Val(symI))

#  I -> -1..1
function _normalize(a::Taylor1, I::Interval{T}, ::Val{true}) where {T<:NumTypes}
    order = get_order(a)
    t = Taylor1(TS.numtype(a), order)
    tnew = mid(I) + t*radius(I)
    return a(tnew)
end

function _normalize(a::Taylor1{Interval{T}}, I::Interval{T}, ::Val{true}) where
        {T<:NumTypes}
    order = get_order(a)
    t = Taylor1(TS.numtype(a), order)
    tnew = interval(mid(I)) + t*interval(radius(I))
    return a(tnew)
end

#  I -> 0..1
function _normalize(a::Taylor1, I::Interval{T}, ::Val{false}) where {T<:NumTypes}
    order = get_order(a)
    t = Taylor1(TS.numtype(a), order)
    tnew = inf(I) + t*diam(I)
    return a(tnew)
end

function _normalize(a::Taylor1{Interval{T}}, I::Interval{T}, ::Val{false}) where
        {T<:NumTypes}
    order = get_order(a)
    t = Taylor1(TS.numtype(a), order)
    tnew = interval(inf(I)) + t*interval(diam(I))
    return a(tnew)
end

#  I -> IntervalBox(-1..1, Val(N))
function _normalize(a::TaylorN, I::AbstractVector{Interval{T}},
        ::Val{true}) where {T<:NumTypes}
    order = get_order(a)
    x = Vector{typeof(a)}(undef, length(I))
    @inbounds for ind in eachindex(x)
        x[ind] = mid(I[ind]) + TaylorN(ind, order=order)*radius(I[ind])
    end
    return a(x)
end

function _normalize(a::TaylorN{Interval{T}}, I::AbstractVector{Interval{T}},
        ::Val{true}) where {T<:NumTypes}
    order = get_order(a)
    x = Vector{typeof(a)}(undef, length(I))
    @inbounds for ind in eachindex(x)
        x[ind] = interval(mid(I[ind])) +
            TaylorN(ind, order=order)*interval(radius(I[ind]))
    end
    return a(x)
end

#  I -> IntervalBox(0..1, Val(N))
function _normalize(a::TaylorN, I::AbstractVector{Interval{T}},
        ::Val{false}) where {T<:NumTypes}
    order = get_order(a)
    x = Vector{typeof(a)}(undef, length(I))
    @inbounds for ind in eachindex(x)
        x[ind] = inf(I[ind]) + TaylorN(ind, order=order)*diam(I[ind])
    end
    return a(x)
end

function _normalize(a::TaylorN{Interval{T}}, I::AbstractVector{Interval{T}},
        ::Val{false}) where {T<:NumTypes}
    order = get_order(a)
    x = Vector{typeof(a)}(undef, length(I))
    @inbounds for ind in eachindex(x)
        x[ind] = interval(inf(I[ind])) +
            TaylorN(ind, order=order)*interval(diam(I[ind]))
    end
    return a(x)
end

# Printing-related methods numbr2str
function TS.numbr2str(zz::Interval, ifirst::Bool=false)
    TS._isthinzero(zz) && return string( zz )
    plusmin = ifelse( ifirst, string(""), string("+ ") )
    return string(plusmin, zz)
end

function TS.numbr2str(zz::ComplexI, ifirst::Bool=false)
    zT = zero(zz.re)
    TS._isthinzero(zz) && return string(zT)
    if ifirst
        cadena = string("( ", zz, " )")
    else
        cadena = string("+ ( ", zz, " )")
    end
    return cadena
end

end
