module TaylorSeriesIAExt

using TaylorSeries

import Base: ^, sqrt, log, asin, acos, acosh, atanh, iszero, ==,
        power_by_squaring

import TaylorSeries: _pow, evaluate, #_evaluate, #_evaluate!,
        normalize_taylor, aff_normalize

using IntervalArithmetic

const NumTypes = IntervalArithmetic.NumTypes


# Internal; union type
const _IntervalVals{S} =
    Union{AbstractVector{Interval{S}}, Tuple{Interval{S}, Vararg{Interval{S}}}}


# _defaultsorting -> false for Interval evaluations
TS._defaultsorting(::Type{T}, ::Type{S}) where {T<:Interval, S} = false
TS._defaultsorting(::Type{T}, ::Type{S}) where {T, S<:Interval} = false
TS._defaultsorting(::Type{T}, ::Type{S}) where {T<:Interval, S<:Interval} = false


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
            if order(a) != order(b)
                a, b = TS.fixorder(a, b)
            end
            return all(isequal_interval.(a.coeffs, b.coeffs))
        end

        function ==(a::HomogeneousPolynomial{$I{T}},
                b::HomogeneousPolynomial{$I{S}}) where {T<:NumTypes, S<:NumTypes}
            order(a) == order(b) &&
                return all(isequal_interval.(a.coeffs, b.coeffs))
            return all(TS._isthinzero, a.coeffs) && all(TS._isthinzero, b.coeffs)
        end

        iszero(a::Taylor1{$I{T}}) where {T<:NumTypes} = all(TS._isthinzero, a.coeffs)

        iszero(a::HomogeneousPolynomial{$I{T}}) where {T<:NumTypes} =
            all(TS._isthinzero, a.coeffs)
    end
end

# Methods related to power, sqr, sqrt, ...
for T in (:Taylor1, :TaylorN)
    @eval begin
        function ^(a::$T{Interval{T}}, n::S) where {T<:NumTypes, S<:Integer}
            n == 0 && return one(a)
            n == 1 && return TS._copy_series(a)
            n == 2 && return TS.square(a)
            n < 0 && return a^float(n)
            return power_by_squaring(a, n)
        end

        ^(a::$T{Interval{T}}, r::Rational) where {T<:NumTypes} = a^float(r)

        # _pow
        function _pow(a::$T{Interval{S}}, n::Integer) where {S<:NumTypes}
            n < 0 && return _pow(a, float(n))
            return power_by_squaring(a, n)
        end
    end
end

function ^(a::Taylor1{Interval{T}}, r::S) where {T<:NumTypes, S<:Real}
    isinteger(r) && r >= 0 && return power_by_squaring(a, Integer(r))
    a0 = _intersect_domain_nonstd(constant_term(a), interval(zero(T), T(Inf)))
    @assert !isempty_interval(a0)
    aux = one(a0^r)
    aa = one(aux) * a
    r == 0.5 && return sqrt(aa)
    return _pow(aa, r)
end

function ^(a::TaylorN{Interval{T}}, r::S) where {T<:NumTypes, S<:Real}
    isinteger(r) && r >= 0 && return power_by_squaring(a, Integer(r))
    a0 = _intersect_domain_nonstd(constant_term(a), interval(zero(T), T(Inf)))
    @assert !isempty_interval(a0)
    aux = one(a0^r)
    aa = one(aux) * a
    r == 0.5 && return sqrt(aa)
    if TS._isthinzero(a0)
        throw(DomainError(aa,
        """The 0-th order TaylorN coefficient must be non-zero
        in order to expand `^` around 0."""))
    end
    return _pow(aa, r)
end

function _pow(a::Taylor1{Interval{T}}, r::S) where {T<:NumTypes, S<:Real}
    isinteger(r) && r >= 0 && return power_by_squaring(a, Integer(r))
    a0 = _intersect_domain_nonstd(constant_term(a), interval(zero(T), T(Inf)))
    @assert !isempty_interval(a0)
    aux = one(a0^r)
    aa = TS._copy_series(a)
    aa[0] = aux * a0
    r == 0.5 && return sqrt(aa)
    a_order = order(aa)
    l0 = findfirst(aa)
    r == 0.5 && return sqrt(a)
    a_order = order(a)
    l0 = findfirst(a)
    # Index of first non-zero coefficient of the result; must be integer
    !isinteger(r*l0) && throw(DomainError(aa,
        """The 0-th order Taylor1 coefficient must be non-zero
        to raise the Taylor1 polynomial to a non-integer exponent."""))
    lnull = trunc(Int, r*l0 )
    (lnull > a_order) && return Taylor1( zero(aux), a_order)
    c_order = l0 == 0 ? a_order : min(a_order, trunc(Int, r*a_order))
    c = Taylor1(zero(aux), c_order)
    aux0 = zero(c)
    for k in eachindex(c)
        TS.pow!(c, aa, aux0, r, k)
    end
    return c
end

function _pow(a::TaylorN{Interval{T}}, r::S) where {T<:NumTypes, S<:Real}
    isinteger(r) && r >= 0 && return power_by_squaring(a, Integer(r))
    a0 = _intersect_domain_nonstd(constant_term(a), interval(zero(T), T(Inf)))
    @assert !isempty_interval(a0)
    aux = one(a0^r)
    # work on a copy, so the caller's series is not modified
    aa = TS._copy_series(a)
    aa[0] = aux * a0
    r == 0.5 && return sqrt(aa)
    a_order = order(aa)
    if TS._isthinzero(a0)
        throw(DomainError(aa,
            """The 0-th order TaylorN coefficient must be non-zero
            in order to expand `^` around 0."""))
    end
    c = TaylorN(space(aa), zero(aux), a_order)
    aux0 = zero(c)
    for k in eachindex(c)
        TS.pow!(c, aa, aux0, r, k)
    end
    return c
end

# sqr!
function TS.sqr!(c::Taylor1{Interval{T}}, a::Taylor1{Interval{T}},
        ::Interval{T}, k::Int) where {T<:NumTypes}
    if k == 0
        TS.sqr_orderzero!(c, a)
        return nothing
    end
    # Sanity
    TS.zero!(c, k)
    # Recursion formula
    kodd = k%2
    kend = (k - 2 + kodd) >> 1
    @inbounds for i = 0:kend
        c[k] += a[i] * a[k-i]
    end
    @inbounds TS.mul!(c, interval(T(2)), c, k)
    kodd == 1 && return nothing
    @inbounds c[k] += a[k >> 1]^2
    return nothing
end

function TS.sqr!(c::TaylorN{Interval{T}}, a::TaylorN{Interval{T}},
        ::Interval{T}, k::Int) where {T<:NumTypes}
    TS._check_same_space(c, a)
    if k == 0
        TS.sqr_orderzero!(c, a)
        return nothing
    end
    # Sanity
    TS.zero!(c, k)
    # Recursion formula
    kodd = k%2
    kend = (k - 2 + kodd) >> 1
    @inbounds for i = 0:kend
        TS.mul!(c[k], a[i], a[k-i])
    end
    @inbounds TS.mul!(c, interval(T(2)), c, k)
    kodd == 1 && return nothing
    TS.accsqr!(c[k], a[k >> 1])
    return nothing
end

function TS.sqr!(c::Taylor1{Interval{T}}, k::Int) where {T<:NumTypes}
    if k == 0
        TS.sqr_orderzero!(c, c)
        return nothing
    end
    # Recursion formula
    kodd = k%2
    kend = (k - 2 + kodd) >> 1
    (kend >= 0) && ( @inbounds c[k] = c[0] * c[k] )
    @inbounds for i = 1:kend
        c[k] += c[i] * c[k-i]
    end
    @inbounds c[k] = interval(T(2)) * c[k]
    (kodd == 0) && ( @inbounds c[k] += c[k >> 1]^2 )
    return nothing
end

function TS.sqr!(c::TaylorN{Interval{T}}, k::Int) where {T<:NumTypes}
    if k == 0
        TS.sqr_orderzero!(c, c)
        return nothing
    end
    # Recursion formula
    kodd = k%2
    kend = (k - 2 + kodd) >> 1
    (kend >= 0) && ( @inbounds TS.mul!(c, c[0][1], c, k) )
    @inbounds for i = 1:kend
        TS.mul!(c[k], c[i], c[k-i])
    end
    @inbounds TS.mul!(c, interval(T(2)), c, k)
    if (kodd == 0)
        TS.accsqr!(c[k], c[k >> 1])
    end
    return nothing
end

function TS.accsqr!(c::HomogeneousPolynomial{Interval{T}},
        a::HomogeneousPolynomial{Interval{T}}) where {T<:NumTypes}
    iszero(a) && return nothing
    TS._check_same_space(c, a)
    sp = a.space
    @inbounds num_coeffs_a = sp.size_table[order(a)+1]
    @inbounds posTb = sp.pos_table[order(c)+1]
    @inbounds idxTb = sp.index_table[order(a)+1]
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
            c[pos] += interval(T(2)) * ca * cb
        end
    end
    return nothing
end


# ----------------------------------------------------------------------------
# Midpoint-radius products of `HomogeneousPolynomial{Interval}` / `TaylorN{Interval}`
#
# The generic kernels multiply and add `Interval`s, whose directed rounding is emulated
# in software (`RoundingEmulator`) and dominates the cost when the coefficients are
# (almost) thin. Here, each factor `x` is enclosed as `x ⊆ [m-r, m+r]` (hardware
# arithmetic), and the accumulated coefficient `Σ_k x_k y_k` is computed as
#
#     M ± rad,   rad = (E + R) (1 + 4(n+4)u)  [+ 16 n η if underflow is possible],
#
# where, for the `n` nonzero terms,
#   * `M` is the floating-point value of `Σ m_k m'_k` (recursive summation) and
#     `E = Σ (|e_k| + |t_k|)`, with `e_k` the exact error of each product (`fma`) and
#     `t_k` the exact error of each addition (TwoSum, Knuth): error-free transformations, so
#     `Σ m_k m'_k = M + Σ (e_k + t_k)` exactly. Hence `E` bounds the rounding error without
#     any `γ_n` estimate, and it is zero when the operations are exact;
#   * `R = Σ (r_k |m'_k| + |m_k| s_k + r_k s_k)` bounds the effect of the radii
#     (`m_k, r_k` and `m'_k, s_k` midpoints and radii of the factors):
#     `|Σ x_k y_k - Σ m_k m'_k| ≤ R` for all `x_k, y_k` in the intervals;
#   * `u = eps/2`, `η = nextfloat(0)`; the factor `(1 + 4(n+4)u)` absorbs the rounding of the
#     sums `E`, `R` of nonnegative terms and of the final product.
# The endpoints are `M ∓ rad` rounded outwards only if the subtraction/addition is inexact
# (TwoSum), so the result is exactly thin (`[M, M]`) when `rad = 0`. Terms with a thin
# zero factor are skipped (their product is exactly 0). The result is a rigorous enclosure
# whose width is of a few ulps for (almost) thin factors, comparable to interval arithmetic;
# for wide factors it can be wider (midpoint-radius product). Because of this, `*` and `mul!`
# are NOT modified: the kernel is only used through the explicit functions `TS.mul_midrad`
# and `TS.mul_midrad!` (for `TaylorN` and `HomogeneousPolynomial`). It is used only if all
# the coefficients of both factors have finite endpoints and a safe decoration (`def`, `dac`
# or `com`); otherwise (or if an overflow occurs) the generic interval arithmetic is used.
# ----------------------------------------------------------------------------
all(f -> isdefined(TS, f), (:mul_midrad, :mul_midrad!, :midrad_poly, :midrad_usable,
        :midrad_acc, :midrad_reset!, :mul_midrad_acc!, :midrad_finalize!)) || error("""
    TaylorSeriesIAExt needs, in TaylorSeries, the stubs
        function mul_midrad end
        function mul_midrad! end
        function midrad_poly end
        function midrad_usable end
        function midrad_acc end
        function midrad_reset! end
        function mul_midrad_acc! end
        function midrad_finalize! end""")

# midpoint and radius of a bounded interval: x ⊆ [m-r, m+r]
@inline function _midrad(x::Interval{T}) where {T<:Base.IEEEFloat}
    lo = inf(x)
    hi = sup(x)
    m = lo/2 + hi/2
    r = max(hi - m, m - lo)
    # upper bound of the exact radius: r is within a relative error u of it, and
    # r*(1+2eps) rounds to something at least r(1+3u) (r = 0 stays 0: thin factors)
    return m, r * (one(T) + 2*eps(T))
end

# (usable, decoration, guaranteed) of the factors for the midpoint-radius kernel: the endpoints
# must be finite (so nonempty and bounded) and the decoration safe (`def`, `dac` or `com`;
# not `trv`/`ill`). The decoration of the result is the minimum of those of the factors,
# as for the interval product.
@inline function _mr_info(x::Interval)
    d = decoration(x)
    return (isfinite(inf(x)) & isfinite(sup(x)) & (d >= def)), d, isguaranteed(x)
end

@inline function _mr_info(coeffs::AbstractVector{<:Interval})
    ok = true
    d = com
    g = true
    @inbounds for x in coeffs
        dx = decoration(x)
        ok &= isfinite(inf(x)) & isfinite(sup(x)) & (dx >= def)
        d = min(d, dx)
        g &= isguaranteed(x)
    end
    return ok, d, g
end

@inline _mr_join(a, b) = (a[1] & b[1], min(a[2], b[2]), a[3] & b[3])

# M - rad rounded downwards / M + rad rounded upwards, only if inexact (TwoSum)
@inline function _mr_down(M::T, rad::T) where {T<:Base.IEEEFloat}
    y = -rad
    s = M + y
    bb = s - M
    t = (M - (s - bb)) + (y - bb)
    return t < 0 ? prevfloat(s) : s
end
@inline function _mr_up(M::T, rad::T) where {T<:Base.IEEEFloat}
    s = M + rad
    bb = s - M
    t = (M - (s - bb)) + (rad - bb)
    return t > 0 ? nextfloat(s) : s
end

# Interval enclosing Σ x_k y_k from M, E, R of `n` terms (see above); `uf`: underflow
# is possible; `d`, `g`: decoration and guarantee flag of the result; `nothing` if it is
# not finite
@inline function _mr_enclose(M::T, E::T, R::T, n::Int, uf::Bool, d::Decoration,
        g::Bool) where {T<:Base.IEEEFloat}
    u = eps(T)/2
    rad = (E + R) * (1 + 4*(n + 4)*u)
    uf && (rad += 16 * n * nextfloat(zero(T)))
    if iszero(rad)
        isfinite(M) || return nothing
        lo = hi = M
    else
        lo = _mr_down(M, rad)
        hi = _mr_up(M, rad)
        (isfinite(lo) && isfinite(hi)) || return nothing
    end
    return IntervalArithmetic._unsafe_interval(
        IntervalArithmetic._unsafe_bareinterval(T, lo, hi), d, g)
end

# Generic kernel (same as in TaylorSeries)
@inline function _mul_output_major_generic!(c::HomogeneousPolynomial, a::HomogeneousPolynomial,
        b::HomogeneousPolynomial, table)
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

# Midpoint-radius kernel; `d`, `g` are the decoration and the `isguaranteed` flag of the result
@inline function _mul_output_major_midrad!(c::HomogeneousPolynomial{Interval{T}},
        a::HomogeneousPolynomial{Interval{T}}, b::HomogeneousPolynomial{Interval{T}},
        table, d::Decoration, g::Bool) where {T<:Base.IEEEFloat}
    offsets = table.output_offsets
    output_pairs = table.output_pairs
    num_right = table.num_right
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    b_coeffs = b.coeffs
    thr = floatmin(T) / eps(T)                  # below this, underflow is possible
    @inbounds for pos in 1:length(offsets)-1
        M = zero(T)
        E = zero(T)
        R = zero(T)
        n = 0
        uf = false
        for csr_pos in offsets[pos]:(offsets[pos+1]-1)
            pair = Int(output_pairs[csr_pos]) - 1
            na = pair ÷ num_right + 1
            nb = pair - (na-1) * num_right + 1
            ma, ra = _midrad(a_coeffs[na])
            mb, rb = _midrad(b_coeffs[nb])
            ((iszero(ma) & iszero(ra)) | (iszero(mb) & iszero(rb))) && continue   # exact 0
            p = ma * mb
            e = fma(ma, mb, -p)                    # ma*mb = p + e exactly
            s = M + p
            bb = s - M
            tt = (M - (s - bb)) + (p - bb)         # M + p = s + tt exactly
            M = s
            E += abs(e) + abs(tt)
            t123 = ra * abs(mb) + abs(ma) * rb + ra * rb
            R += t123
            uf |= (abs(p) < thr) | ((t123 < thr) & (!iszero(ra) | !iszero(rb)))
            n += 1
        end
        n == 0 && continue
        prod = _mr_enclose(M, E, R, n, uf, d, g)
        if prod === nothing                       # overflow: generic interval arithmetic
            acc = c_coeffs[pos]
            for csr_pos in offsets[pos]:(offsets[pos+1]-1)
                pair = Int(output_pairs[csr_pos]) - 1
                na = pair ÷ num_right + 1
                nb = pair - (na-1) * num_right + 1
                acc += a_coeffs[na] * b_coeffs[nb]
            end
            c_coeffs[pos] = acc
        else
            c_coeffs[pos] += prod
        end
    end
    return nothing
end

@inline function _mul_output_major_dispatch!(c::HomogeneousPolynomial{Interval{T}},
        a::HomogeneousPolynomial{Interval{T}}, b::HomogeneousPolynomial{Interval{T}}) where
        {T<:Base.IEEEFloat}
    (TS._isthinzero(b) || TS._isthinzero(a)) && return nothing
    degree_a = order(a)
    degree_b = order(b)
    degree_a == 0 && return _muladd_scalar_dispatch!(c, a.coeffs[1], b)
    degree_b == 0 && return _muladd_scalar_dispatch!(c, b.coeffs[1], a)
    table = TS._init_output_major_product_table!(c.space, degree_a, degree_b)
    ok, d, g = _mr_join(_mr_info(a.coeffs), _mr_info(b.coeffs))
    if ok
        _mul_output_major_midrad!(c, a, b, table, d, g)
    else
        _mul_output_major_generic!(c, a, b, table)
    end
    return nothing
end

# c += scalar * a  (one of the factors is of degree 0)
@inline function _muladd_scalar_dispatch!(c::HomogeneousPolynomial{Interval{T}},
        scalar::Interval{T}, a::HomogeneousPolynomial{Interval{T}}) where
        {T<:Base.IEEEFloat}
    TS._isthinzero(scalar) && return nothing
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    ok, d, g = _mr_join(_mr_info(scalar), _mr_info(a_coeffs))
    if !ok
        @inbounds for i in eachindex(c_coeffs)
            ai = a_coeffs[i]
            TS._isthinzero(ai) && continue
            c_coeffs[i] += scalar * ai
        end
        return nothing
    end
    ms, rs = _midrad(scalar)
    thr = floatmin(T) / eps(T)
    @inbounds for i in eachindex(c_coeffs)
        ai = a_coeffs[i]
        TS._isthinzero(ai) && continue
        ma, ra = _midrad(ai)
        p = ms * ma
        e = fma(ms, ma, -p)                        # ms*ma = p + e exactly
        t123 = rs * abs(ma) + abs(ms) * ra + rs * ra
        uf = (abs(p) < thr) | ((t123 < thr) & (!iszero(rs) | !iszero(ra)))
        prod = _mr_enclose(p, abs(e), t123, 1, uf, d, g)
        c_coeffs[i] += prod === nothing ? scalar * ai : prod
    end
    return nothing
end

# c += a * b with the midpoint-radius kernel (homogeneous polynomials)
TS.mul_midrad!(c::HomogeneousPolynomial{Interval{T}}, a::HomogeneousPolynomial{Interval{T}},
        b::HomogeneousPolynomial{Interval{T}}) where {T<:Base.IEEEFloat} =
    _mul_output_major_dispatch!(c, a, b)

# c += a * b with the midpoint-radius kernel (TaylorN)
function TS.mul_midrad!(c::TaylorN{Interval{T}}, a::TaylorN{Interval{T}},
        b::TaylorN{Interval{T}}) where {T<:Base.IEEEFloat}
    TS._check_same_space(c, a, b)
    for k in eachindex(c)
        kk = k + 1
        @inbounds for i = 0:k
            _mul_output_major_dispatch!(c.coeffs[kk], a.coeffs[i+1], b.coeffs[kk-i])
        end
    end
    return nothing
end

# a * b with the midpoint-radius kernel
function TS.mul_midrad(a::TaylorN{Interval{T}}, b::TaylorN{Interval{T}}) where
        {T<:Base.IEEEFloat}
    TS._check_same_space(a, b)
    if TS.order(a) != TS.order(b)
        a, b = TS.fixorder(a, b)
    end
    c = zero(a)
    TS.mul_midrad!(c, a, b)
    return c
end


# ----------------------------------------------------------------------------
# Accumulation of products of whole polynomials in midpoint-radius form
#
# `TS.midrad_poly(p)` converts a `TaylorN{Interval}` once to midpoints and radii (per degree);
# `TS.midrad_acc(p, maxdeg)` creates accumulators `M, E, R, n` (see above) for each monomial of
# the degrees `0:maxdeg` of the products; `TS.mul_midrad_acc!(acc, pa, pb[, factor])` accumulates
# `factor * pa * pb` (`factor` a power of two: exact) for all the pairs of degrees with
# `u+v ≤ maxdeg`, with no interval arithmetic at all; `TS.midrad_finalize!(dest, acc, d)`
# writes into the homogeneous polynomial `dest` of degree `d` one interval per monomial (one
# `_mr_enclose` for the whole accumulation, so a single rigorous bound of the sum of all the
# products accumulated); `TS.midrad_reset!(acc)` clears the accumulators. The error-free
# transformations make the order of the accumulation irrelevant for the validity of the bound.
# The decoration of the result is the minimum of those of all the factors accumulated and its
# `isguaranteed` flag their conjunction. `midrad_finalize!` returns `false` (and the caller must
# use the interval arithmetic) if some result is not finite.
# ----------------------------------------------------------------------------
struct MidRadPoly{T<:Base.IEEEFloat,S}
    mid::Vector{Vector{T}}
    rad::Vector{Vector{T}}
    ok::Bool
    dec::Decoration
    g::Bool
    space::S
end

mutable struct MidRadAcc{T<:Base.IEEEFloat}
    M::Vector{Vector{T}}
    E::Vector{Vector{T}}
    R::Vector{Vector{T}}
    n::Vector{Vector{Int}}
    uf::Vector{Vector{Bool}}
    maxdeg::Int
    dec::Decoration
    g::Bool
end

function TS.midrad_poly(p::TaylorN{Interval{T}}) where {T<:Base.IEEEFloat}
    nd = length(p.coeffs)
    mid = Vector{Vector{T}}(undef, nd)
    rad = Vector{Vector{T}}(undef, nd)
    ok = true
    d = com
    g = true
    for k in 1:nd
        cs = p.coeffs[k].coeffs
        m = Vector{T}(undef, length(cs))
        r = Vector{T}(undef, length(cs))
        @inbounds for i in eachindex(cs)
            x = cs[i]
            dx = decoration(x)
            ok &= isfinite(inf(x)) & isfinite(sup(x)) & (dx >= def)
            d = min(d, dx)
            g &= isguaranteed(x)
            m[i], r[i] = _midrad(x)
        end
        mid[k] = m
        rad[k] = r
    end
    sp = p.coeffs[1].space
    return MidRadPoly{T,typeof(sp)}(mid, rad, ok, d, g, sp)
end

TS.midrad_usable(mp::MidRadPoly) = mp.ok

function TS.midrad_acc(::TaylorN{Interval{T}}, maxdeg::Int) where {T<:Base.IEEEFloat}
    mk() = [T[] for _ in 0:maxdeg]
    return MidRadAcc{T}(mk(), mk(), mk(), [Int[] for _ in 0:maxdeg],
        [Bool[] for _ in 0:maxdeg], maxdeg, com, true)
end

function TS.midrad_reset!(acc::MidRadAcc{T}) where {T}
    for d in 1:acc.maxdeg+1
        fill!(acc.M[d], zero(T))
        fill!(acc.E[d], zero(T))
        fill!(acc.R[d], zero(T))
        fill!(acc.n[d], 0)
        fill!(acc.uf[d], false)
    end
    acc.dec = com
    acc.g = true
    return nothing
end

@inline function _acc_arrays!(acc::MidRadAcc{T}, d::Int, nout::Int) where {T}
    Md = acc.M[d+1]
    if length(Md) != nout
        resize!(Md, nout); fill!(Md, zero(T))
        resize!(acc.E[d+1], nout); fill!(acc.E[d+1], zero(T))
        resize!(acc.R[d+1], nout); fill!(acc.R[d+1], zero(T))
        resize!(acc.n[d+1], nout); fill!(acc.n[d+1], 0)
        resize!(acc.uf[d+1], nout); fill!(acc.uf[d+1], false)
    end
    return Md, acc.E[d+1], acc.R[d+1], acc.n[d+1], acc.uf[d+1]
end

# accumulate the term (ma ± ra)(mb ± rb) in the position `pos` (see the section above)
@inline function _acc_term!(Md, Ed, Rd, nd, ufd, pos::Int, ma::T, ra::T, mb::T, rb::T,
        thr::T) where {T}
    ((iszero(ma) & iszero(ra)) | (iszero(mb) & iszero(rb))) && return nothing
    p = ma * mb
    e = fma(ma, mb, -p)
    @inbounds begin
        M = Md[pos]
        s = M + p
        bb = s - M
        tt = (M - (s - bb)) + (p - bb)
        Md[pos] = s
        Ed[pos] += abs(e) + abs(tt)
        t123 = ra * abs(mb) + abs(ma) * rb + ra * rb
        Rd[pos] += t123
        ufd[pos] |= (abs(p) < thr) | ((t123 < thr) & (!iszero(ra) | !iszero(rb)))
        nd[pos] += 1
    end
    return nothing
end

function TS.mul_midrad_acc!(acc::MidRadAcc{T}, pa::MidRadPoly{T}, pb::MidRadPoly{T},
        factor::T = one(T)) where {T<:Base.IEEEFloat}
    acc.dec = min(acc.dec, pa.dec, pb.dec)
    acc.g &= pa.g & pb.g
    thr = floatmin(T) / eps(T)
    for u in 0:length(pa.mid)-1, v in 0:length(pb.mid)-1
        d = u + v
        d > acc.maxdeg && continue
        mau, rau = pa.mid[u+1], pa.rad[u+1]
        mbv, rbv = pb.mid[v+1], pb.rad[v+1]
        if u == 0
            Md, Ed, Rd, nd, ufd = _acc_arrays!(acc, d, length(mbv))
            msc, rsc = factor * mau[1], factor * rau[1]
            @inbounds for i in eachindex(mbv)
                _acc_term!(Md, Ed, Rd, nd, ufd, i, msc, rsc, mbv[i], rbv[i], thr)
            end
        elseif v == 0
            Md, Ed, Rd, nd, ufd = _acc_arrays!(acc, d, length(mau))
            msc, rsc = mbv[1], rbv[1]
            @inbounds for i in eachindex(mau)
                _acc_term!(Md, Ed, Rd, nd, ufd, i, factor * mau[i], factor * rau[i],
                    msc, rsc, thr)
            end
        else
            table = TS._init_output_major_product_table!(pa.space, u, v)
            offsets = table.output_offsets
            output_pairs = table.output_pairs
            num_right = table.num_right
            Md, Ed, Rd, nd, ufd = _acc_arrays!(acc, d, length(offsets) - 1)
            @inbounds for pos in 1:length(offsets)-1
                for csr_pos in offsets[pos]:(offsets[pos+1]-1)
                    pair = Int(output_pairs[csr_pos]) - 1
                    na = pair ÷ num_right + 1
                    nb = pair - (na-1) * num_right + 1
                    _acc_term!(Md, Ed, Rd, nd, ufd, pos, factor * mau[na], factor * rau[na],
                        mbv[nb], rbv[nb], thr)
                end
            end
        end
    end
    return nothing
end

function TS.midrad_finalize!(dest::HomogeneousPolynomial{Interval{T}}, acc::MidRadAcc{T},
        d::Int) where {T<:Base.IEEEFloat}
    (d > acc.maxdeg || isempty(acc.M[d+1])) && return true       # nothing accumulated
    Md, Ed, Rd, nd, ufd = acc.M[d+1], acc.E[d+1], acc.R[d+1], acc.n[d+1], acc.uf[d+1]
    @assert length(Md) == length(dest.coeffs)
    @inbounds for pos in eachindex(Md)
        n = nd[pos]
        n == 0 && continue
        prod = _mr_enclose(Md[pos], Ed[pos], Rd[pos], n, ufd[pos], acc.dec, acc.g)
        prod === nothing && return false
        dest.coeffs[pos] = prod
    end
    return true
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
    order = TS.order(a)
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
        TS.sqrt!(c, aa, zero(a0), k, lnull)
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
    order = TS.order(a)
    c = TaylorN(space(a), zero(aux), order)
    for k in eachindex(aa)
        TS.sqrt!(c, aa, zero(a0), k)
    end
    return c
end

function TS.sqrt!(c::Taylor1{Interval{T}}, a::Taylor1{Interval{T}},
        ::Interval{T}, k::Int, k0::Int=0) where {T<:NumTypes}
    if k == k0
        @inbounds c[k] = sqrt(a[2*k0])
        return nothing
    end
    kodd = (k - k0)%2
    kend = div(k - k0 - 2 + kodd, 2)
    a_order = order(a)
    imax = min(k0+kend, a_order)
    imin = max(k0+1, k+k0-a_order)
    imin ≤ imax && ( @inbounds c[k] = c[imin] * c[k+k0-imin] )
    @inbounds for i = imin+1:imax
        c[k] += c[i] * c[k+k0-i]
    end
    intvl2 = interval(T(2))
    if k+k0 ≤ a_order
        @inbounds aux = a[k+k0] - intvl2 * c[k]
    else
        @inbounds aux = - intvl2 * c[k]
    end
    if kodd == 0
        @inbounds aux = aux - c[kend+k0+1]^2
    end
    @inbounds c[k] = aux / (intvl2 * c[k0])
    return nothing
end

function TS.sqrt!(c::TaylorN{Interval{T}}, a::TaylorN{Interval{T}},
        ::Interval{T}, k::Int) where {T<:NumTypes}
    TS._check_same_space(c, a)
    if k == 0
        @inbounds c[0][1] = sqrt( constant_term(a) )
        return nothing
    end
    # Recursion formula
    kodd = k%2
    kend = (k - 2 + kodd) >> 1
    # c[k] <- a[k]
    @inbounds for i in eachindex(c[k])
        c[k][i] = a[k][i]
    end
    if kodd == 0
        # @inbounds c[k] <- c[k] - (c[kend+1])^2
        @inbounds TS.mul_scalar!(c[k], -interval(T(1)), c[kend+1], c[kend+1])
    end
    intvl2 = interval(T(2))
    @inbounds for i = 1:kend
        # c[k] <- c[k] - 2*c[i]*c[k-i]
        TS.mul_scalar!(c[k], -intvl2, c[i], c[k-i])
    end
    # @inbounds c[k] <- c[k] / (2*c[0])
    TS.div!(c[k], c[k], intvl2*constant_term(c))

    return nothing
end

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
            # aa = one(aux) * a
            aa[0] = one(aux) * a0
            order = TS.order(a)
            c = TS._constant_series_like(a, aux, order)
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
            order = TS.order(a)
            aa = convert($T{typeof(aux)}, a)
            aa[0] = uno * a0
            c = TS._constant_series_like(a, aux, order)
            r = TS._constant_series_like(a, sqrt(uno - a0sqr), order)
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
            order = TS.order(a)
            aa = convert($T{typeof(aux)}, a)
            aa[0] = uno * a0
            c = TS._constant_series_like(a, aux, order)
            r = TS._constant_series_like(a, sqrt(uno - a0sqr), order)
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
            order = TS.order(a)
            aa = convert($T{typeof(aux)}, a)
            aa[0] = uno * a0
            c = TS._constant_series_like(a, aux, order)
            r = TS._constant_series_like(a, sqrt(a0sqr - uno), order)
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
            order = TS.order(a)
            uno = one(a0)
            aa = convert($T{typeof(aux)}, a)
            aa[0] = uno * a0
            c = TS._constant_series_like(a, aux, order)
            r = TS._constant_series_like(a, uno - a0^2, order)
            TS._isthinzero(constant_term(r)) && throw(DomainError(a,
                """Series expansion of atanh(x) diverges at x = ±1."""))
            for k in eachindex(a)
                TS.atanh!(c, aa, r, k)
            end
            return c
        end

    end
end

# Some internal functions
function TS.exp!(c::Taylor1{Interval{T}}, a::Taylor1{Interval{T}}, k::Int) where {T<:NumTypes}
    if k == 0
        @inbounds c[0] = exp(constant_term(a))
        return nothing
    end
    intvlk = interval(T(k))
    @inbounds c[k] = intvlk * a[k] * c[0]
    @inbounds for i = 1:k-1
        c[k] += interval(T(k-i)) * a[k-i] * c[i]
    end
    @inbounds c[k] = c[k] / intvlk
    return nothing
end

function TS.exp!(c::TaylorN{Interval{T}}, a::TaylorN{Interval{T}}, k::Int) where {T<:NumTypes}
    if k == 0
        @inbounds c[0] = exp(constant_term(a))
        return nothing
    end
    intvlk = interval(T(k))
    @inbounds TS.mul!(c[k], intvlk * a[k], c[0])
    @inbounds for i = 1:k-1
        TS.mul!(c[k], interval(T(k-i)) * a[k-i], c[i])
    end
    @inbounds c[k] = c[k] / intvlk
    return nothing
end

function TS.expm1!(c::Taylor1{Interval{T}}, a::Taylor1{Interval{T}}, k::Int) where {T<:NumTypes}
    if k == 0
        @inbounds c[0] = expm1(constant_term(a))
        return nothing
    end
    c0 = c[0]+one(c[0])
    intvlk = interval(T(k))
    @inbounds c[k] = intvlk * a[k] * c0
    @inbounds for i = 1:k-1
        c[k] += interval(T(k-i)) * a[k-i] * c[i]
    end
    @inbounds c[k] = c[k] / intvlk
    return nothing
end

function TS.expm1!(c::TaylorN{Interval{T}}, a::TaylorN{Interval{T}}, k::Int) where {T<:NumTypes}
    if k == 0
        @inbounds c[0] = expm1(constant_term(a))
        return nothing
    end
    c0 = c[0]+one(c[0])
    intvlk = interval(T(k))
    @inbounds TS.mul!(c[k], intvlk * a[k], c0)
    @inbounds for i = 1:k-1
        TS.mul!(c[k], interval(T(k-i)) * a[k-i], c[i])
    end
    @inbounds c[k] = c[k] / intvlk
    return nothing
end

function TS.log!(c::Taylor1{Interval{T}}, a::Taylor1{Interval{T}}, k::Int) where {T<:NumTypes}
    if k == 0
        @inbounds c[0] = log(constant_term(a))
        return nothing
    elseif k == 1
        @inbounds c[1] = a[1] / constant_term(a)
        return nothing
    end
    @inbounds c[k] = interval(T(k-1)) * a[1] * c[k-1]
    @inbounds for i = 2:k-1
        c[k] += interval(T(k-i)) * a[i] * c[k-i]
    end
    @inbounds c[k] = (a[k] - c[k]/interval(T(k))) / constant_term(a)
    return nothing
end

function TS.log!(c::TaylorN{Interval{T}}, a::TaylorN{Interval{T}}, k::Int) where {T<:NumTypes}
    if k == 0
        @inbounds c[0] = log(constant_term(a))
        return nothing
    elseif k == 1
        @inbounds c[1] = a[1] / constant_term(a)
        return nothing
    end
    @inbounds TS.mul!(c[k], interval(T(k-1))*a[1], c[k-1])
    @inbounds for i = 2:k-1
        TS.mul!(c[k], interval(T(k-i))*a[i], c[k-i])
    end
    @inbounds c[k] = (a[k] - c[k]/interval(T(k))) / constant_term(a)
    return nothing
end

function TS.log1p!(c::Taylor1{Interval{T}}, a::Taylor1{Interval{T}}, k::Int) where {T<:NumTypes}
    a0 = constant_term(a)
    a0p1 = a0+one(a0)
    if k == 0
        @inbounds c[0] = log1p(a0)
        return nothing
    elseif k == 1
        @inbounds c[1] = a[1] / a0p1
        return nothing
    end
    @inbounds c[k] = interval(T(k-1)) * a[1] * c[k-1]
    @inbounds for i = 2:k-1
        c[k] += interval(T(k-i)) * a[i] * c[k-i]
    end
    @inbounds c[k] = (a[k] - c[k]/interval(T(k))) / a0p1
    return nothing
end

function TS.log1p!(c::TaylorN{Interval{T}}, a::TaylorN{Interval{T}}, k::Int) where {T<:NumTypes}
    a0 = constant_term(a)
    a0p1 = a0+one(a0)
    if k == 0
        @inbounds c[0] = log1p(a0)
        return nothing
    elseif k == 1
        @inbounds c[1] = a[1] / a0p1
        return nothing
    end
    @inbounds TS.mul!(c[k], interval(T(k-1))*a[1], c[k-1])
    @inbounds for i = 2:k-1
        TS.mul!(c[k], interval(T(k-i))*a[i], c[k-i])
    end
    @inbounds c[k] = (a[k] - c[k]/interval(T(k))) / a0p1
    return nothing
end

function TS.sincos!(s::Taylor1{Interval{T}}, c::Taylor1{Interval{T}},
        a::Taylor1{Interval{T}}, k::Int) where {T<:NumTypes}
    if k == 0
        a0 = constant_term(a)
        @inbounds s[0], c[0] = sincos( a0 )
        return nothing
    end
    x = a[1]
    @inbounds s[k] = x * c[k-1]
    @inbounds c[k] = -x * s[k-1]
    @inbounds for i = 2:k
        x = interval(T(i)) * a[i]
        s[k] += x * c[k-i]
        c[k] -= x * s[k-i]
    end
    intvlk = interval(T(k))
    @inbounds s[k] = s[k] / intvlk
    @inbounds c[k] = c[k] / intvlk
    return nothing
end

function TS.sincos!(s::TaylorN{Interval{T}}, c::TaylorN{Interval{T}},
        a::TaylorN{Interval{T}}, k::Int) where {T<:NumTypes}
    if k == 0
        a0 = constant_term(a)
        @inbounds s[0], c[0] = sincos( a0 )
        return nothing
    end
    x = a[1]
    TS.mul!(s[k], x, c[k-1])
    TS.mul!(c[k], -x, s[k-1])
    @inbounds for i = 2:k
        x = interval(T(i)) * a[i]
        TS.mul!(s[k], x, c[k-i])
        TS.mul!(c[k], -x, s[k-i])
    end
    intvlk = interval(T(k))
    @inbounds s[k] = s[k] / intvlk
    @inbounds c[k] = c[k] / intvlk
    return nothing
end

function TS.tan!(c::Taylor1{Interval{T}}, a::Taylor1{Interval{T}},
        c2::Taylor1{Interval{T}}, k::Int) where {T<:NumTypes}
    if k == 0
        @inbounds aux = tan( constant_term(a) )
        @inbounds c[0] = aux
        @inbounds c2[0] = aux^2
        return nothing
    end
    intvlk = interval(T(k))
    @inbounds c[k] = intvlk * a[k] * c2[0]
    @inbounds for i = 1:k-1
        c[k] += interval(T(k-i)) * a[k-i] * c2[i]
    end
    @inbounds c[k] = a[k] + c[k]/intvlk
    TS.sqr!(c2, c, zero(c[0]), k)
    return nothing
end

function TS.tan!(c::TaylorN{Interval{T}}, a::TaylorN{Interval{T}},
        c2::TaylorN{Interval{T}}, k::Int) where {T<:NumTypes}
    if k == 0
        @inbounds aux = tan( constant_term(a) )
        @inbounds c[0] = aux
        @inbounds c2[0] = aux^2
        return nothing
    end
    intvlk = interval(T(k))
    @inbounds TS.mul!(c[k], intvlk * a[k], c2[0])
    @inbounds for i = 1:k-1
        TS.mul!(c[k], interval(T(k-i)) * a[k-i], c2[i])
    end
    @inbounds c[k] = a[k] + c[k]/intvlk
    TS.sqr!(c2, c, zero(c[0][1]), k)
    return nothing
end

function TS.asin!(c::Taylor1{Interval{T}}, a::Taylor1{Interval{T}},
        r::Taylor1{Interval{T}}, k::Int) where {T<:NumTypes}
    a0 = constant_term(a)
    if k == 0
        @inbounds c[0] = asin( a0 )
        @inbounds r[0] = sqrt( one(a0) - a0^2 )
        return nothing
    end
    @inbounds c[k] = interval(T(k-1)) * r[1] * c[k-1]
    @inbounds for i in 2:k-1
        c[k] += interval(T(k-i)) * r[i] * c[k-i]
    end
    TS.sqrt!(r, one(a[0])-a^2, zero(a0), k)
    @inbounds c[k] = (a[k] - c[k]/interval(T(k))) / constant_term(r)
    return nothing
end

function TS.asin!(c::TaylorN{Interval{T}}, a::TaylorN{Interval{T}},
        r::TaylorN{Interval{T}}, k::Int) where {T<:NumTypes}
    a0 = constant_term(a)
    if k == 0
        @inbounds c[0] = asin( a0 )
        @inbounds r[0] = sqrt( one(a0) - a0^2 )
        return nothing
    end
    @inbounds TS.mul!(c[k], interval(T(k-1)) * r[1], c[k-1])
    @inbounds for i in 2:k-1
        TS.mul!(c[k], interval(T(k-i)) * r[i], c[k-i])
    end
    TS.sqrt!(r, one(a[0])-a^2, zero(a0), k)
    @inbounds c[k] = (a[k] - c[k]/interval(T(k))) / constant_term(r)
    return nothing
end

function TS.acos!(c::Taylor1{Interval{T}}, a::Taylor1{Interval{T}},
        r::Taylor1{Interval{T}}, k::Int) where {T<:NumTypes}
    a0 = constant_term(a)
    if k == 0
        @inbounds c[0] = acos( a0 )
        @inbounds r[0] = sqrt( one(a0) - a0^2 )
        return nothing
    end
    @inbounds c[k] = interval(T(k-1)) * r[1] * c[k-1]
    @inbounds for i in 2:k-1
        c[k] += interval(T(k-i)) * r[i] * c[k-i]
    end
    TS.sqrt!(r, one(a[0])-a^2, zero(a0), k)
    @inbounds c[k] = -(a[k] + c[k]/interval(T(k))) / constant_term(r)
    return nothing
end

function TS.acos!(c::TaylorN{Interval{T}}, a::TaylorN{Interval{T}},
        r::TaylorN{Interval{T}}, k::Int) where {T<:NumTypes}
    a0 = constant_term(a)
    if k == 0
        @inbounds c[0] = acos( a0 )
        @inbounds r[0] = sqrt( one(a0) - a0^2 )
        return nothing
    end
    @inbounds TS.mul!(c[k], interval(T(k-1)) * r[1], c[k-1])
    @inbounds for i in 2:k-1
        TS.mul!(c[k], interval(T(k-i)) * r[i], c[k-i])
    end
    TS.sqrt!(r, one(a[0])-a^2, zero(a0), k)
    @inbounds c[k] = -(a[k] + c[k]/interval(T(k))) / constant_term(r)
    return nothing
end

function TS.atan!(c::Taylor1{Interval{T}}, a::Taylor1{Interval{T}},
        r::Taylor1{Interval{T}}, k::Int) where {T<:NumTypes}
    if k == 0
        a0 = constant_term(a)
        @inbounds c[0] = atan( a0 )
        @inbounds r[0] = one(a0) + a0^2
        return nothing
    end
    @inbounds c[k] = interval(T(k-1)) * r[1] * c[k-1]
    @inbounds for i in 2:k-1
        c[k] += interval(T(k-i)) * r[i] * c[k-i]
    end
    TS.sqr!(r, a, zero(a[0]), k)
    @inbounds c[k] = (a[k] - c[k]/interval(T(k))) / constant_term(r)
    return nothing
end

function TS.atan!(c::TaylorN{Interval{T}}, a::TaylorN{Interval{T}},
        r::TaylorN{Interval{T}}, k::Int) where {T<:NumTypes}
    if k == 0
        a0 = constant_term(a)
        @inbounds c[0] = atan( a0 )
        @inbounds r[0] = one(a0) + a0^2
        return nothing
    end
    @inbounds TS.mul!(c[k], interval(T(k-1)) * r[1], c[k-1])
    @inbounds for i in 2:k-1
        TS.mul!(c[k], interval(T(k-i)) * r[i], c[k-i])
    end
    TS.sqr!(r, a, zero(a[0][1]), k)
    @inbounds c[k] = (a[k] - c[k]/interval(T(k))) / constant_term(r)
    return nothing
end

function TS.sinhcosh!(s::Taylor1{Interval{T}}, c::Taylor1{Interval{T}},
        a::Taylor1{Interval{T}}, k::Int) where {T<:NumTypes}
    if k == 0
        @inbounds s[0] = sinh( constant_term(a) )
        @inbounds c[0] = cosh( constant_term(a) )
        return nothing
    end
    x = a[1]
    @inbounds s[k] = x * c[k-1]
    @inbounds c[k] = x * s[k-1]
    @inbounds for i = 2:k
        x = interval(T(i)) * a[i]
        s[k] += x * c[k-i]
        c[k] += x * s[k-i]
    end
    intvlk = interval(T(k))
    s[k] = s[k] / intvlk
    c[k] = c[k] / intvlk
    return nothing
end

function TS.sinhcosh!(s::TaylorN{Interval{T}}, c::TaylorN{Interval{T}},
        a::TaylorN{Interval{T}}, k::Int) where {T<:NumTypes}
    if k == 0
        @inbounds s[0] = sinh( constant_term(a) )
        @inbounds c[0] = cosh( constant_term(a) )
        return nothing
    end
    x = a[1]
    @inbounds TS.mul!(s[k], x, c[k-1])
    @inbounds TS.mul!(c[k], x, s[k-1])
    @inbounds for i = 2:k
        x = interval(T(i)) * a[i]
        TS.mul!(s[k], x, c[k-i])
        TS.mul!(c[k], x, s[k-i])
    end
    intvlk = interval(T(k))
    s[k] = s[k] / intvlk
    c[k] = c[k] / intvlk
    return nothing
end

function TS.tanh!(c::Taylor1{Interval{T}}, a::Taylor1{Interval{T}},
        c2::Taylor1{Interval{T}}, k::Int) where {T<:NumTypes}
    if k == 0
        @inbounds aux = tanh( constant_term(a) )
        @inbounds c[0] = aux
        @inbounds c2[0] = aux^2
        return nothing
    end
    @inbounds c[k] = k * a[k] * c2[0]
    @inbounds for i = 1:k-1
        c[k] += (k-i) * a[k-i] * c2[i]
    end
    @inbounds c[k] = a[k] - c[k]/k
    TS.sqr!(c2, c, zero(c[0]), k)
    return nothing
end

function TS.tanh!(c::TaylorN{Interval{T}}, a::TaylorN{Interval{T}},
        c2::TaylorN{Interval{T}}, k::Int) where {T<:NumTypes}
    if k == 0
        @inbounds aux = tanh( constant_term(a) )
        @inbounds c[0] = aux
        @inbounds c2[0] = aux^2
        return nothing
    end
    @inbounds TS.mul!(c[k], k * a[k], c2[0])
    @inbounds for i = 1:k-1
        TS.mul!(c[k], (k-i) * a[k-i], c2[i])
    end
    @inbounds c[k] = a[k] - c[k]/k
    TS.sqr!(c2, c, zero(c[0][1]), k)
    return nothing
end

function TS.asinh!(c::Taylor1{Interval{T}}, a::Taylor1{Interval{T}},
        r::Taylor1{Interval{T}}, k::Int) where {T<:NumTypes}
    a0 = constant_term(a)
    if k == 0
        @inbounds c[0] = asinh( a0 )
        @inbounds r[0] = sqrt( a0^2 + one(a0) )
        return nothing
    end
    @inbounds c[k] = interval(T(k-1)) * r[1] * c[k-1]
    @inbounds for i in 2:k-1
        c[k] += interval(T(k-i)) * r[i] * c[k-i]
    end
    TS.sqrt!(r, a^2+one(a[0]), zero(a0), k)
    @inbounds c[k] = (a[k] - c[k]/interval(T(k))) / constant_term(r)
    return nothing
end

function TS.asinh!(c::TaylorN{Interval{T}}, a::TaylorN{Interval{T}},
        r::TaylorN{Interval{T}}, k::Int) where {T<:NumTypes}
    a0 = constant_term(a)
    if k == 0
        @inbounds c[0] = asinh( a0 )
        @inbounds r[0] = sqrt( a0^2 + one(a0) )
        return nothing
    end
    @inbounds TS.mul!(c[k], interval(T(k-1)) * r[1], c[k-1])
    @inbounds for i in 2:k-1
        TS.mul!(c[k], interval(T(k-i)) * r[i], c[k-i])
    end
    TS.sqrt!(r, a^2+one(a[0]), zero(a0), k)
    @inbounds c[k] = (a[k] - c[k]/interval(T(k))) / constant_term(r)
    return nothing
end

function TS.acosh!(c::Taylor1{Interval{T}}, a::Taylor1{Interval{T}},
        r::Taylor1{Interval{T}}, k::Int) where {T<:NumTypes}
    a0 = constant_term(a)
    if k == 0
        @inbounds c[0] = acosh( a0 )
        @inbounds r[0] = sqrt( a0^2 - one(a0) )
        return nothing
    end
    @inbounds c[k] = interval(T(k-1)) * r[1] * c[k-1]
    @inbounds for i in 2:k-1
        c[k] += interval(T(k-i)) * r[i] * c[k-i]
    end
    TS.sqrt!(r, a^2-one(a[0]), zero(a0), k)
    @inbounds c[k] = (a[k] - c[k]/interval(T(k))) / constant_term(r)
    return nothing
end

function TS.acosh!(c::TaylorN{Interval{T}}, a::TaylorN{Interval{T}},
        r::TaylorN{Interval{T}}, k::Int) where {T<:NumTypes}
    a0 = constant_term(a)
    if k == 0
        @inbounds c[0] = acosh( a0 )
        @inbounds r[0] = sqrt( a0^2 - one(a0) )
        return nothing
    end
    @inbounds TS.mul!(c[k], interval(T(k-1)) * r[1], c[k-1])
    @inbounds for i in 2:k-1
        TS.mul!(c[k], interval(T(k-i)) * r[i], c[k-i])
    end
    TS.sqrt!(r, a^2-one(a[0]), zero(a0), k)
    @inbounds c[k] = (a[k] - c[k]/interval(T(k))) / constant_term(r)
    return nothing
end

function TS.atanh!(c::Taylor1{Interval{T}}, a::Taylor1{Interval{T}},
        r::Taylor1{Interval{T}}, k::Int) where {T<:NumTypes}
    a0 = constant_term(a)
    if k == 0
        @inbounds c[0] = atanh( a0 )
        @inbounds r[0] = one(a0) - a0^2
        return nothing
    end
    @inbounds c[k] = interval(T(k-1)) * r[1] * c[k-1]
    @inbounds for i in 2:k-1
        c[k] += interval(T(k-i)) * r[i] * c[k-i]
    end
    TS.sqr!(r, a, zero(a0), k)
    @inbounds c[k] = (a[k] + c[k]/interval(T(k))) / constant_term(r)
    return nothing
end

function TS.atanh!(c::TaylorN{Interval{T}}, a::TaylorN{Interval{T}},
        r::TaylorN{Interval{T}}, k::Int) where {T<:NumTypes}
    a0 = constant_term(a)
    if k == 0
        @inbounds c[0] = atanh( a0 )
        @inbounds r[0] = one(a0) - a0^2
        return nothing
    end
    @inbounds TS.mul!(c[k], interval(T(k-1)) * r[1], c[k-1])
    @inbounds for i in 2:k-1
        TS.mul!(c[k], interval(T(k-i)) * r[i], c[k-i])
    end
    TS.sqr!(r, a, zero(a0), k)
    @inbounds c[k] = (a[k] + c[k]/interval(T(k))) / constant_term(r)
    return nothing
end


function TS.evaluate(a::Taylor1{Interval{T}}, dx::S) where {T<:NumTypes, S<:NumTypes}
    dxI, _ = promote(dx, constant_term(a))
    return evaluate(a, dxI)
end

function TS.evaluate(a::Taylor1{T}, dx::Interval{T}) where {T<:NumTypes}
    order = TS.order(a)
    order == 0 && return a[0] * one(dx)
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

function TS.evaluate(a::Taylor1{Interval{T}}, dx::Interval{T}) where {T<:NumTypes}
    order = TS.order(a)
    order == 0 && return a[0] * one(dx)
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

function TS.evaluate(a::Taylor1{TaylorN{T}}, dx::Interval{S}) where {T<:Real, S<:Real}
    order = TS.order(a)
    order == 0 && return a[0] * one(dx)
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


TS._evaluate(a::HomogeneousPolynomial{T}, dx::AbstractVector{Interval{S}}) where
    {T<:Real, S<:NumTypes} = _evaluate_hp_interval(a, dx)

TS._evaluate(a::HomogeneousPolynomial{T}, dx::Tuple{Interval{S}, Vararg{Interval{S}}}) where
    {T<:Real, S<:NumTypes} = _evaluate_hp_interval(a, dx)

# Use the specialized methods when the box is `[-1,1]^n` or `[0,1]^n`
function _evaluate_hp_interval(a::HomogeneousPolynomial, dx::_IntervalVals{S}) where
        {S<:NumTypes}
    Isym = interval(-one(S), one(S))
    all(x -> isequal_interval(x, Isym), dx) && return TS._evaluate(a, dx, Val(true))
    Ipos = interval(zero(S), one(S))
    all(x -> isequal_interval(x, Ipos), dx) && return TS._evaluate(a, dx, Val(false))
    return _evaluate_hp_generic(a, dx)
end

function _evaluate_hp_generic(a::HomogeneousPolynomial, dx::_IntervalVals{S}) where
        {S<:NumTypes}
    order(a) == 0 && return a[1] + interval(zero(S))
    ct = a.space.coeff_table[order(a)+1]
    suma = zero(a[1]) * interval(zero(S))
    for (i, a_coeff) in enumerate(a.coeffs)
        TS._isthinzero(a_coeff) && continue
        term = a_coeff * one(dx[1])
        @inbounds for (j, x) in enumerate(dx)
            exponent = ct[i][j]
            exponent == 0 && continue
            term *= Base.literal_pow(^, x, Val(exponent))
        end
        suma += term
    end
    return suma
end

# Case `[-1,1]^n`: If the total odd, each monomial ranges over `[-1,1]`;
# if the order is even, its range is `[0,1]` when all exponents are
# even, and `[-1,1]` otherwise.
# Here `Val(true)`/`Val(false)` select the box, not sorting.
function TS._evaluate(a::HomogeneousPolynomial{T}, dx::_IntervalVals{S},
        ::Val{true}) where {T<:Real, S<:NumTypes}
    order(a) == 0 && return a[1] + interval(zero(T))
    ct = a.space.coeff_table[order(a)+1]
    suma = a[1] * interval(zero(S))
    Isym = dx[1]
    Ieven = interval(zero(S), one(S))
    odd_order = isodd(order(a))
    for (i, a_coeff) in enumerate(a.coeffs)
        TS._isthinzero(a_coeff) && continue
        # if isodd(sum(ct[i]))
        #     suma += sum(a_coeff) * dx[1]
        #     continue
        # end
        # @inbounds tmp = iseven(ct[i][1]) ? Ieven : dx[1]
        # for n in 2:length(dx)
        #     @inbounds vv = iseven(ct[i][n]) ? Ieven : dx[1]
        #     tmp *= vv
        # end
        tmp = (odd_order || !all(iseven, ct[i])) ? Isym : Ieven
        suma += a_coeff * tmp
    end
    return suma
end

# `[0,1]^n`: every monomial ranges over `[0,1]`
function TS._evaluate(a::HomogeneousPolynomial{T}, dx::_IntervalVals{S},
        ::Val{false}) where {T<:Real, S<:NumTypes}
    order(a) == 0 && return a[1] + interval(zero(S))
    suma = a[1] * dx[1]
    @inbounds for a_coeff in a.coeffs
        suma += a_coeff * dx[1]
    end
    return suma
end

function TS._evaluate(a::TaylorN{T}, vals::NTuple{N,TaylorN{Interval{S}}}) where
        {N, T<:Real, S<:NumTypes}
    @assert get_numvars(a.space) == N
    TS._check_same_space_all(a, vals)
    R = promote_type(TS.numtype(a), typeof(vals[1]))
    a_length = length(a)
    suma = Vector{R}(undef, a_length)
    @inbounds for homPol in 1:a_length
        suma[homPol] = TS._evaluate(a.coeffs[homPol], vals)
    end
    return suma
end

function TS._evaluate(a::TaylorN{Interval{T}},
        vals::NTuple{N,TaylorN{Interval{T}}}) where {N, T<:NumTypes}
    @assert get_numvars(a.space) == N
    TS._check_same_space_all(a, vals)
    a_length = length(a)
    suma = Vector{TaylorN{Interval{T}}}(undef, a_length)
    @inbounds for homPol in 1:a_length
        suma[homPol] = TS._evaluate(a.coeffs[homPol], vals)
    end
    return suma
end

function TS._evaluate(a::HomogeneousPolynomial{T},
        vals::NTuple{N,TaylorN{Interval{S}}}) where {N, T<:Real, S<:NumTypes}
    TS._check_same_space_all(a, vals)
    ct = a.space.coeff_table[order(a)+1]
    suma = zero(a[1])*vals[1]
    for (i, a_coeff) in enumerate(a.coeffs)
        TS._isthinzero(a_coeff) && continue
        term = Base.literal_pow(^, vals[1], Val(0))
        @inbounds for j in eachindex(vals)
            exponent = ct[i][j]
            exponent == 0 && continue
            term *= Base.literal_pow(^, vals[j], Val(exponent))
        end
        suma += a_coeff * term
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

aff_normalize(x, I::Interval, ::Val{true})  = interval(mid(I)) + x * interval(radius(I))
aff_normalize(x, I::Interval, ::Val{false}) = interval(inf(I)) + x * interval(diam(I))

for bb in (:true, :false)
    @eval function _normalize(a::Taylor1, I::Interval{T}, ::Val{$bb}) where {T<:NumTypes}
        S = promote_type(TS.numtype(a), Interval{T})
        z = zero(convert(S, constant_term(a)))
        t = Taylor1([z, one(z)], order(a))
        return a(aff_normalize(t, I, Val($bb)))
    end

    @eval function _normalize(a::TaylorN, I::AbstractVector{Interval{T}},
            ::Val{$bb}) where {T<:NumTypes}
        order = TS.order(a)
        S = promote_type(TS.numtype(a), Interval{T})
        x = Vector{TaylorN{S}}(undef, length(I))
        @inbounds for ind in eachindex(x)
            # x[ind] = mid(I[ind]) + TaylorN(ind, order=order)*radius(I[ind])
            x[ind] = aff_normalize(
                    TaylorN(space(a), S, ind, order=order), I[ind], Val($bb))
        end
        aa = convert(TaylorN{S}, a)
        return evaluate(aa, x)
    end
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