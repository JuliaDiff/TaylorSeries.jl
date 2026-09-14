# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

## Evaluating ##
"""
    _evaluate(a::Taylor1, dx::NumberNotSeries)

Evaluate `a` at the ordinary numeric scalar `dx` using Horner's rule without
creating intermediate containers. This is the scalar kernel shared by
`evaluate` and the corresponding array `evaluate!` method.
"""
@inline function _evaluate(a::Taylor1{T}, dx::S) where
        {T<:NumberNotSeries, S<:NumberNotSeries}
    a_coeffs = a.coeffs
    @inbounds suma = zero(a_coeffs[end])*dx
    @inbounds for k in reverse(eachindex(a_coeffs))
        suma = suma * dx + a_coeffs[k]
    end
    return suma
end

"""
    evaluate(a, [dx])

Evaluate a `Taylor1` polynomial using Horner's rule (hand coded). If `dx` is
omitted, its value is considered as zero. Note that the syntax `a(dx)` is
equivalent to `evaluate(a,dx)`, and `a()` is equivalent to `evaluate(a)`.
"""
evaluate(a::Taylor1{T}, dx::S) where
    {T<:NumberNotSeries, S<:NumberNotSeries} = _evaluate(a, dx)

function evaluate(a::Taylor1{T}, dx::S) where {T<:Number, S<:Number}
    a_coeffs = a.coeffs
    suma = a_coeffs[end]*zero(dx)
    @inbounds for k in reverse(eachindex(a_coeffs))
        suma = suma * dx + a_coeffs[k]
    end
    return suma
end

evaluate(a::Taylor1{T}) where {T<:Number} = a.coeffs[1]


"""
    evaluate(x, δt)

Evaluates each element of `x::AbstractArray{Taylor1{T}}`,
representing the dependent variables of an ODE, at *time* δt. Note that the
syntax `x(δt)` is equivalent to `evaluate(x, δt)`, and `x()`
is equivalent to `evaluate(x)`.
"""
# TODO: Preserve a concrete promoted element type for empty arrays without
# changing the container semantics of static arrays.
evaluate(x::AbstractArray{Taylor1{T}}, δt::S) where
    {T<:Number, S<:Number} = evaluate.(x, δt)

evaluate(a::AbstractArray{Taylor1{T}}) where {T<:Number} = getcoeff.(a, 0)


"""
    evaluate(a::Taylor1, x::Taylor1)

Substitute `x::Taylor1` as independent variable in a `a::Taylor1` polynomial.
Note that the syntax `a(x)` is equivalent to `evaluate(a, x)`.
"""
evaluate(a::Taylor1{T}, x::Taylor1{S}) where {T<:Number, S<:Number} =
    evaluate(promote(a, x)...)

function evaluate(a::Taylor1{T}, x::Taylor1{T}) where {T<:NumberNotSeries}
    if order(a) != order(x)
        a, x = fixorder(a, x)
    end
    suma = zero(x)
    aux = zero(x)
    _horner!(suma, a, x, aux)
    return suma
end

function evaluate(a::Taylor1{T}, x::Taylor1{T}) where {T<:Number}
    if order(a) != order(x)
        a, x = fixorder(a, x)
    end
    @inbounds suma = a.coeffs[end]*zero(x)
    _horner!(suma, a, x, zero(suma))
    return suma
end

"""
    _evaluation_order(a, x)

Return the minimum positive order represented by `x` and the coefficients of
`a`, or zero if they all have order zero. Order-zero series represent exact
constants, so they do not truncate a positive-order composition.
"""
function _evaluation_order(a::Taylor1{<:AbstractSeries}, x::AbstractSeries)
    ord = order(x)
    for coeff in a.coeffs
        coeff_order = order(coeff)
        # Only positive orders limit the information available in the result.
        if ord == 0
            ord = coeff_order
        elseif coeff_order > 0
            ord = min(ord, coeff_order)
        end
    end
    return ord
end

function evaluate(a::Taylor1{Taylor1{T}}, x::Taylor1{T}) where {T<:NumberNotSeriesN}
    suma = Taylor1(zero(x[0]), _evaluation_order(a, x))
    _horner!(suma, a, x, zero(suma))
    return suma
end

function evaluate(a::Taylor1{T}, x::Taylor1{Taylor1{T}}) where {T<:NumberNotSeriesN}
    @inbounds suma = a[end]*zero(x)
    _horner!(suma, a, x, zero(suma))
    return suma
end

evaluate(p::Taylor1{T}, x::AbstractArray{S}) where {T<:Number, S<:Number} =
    evaluate.(Ref(p), x)

# Substitute a TaylorN into a Taylor1
function evaluate(a::Taylor1{T}, dx::TaylorN{T}) where {T<:NumberNotSeries}
    suma = TaylorN(dx.space, zero(T), order(dx))
    _horner!(suma, a, dx, zero(suma))
    return suma
end

function evaluate(a::Taylor1{T}, dx::Taylor1{TaylorN{T}}) where
        {T<:NumberNotSeries}
    if order(a) != order(dx)
        a, dx = fixorder(a, dx)
    end
    suma = Taylor1( zero(dx[0]), order(a))
    aux  = zero(suma)
    _horner!(suma, a, dx, aux)
    return suma
end

# Evaluate a Taylor1{TaylorN{T}} on Vector{T} (or Vector{TaylorN{T}}) is interpreted
# as a substitution on the TaylorN vars
function evaluate(a::Taylor1{TaylorN{T}}, dx::AbstractVector{S}) where
        {T<:NumberNotSeries, S<:NumberNotSeries}
    @assert length(dx) == get_numvars(a[0])
    suma = Taylor1( zero(a[0][0][1])*one(dx[1]), order(a))
    suma.coeffs .= evaluate.(a[:], Ref(dx))
    return suma
end

function evaluate(a::Taylor1{TaylorN{T}}, dx::AbstractVector{TaylorN{T}}) where
        {T<:NumberNotSeries}
    @assert length(dx) == get_numvars(a[0])
    _check_same_space(a[0], dx[1])
    suma = Taylor1( zero(a[0]), order(a))
    suma.coeffs .= evaluate.(a[:], Ref(dx))
    return suma
end

function evaluate(a::Taylor1{TaylorN{T}}, ind::Int, dx::T) where
        {T<:NumberNotSeries}
    @assert (1 ≤ ind ≤ get_numvars(a[0])) "Invalid `ind`; it must be between 1 and `get_numvars()`"
    suma = Taylor1( zero(a[0]), order(a))
    for ord in eachindex(suma)
        for ordQ in eachindex(a[0])
            _evaluate!(suma[ord], a[ord][ordQ], ind, dx)
        end
    end
    return suma
end

function evaluate(a::Taylor1{TaylorN{T}}, ind::Int, dx::TaylorN{T}) where
        {T<:NumberNotSeries}
    @assert (1 ≤ ind ≤ get_numvars(a[0])) "Invalid `ind`; it must be between 1 and `get_numvars()`"
    _check_same_space(a[0], dx)
    suma = Taylor1( zero(a[0]), order(a))
    aux = zero(dx)
    for ord in eachindex(suma)
        for ordQ in eachindex(a[0])
            _evaluate!(suma[ord], a[ord][ordQ], ind, dx, aux)
        end
    end
    return suma
end


#function-like behavior for Taylor1
(p::Taylor1)(x) = evaluate(p, x)
(p::Taylor1)()  = evaluate(p)

#function-like behavior for AbstractArray{Taylor1{T}} (asumes Julia version >= 1.6)
(p::AbstractArray{Taylor1{T}})(x) where {T<:Number} = evaluate.(p, x)
(p::AbstractArray{Taylor1{T}})() where {T<:Number} = evaluate.(p)


"""
    evaluate(a, [vals])

Evaluate a `HomogeneousPolynomial` polynomial at `vals`. If `vals` is omitted,
it's evaluated at zero. Note that the syntax `a(vals)` is equivalent to
`evaluate(a, vals)`; and `a()` is equivalent to `evaluate(a)`.
"""
function evaluate(a::HomogeneousPolynomial, vals::NTuple{N,<:Number}) where {N}
    @assert length(vals) == get_numvars(a)
    return _evaluate(a, vals)
end

evaluate(a::HomogeneousPolynomial{T}, vals::AbstractArray{S,1} ) where
    {T<:Number,S<:NumberNotSeriesN} = evaluate(a, (vals...,))

evaluate(a::HomogeneousPolynomial, v, vals::Vararg{Number,N}) where {N} =
    evaluate(a, promote(v, vals...,))

evaluate(a::HomogeneousPolynomial, v) = evaluate(a, promote(v...,))

function evaluate(a::HomogeneousPolynomial{T}) where {T}
    order(a) == 0 && return a[1]
    return zero(a[1])
end

"""
    _evaluate(a::HomogeneousPolynomial, vals)

Evaluate one homogeneous polynomial directly at `vals` without constructing
vectors of powers or monomial terms. The caller must supply a nonempty
collection with one value per variable and, for series values, compatible
series spaces.
"""
function _evaluate(a::HomogeneousPolynomial{T},
        vals) where {T}
    order(a) == 0 && return a[1]*one(vals[1])
    ct = a.space.coeff_table[order(a)+1]
    suma = zero(a[1]*one(vals[1]))
    for (i, a_coeff) in enumerate(a.coeffs)
        TS._isthinzero(a_coeff) && continue
        term = a_coeff * one(vals[1])
        @inbounds for j in eachindex(vals)
            exponent = ct[i][j]
            exponent == 0 && continue
            term *= Base.literal_pow(^, vals[j], Val(exponent))
        end
        suma += term
    end
    return suma
end

# Ordinary scalar vectors can use type-level promotion. The generic method
# above derives zeros and ones from instances to retain series-space metadata.
function _evaluate(a::HomogeneousPolynomial{T},
        vals::AbstractVector{S}) where
        {T<:NumberNotSeries,S<:NumberNotSeries}
    R = promote_type(T, S)
    order(a) == 0 && return convert(R, a[1])
    ct = a.space.coeff_table[order(a)+1]
    suma = zero(R)
    for (i, a_coeff) in enumerate(a.coeffs)
        TS._isthinzero(a_coeff) && continue
        term = convert(R, a_coeff)
        @inbounds for j in eachindex(vals)
            exponent = ct[i][j]
            exponent == 0 && continue
            term *= vals[j]^exponent
        end
        suma += term
    end
    return suma
end

function _evaluate!(res::TaylorN{T}, a::HomogeneousPolynomial{T},
        vals::NTuple{N,<:TaylorN{T}}, valscache::Vector{TaylorN{T}},
        aux::TaylorN{T}) where {N,T<:NumberNotSeries}
    _check_same_space(res, a)
    _check_same_space(a, vals[1])
    ct = a.space.coeff_table[order(a)+1]
    for el in eachindex(valscache)
        power_by_squaring!(valscache[el], vals[el], aux, ct[1][el])
    end
    for (i, a_coeff) in enumerate(a.coeffs)
        TS._isthinzero(a_coeff) && continue
        # valscache .= vals .^ ct[i]
        @inbounds for el in eachindex(valscache)
            power_by_squaring!(valscache[el], vals[el], aux, ct[i][el])
        end
        # aux = one(valscache[1])
        for ord in eachindex(aux)
            @inbounds one!(aux, valscache[1], ord)
        end
        for j in eachindex(valscache)
            # aux *= valscache[j]
            mul!(aux, valscache[j])
        end
        # res += a_coeff * aux
        for ord in eachindex(aux)
            muladd!(res, a_coeff, aux, ord)
        end
    end
    return nothing
end

function _evaluate(a::HomogeneousPolynomial{T},
        vals::NTuple{N,<:TaylorN{T}}) where {N,T<:NumberNotSeries}
    # @assert length(vals) == get_numvars()
    order(a) == 0 && return a[1]*one(vals[1])
    _check_same_space(a, vals[1])
    suma = TaylorN(vals[1].space, zero(T), order(vals[1]))
    valscache = [zero(val) for val in vals]
    aux = zero(suma)
    _evaluate!(suma, a, vals, valscache, aux)
    return suma
end

function _evaluate(a::HomogeneousPolynomial{T}, ind::Int, val::T) where
        {T<:NumberNotSeries}
    suma = TaylorN(a.space, zero(T), order(a.space))
    _evaluate!(suma, a, ind, val)
    return suma
end


#function-like behavior for HomogeneousPolynomial
(p::HomogeneousPolynomial)(x) = evaluate(p, x)
(p::HomogeneousPolynomial)(x, v::Vararg{Number,N}) where {N} =
    evaluate(p, promote(x, v...,))
(p::HomogeneousPolynomial)() = evaluate(p)


"""
    evaluate(a, [vals]; sorting::Bool=true)

Evaluate the `TaylorN` polynomial `a` at `vals`.
If `vals` is omitted, it's evaluated at zero. The
keyword parameter `sorting` can be used to avoid
sorting (in increasing order by `abs2`) the
terms that are added.

Note that the syntax `a(vals)` is equivalent to
`evaluate(a, vals)`; and `a()` is equivalent to
`evaluate(a)`; use a(b::Bool, x) corresponds to
evaluate(a, x, sorting=b).
"""
function evaluate(a::TaylorN, vals::NTuple{N,<:Number};
        sorting::Bool=true) where {N}
    @assert get_numvars(a) == N
    return _evaluate(a, vals, Val(sorting))
end

function evaluate(a::TaylorN, vals::NTuple{N,<:AbstractSeries};
        sorting::Bool=false) where {N}
    @assert get_numvars(a) == N
    return _evaluate(a, vals, Val(sorting))
end

evaluate(a::TaylorN{T}, vals::AbstractVector{<:Number}; sorting::Bool=true) where
    {T<:NumberNotSeries} = evaluate(a, (vals...,); sorting)

evaluate(a::TaylorN{T}, vals::AbstractVector{<:AbstractSeries}; sorting::Bool=false) where
    {T<:NumberNotSeries} = evaluate(a, (vals...,); sorting)

evaluate(a::TaylorN{Taylor1{T}}, vals::AbstractVector{S};
    sorting::Bool=false) where {T, S} = evaluate(a, (vals...,); sorting)

function evaluate(a::TaylorN{T}, s::Symbol, val::S) where
        {T<:Number, S<:NumberNotSeriesN}
    ind = lookupvar(a.space, s)
    @assert (1 ≤ ind ≤ get_numvars(a)) "Symbol is not a TaylorN variable; see `get_variable_names()`"
    return evaluate(a, ind, val)
end

function evaluate(a::TaylorN{T}, ind::Int, val::S) where {T<:Number,
        S<:NumberNotSeriesN}
    @assert (1 ≤ ind ≤ get_numvars(a)) "Invalid `ind`; it must be between 1 and `get_numvars()`"
    R = promote_type(T,S)
    return _evaluate(convert(TaylorN{R}, a), ind, convert(R, val))
end

function evaluate(a::TaylorN{T}, s::Symbol, val::TaylorN) where {T<:Number}
    ind = lookupvar(a.space, s)
    @assert (1 ≤ ind ≤ get_numvars(a)) "Symbol is not a TaylorN variable; see `get_variable_names()`"
    return evaluate(a, ind, val)
end

function evaluate(a::TaylorN{T}, ind::Int, val::TaylorN) where {T<:Number}
    @assert (1 ≤ ind ≤ get_numvars(a)) "Invalid `ind`; it must be between 1 and `get_numvars()`"
    _check_same_space(a, val)
    a, val = fixorder(a, val)
    a, val = promote(a, val)
    return _evaluate(a, ind, val)
end

evaluate(a::TaylorN{T}, x::Pair{Symbol,S}) where {T, S} =
    evaluate(a, first(x), last(x))

evaluate(a::TaylorN{T}) where {T<:Number} = constant_term(a)

"""
    _evaluate(a::TaylorN, vals::Tuple, ::Val{true})
    _evaluate(a::TaylorN, vals::Tuple, ::Val{false})

Evaluate `a` at `vals` and return the sum of its contributions.
`Val(true)` sorts the contributions from smallest to largest magnitude before
adding them, to reduce rounding error. `Val(false)` skips this sorting step.

The third argument is needed because `_evaluate(a::TaylorN, vals::Tuple)`
returns an array with one contribution per polynomial degree, rather than
their sum. The sorted method needs those separate values so it can reorder
them before adding them. Using `Val(true)` or `Val(false)` selects the method
for the requested summation without changing the two-argument method's result.

Sorting currently allocates temporary storage. Reusing a caller-provided
buffer could avoid that allocation while keeping the same sorting behavior.
"""
_evaluate(a::TaylorN{T}, vals::NTuple, ::Val{true}) where
    {T<:NumberNotSeries} = sum( sort!(_evaluate(a, vals), by=abs2) )

_evaluate(a::TaylorN{T}, vals::NTuple, ::Val{false}) where {T<:Number} =
    sum( _evaluate(a, vals) )

function _evaluate(a::TaylorN{T}, vals::NTuple{N,<:TaylorN}, ::Val{false}) where
        {N,T<:Number}
    R = promote_type(T, TS.numtype(vals[1]))
    _check_same_space(a, vals[1])
    res = TaylorN(vals[1].space, zero(R), order(vals[1]))
    vvals = ntuple(i -> convert(TaylorN{R}, vals[i]), length(vals))
    valscache = [zero(val) for val in vvals]
    aux = zero(res)
    @inbounds for homPol in eachindex(a)
        _evaluate!(res, a[homPol], vvals, valscache, aux)
    end
    return res
end

"""
    _evaluate(a::TaylorN, vals::Tuple)

Return an array containing the evaluated contribution from each polynomial
degree of `a`, starting with its constant term. The entries have not been
added together. For the complete result, use the three-argument form with
`Val(true)` to sort before adding, or `Val(false)` to add without sorting.
"""
function _evaluate(a::TaylorN{T}, vals::NTuple{N,<:Number}) where {N,T<:Number}
    R = promote_type(T, typeof(vals[1]))
    suma = zeros(R, length(a))
    @inbounds for homPol in eachindex(a)
        suma[homPol+1] = _evaluate(a[homPol], vals)
    end
    return suma
end

"""
    _evaluate(a::TaylorN, vals::AbstractVector, ::Val{false})

Evaluate `a` at `vals` and add the contributions in increasing polynomial
degree, without sorting by their values. This vector method adds each
contribution directly, without allocating an array to hold them first.

The third argument distinguishes this complete result from the two-argument
`_evaluate(a::TaylorN, vals::Tuple)`, which returns the separate contributions.
"""
function _evaluate(a::TaylorN{T},
        vals::AbstractVector{S}, ::Val{false}) where {T<:Number, S<:Number}
    @assert length(vals) == get_numvars(a)
    suma = zero(a[0][1] * one(vals[1]))
    @inbounds for homPol in eachindex(a)
        suma += _evaluate(a[homPol], vals)
    end
    return suma
end

function _evaluate(a::TaylorN{T}, vals::NTuple{N,<:TaylorN}) where {N,T<:Number}
    R = promote_type(T, TS.numtype(vals[1]))
    _check_same_space(a, vals[1])
    suma = [TaylorN(vals[1].space, zero(R), order(vals[1])) for _ in eachindex(a)]
    valscache = [zero(val) for val in vals]
    aux = zero(suma[1])
    _evaluate!(suma, a, vals, valscache, aux)
    return suma
end


function _evaluate(a::TaylorN{T}, ind::Int, val::T) where {T<:NumberNotSeriesN}
    suma = TaylorN(a.space, zero(a[0]*val), order(a))
    vval = convert(numtype(suma), val)
    suma, a = promote(suma, a)
    @inbounds for ordQ in eachindex(a)
        _evaluate!(suma, a[ordQ], ind, vval)
    end
    return suma
end

function _evaluate(a::TaylorN{T}, ind::Int, val::TaylorN{T}) where
        {T<:NumberNotSeriesN}
    _check_same_space(a, val)
    suma = TaylorN(a.space, zero(a[0]), order(a))
    aux = zero(suma)
    @inbounds for ordQ in eachindex(a)
        _evaluate!(suma, a[ordQ], ind, val, aux)
    end
    return suma
end

function _evaluate!(res::Vector{TaylorN{T}}, a::TaylorN{T},
        vals::NTuple{N,<:TaylorN}, valscache::Vector{TaylorN{T}},
        aux::TaylorN{T}) where {N,T<:Number}
    @inbounds for homPol in eachindex(a)
        _evaluate!(res[homPol+1], a[homPol], vals, valscache, aux)
    end
    return nothing
end

function _evaluate!(suma::TaylorN{T}, a::HomogeneousPolynomial{T}, ind::Int,
        val::T) where {T<:NumberNotSeriesN}
    _check_same_space(suma, a)
    order = TS.order(a)
    orderTN = TS.order(a.space)
    if order == 0
        suma[0][1] = a[1]*one(val)
        return nothing
    end
    vv = Base.literal_pow.(^, val, Val.(0:order))
    vct = zero(a.space.coeff_table[order+1][1])
    zct = zero(a.space.coeff_table[order+1][1])
    for (i, a_coeff) in enumerate(a.coeffs)
        iszero(a_coeff) && continue
        vpow = a.space.coeff_table[order+1][i][ind]
        if vpow == 0
            suma[order][i] += a_coeff
            continue
        end
        vct .= a.space.coeff_table[order+1][i]
        zct[ind] = vpow
        red_order = order - vpow
        kdic = in_base(orderTN, vct - zct)
        zct[ind] = 0
        pos = a.space.pos_table[red_order+1][kdic]
        suma[red_order][pos] += a_coeff * vv[vpow+1]
    end
    return nothing
end

function _evaluate!(suma::TaylorN{T}, a::HomogeneousPolynomial{T}, ind::Int,
        val::TaylorN{T}, aux::TaylorN{T}) where {T<:NumberNotSeriesN}
    _check_same_space(suma, a, val)
    order = TS.order(a)
    if order == 0
        suma[0][1] = a[1]
        return nothing
    end
    vv = zero(suma)
    vvaux = zero(vv)
    za = zero(a)
    for (i, a_coeff) in enumerate(a.coeffs)
        iszero(a_coeff) && continue
        vpow = a.space.coeff_table[order+1][i][ind]
        if vpow == 0
            suma[order][i] += a_coeff
            continue
        end
        # vv = val ^ vpow
        if constant_term(val) == 0
            zero!(vvaux)
            power_by_squaring!(vv, val, vvaux, vpow)
        else
            for ordQ in eachindex(val)
                zero!(vv, ordQ)
                pow!(vv, val, vvaux, vpow, ordQ)
            end
        end
        za[i] = a_coeff
        zero!(aux)
        _evaluate!(aux, za, ind, one(T))
        za[i] = zero(a_coeff)
        for ordQ in eachindex(suma)
            mul!(suma, vv, aux, ordQ)
        end
    end
    return nothing
end


# High-dimensional array evaluation. Keep this allocating interface as a
# broadcast of scalar evaluations so it retains scalar promotion, sorting,
# shape, and StaticArray container behavior while accepting views as inputs.
# TODO: Preserve concrete result element types for empty arrays once this can
# be done without duplicating the scalar promotion rules.
function evaluate(A::AbstractArray{TaylorN{T}}, vals::AbstractVector{S};
        sorting::Bool=!(T <: AbstractSeries || S <: AbstractSeries)) where
        {T<:Number,S<:Number}
    return evaluate.(A, Ref(vals); sorting)
end

function evaluate(A::AbstractArray{TaylorN{T}}, vals::Tuple{S,Vararg{S}};
        sorting::Bool=!(T <: AbstractSeries || S <: AbstractSeries)) where
        {T<:Number,S<:Number}
    return evaluate.(A, Ref(vals); sorting)
end

evaluate(A::AbstractArray{TaylorN{T}}) where {T<:Number} = evaluate.(A)

#function-like behavior for TaylorN
(p::TaylorN)(x) = evaluate(p, x)
(p::TaylorN)() = evaluate(p)
(p::TaylorN)(s::S, x) where {S<:Union{Symbol, Int}} = evaluate(p, s, x)
(p::TaylorN)(x::Pair) = evaluate(p, first(x), last(x))
(p::TaylorN)(x, v::Vararg{T}) where {T} = evaluate(p, (x, v...,))
(p::TaylorN)(b::Bool, x) = evaluate(p, x, sorting=b)
(p::TaylorN)(b::Bool, x, v::Vararg{T}) where {T} = evaluate(p, (x, v...,), sorting=b)

#function-like behavior for AbstractArray{TaylorN{T}}
(p::AbstractArray{TaylorN{T}})(x) where {T<:Number} = evaluate(p, x)
(p::AbstractArray{TaylorN{T}})() where {T<:Number} = evaluate(p)


"""
    evaluate!(x, δt, dest)
    evaluate!(x, vals, dest; sorting=...)

Evaluate a polynomial or each polynomial in `x` and write the results into
`dest`. For ordinary numeric `Taylor1` inputs, evaluation uses a direct Horner
kernel. For `TaylorN` arrays, `sorting=false` sums homogeneous components in
stored order without a temporary component vector; `sorting=true` preserves
magnitude-sorted scalar evaluation and may allocate. Sorting defaults to true
for ordinary scalar types and false when the coefficient or evaluation-value
type is a series.

For series-valued substitutions (i.e., when the evaluation argument is
a Taylor series variable), the `evaluate!` methods that accept scratch buffers
should be used, so that buffers are allocated only once and are reused across
calls.
"""
function evaluate!(x::AbstractArray{Taylor1{T}}, δt::S,
        dest::AbstractArray{R}) where
        {T<:NumberNotSeries, S<:NumberNotSeries, R<:NumberNotSeries}
    @inbounds for i in eachindex(x, dest)
        dest[i] = _evaluate(x[i], δt)
    end
    return nothing
end

function evaluate!(a::Taylor1{Taylor1{T}}, δt::S,
        dest::Taylor1{R}) where
        {T<:NumberNotSeries,S<:NumberNotSeries,R<:NumberNotSeries}
    # A scalar evaluation value lets every inner coefficient be updated in
    # place, without a second series-valued Horner buffer.
    for coeff in a.coeffs
        # Order-zero series are exact constants and may be zero-extended.
        (iszero(order(coeff)) || order(dest) <= order(coeff)) || throw(DimensionMismatch(
            "destination order exceeds a Taylor1 coefficient order"))
    end
    zero!(dest)
    @inbounds for k in reverse(eachindex(a))
        coeff = a[k]
        for ord in eachindex(dest)
            dest[ord] *= δt
            ord <= order(coeff) && (dest[ord] += coeff[ord])
        end
    end
    return nothing
end

function evaluate!(x::AbstractArray{Taylor1{Taylor1{T}}}, δt::S,
        dest::AbstractArray{Taylor1{R}}) where
        {T<:NumberNotSeries,S<:NumberNotSeries,R<:NumberNotSeries}
    @inbounds for i in eachindex(x, dest)
        evaluate!(x[i], δt, dest[i])
    end
    return nothing
end

function evaluate!(a::Taylor1{TaylorN{T}}, δt::S,
        dest::TaylorN{T}) where {T<:Number, S<:NumberNotSeries}
    # As above, scalar multiplication can update each homogeneous component
    # directly, so this path needs no separate outer TaylorN scratch buffer.
    for coeff in a.coeffs
        _check_same_space(dest, coeff)
        # Order-zero series are exact constants and may be zero-extended.
        (iszero(order(coeff)) || order(dest) <= order(coeff)) ||
            throw(DimensionMismatch(
                "destination order exceeds a TaylorN coefficient order"))
    end
    zero!(dest)
    @inbounds for k in reverse(eachindex(a))
        a_coeff = a[k]
        for ordQ in eachindex(dest)
            dest_hp = dest[ordQ].coeffs
            if ordQ <= order(a_coeff)
                a_hp = a_coeff[ordQ].coeffs
                for j in eachindex(dest_hp)
                    dest_hp[j] = dest_hp[j] * δt + a_hp[j]
                end
            else
                for j in eachindex(dest_hp)
                    dest_hp[j] *= δt
                end
            end
        end
    end
    return nothing
end

function evaluate!(x::AbstractArray{Taylor1{TaylorN{T}}}, δt::S,
        dest::AbstractArray{TaylorN{T}}) where {T<:Number, S<:NumberNotSeries}
    @inbounds for i in eachindex(x, dest)
        evaluate!(x[i], δt, dest[i])
    end
    return nothing
end

"""
    _check_series_evaluation(a, δt, dest, aux)

Validate the series-space, order, and non-aliasing requirements of
series-valued `Taylor1` evaluation with explicit workspace. All checks happen
before the destination is mutated.
"""
function _check_series_evaluation(a::Taylor1{T}, δt::T, dest::T,
        aux::T) where {T<:Union{Taylor1,TaylorN}}
    # Horner multiplication requires independent input, output, and scratch
    # storage from compatible series families (and the same TaylorN JetSpace).
    dest === δt && throw(ArgumentError("destination must not alias the evaluation value"))
    dest === aux && throw(ArgumentError("destination and scratch value must not alias"))
    δt === aux && throw(ArgumentError("evaluation and scratch values must not alias"))
    _check_same_space(dest, δt, aux)
    order(dest) == order(aux) || throw(DimensionMismatch(
        "destination and scratch value must have the same order"))
    order(dest) <= order(δt) || throw(DimensionMismatch(
        "destination order exceeds the evaluation value order"))
    for coeff in a.coeffs
        _check_same_space(dest, coeff)
        # Coefficients are read throughout Horner evaluation and must not be
        # overwritten through either mutable work buffer.
        coeff === dest && throw(ArgumentError(
            "destination must not alias a polynomial coefficient"))
        coeff === aux && throw(ArgumentError(
            "scratch value must not alias a polynomial coefficient"))
        # Order-zero coefficients are exact constants; only positive orders
        # constrain the amount of information available in the destination.
        (iszero(order(coeff)) || order(dest) <= order(coeff)) ||
            throw(DimensionMismatch(
                "destination order exceeds a nonconstant coefficient order"))
    end
    return nothing
end

"""
    evaluate!(a, δt, dest, aux)
    evaluate!(x, δt, dest, aux)

Evaluate a `Taylor1` polynomial, or an array of them, at a series-valued
`δt`. The result is written into `dest`, while `aux` is reusable scratch
storage. These explicit-workspace overloads do not create scratch storage.
The destination, evaluation value and scratch value must have the same
concrete series type and, for `TaylorN`, the same `JetSpace`. They must not
alias one another; neither `dest` nor `aux` may alias a mutable coefficient.
`dest` and `aux` must have the same order, no greater than the order of `δt`.
Nonconstant coefficient orders must also be at least the destination order.
"""
function evaluate!(a::Taylor1{T}, δt::T, dest::T,
        aux::T) where {T<:Union{Taylor1,TaylorN}}
    _check_series_evaluation(a, δt, dest, aux)
    zero!(dest)
    _horner!(dest, a, δt, aux)
    return nothing
end

function evaluate!(x::AbstractArray{Taylor1{T}}, δt::T,
        dest::AbstractArray{T}, aux::T) where {T<:Union{Taylor1,TaylorN}}
    @inbounds for i in eachindex(x, dest)
        evaluate!(x[i], δt, dest[i], aux)
    end
    return nothing
end

function evaluate!(x::AbstractArray{Taylor1{T}}, δt::T,
        dest::AbstractArray{T}) where {T<:Union{Taylor1,TaylorN}}
    if isempty(dest)
        isempty(x) && return nothing
        throw(DimensionMismatch("source and destination arrays must have matching indices"))
    end
    # Compatibility method: the four-argument `evaluate!(x, δt, dest, aux)`
    # method does not allocate `aux`; this form creates one auxiliary variable
    # `aux` series which is allocated only once and reused throughout.
    # TODO: Detect heterogeneous destination orders and create correctly sized
    # per-element scratch only for that fallback.
    aux = zero(dest[firstindex(dest)])
    evaluate!(x, δt, dest, aux)
    return nothing
end

## In place evaluation of multivariable arrays
"""
    _evaluate!(x, vals, dest, ::Val{sorting})

Evaluate each polynomial in `x` at `vals` and write its result into `dest`.
`Val(false)` adds contributions directly without a temporary result array.
`Val(true)` sorts each polynomial's contributions by magnitude before adding
them, and may allocate temporary storage for that sorting.
"""
function _evaluate!(x::AbstractArray{TaylorN{T}},
        δx::AbstractVector{S}, dest::AbstractArray{R},
        ::Val{false}) where
        {T<:Number,S<:Number,R<:Number}
    @inbounds for i in eachindex(x, dest)
        dest[i] = _evaluate(x[i], δx, Val(false))
    end
    return nothing
end

function _evaluate!(x::AbstractArray{TaylorN{T}},
        δx::AbstractVector{S}, dest::AbstractArray{R},
        ::Val{true}) where
        {T<:Number,S<:Number,R<:Number}
    @inbounds for i in eachindex(x, dest)
        dest[i] = evaluate(x[i], δx; sorting=true)
    end
    return nothing
end

@inline function evaluate!(x::AbstractArray{TaylorN{T}}, δx::AbstractVector{S},
        dest::AbstractArray{R};
        sorting::Bool=!(T <: AbstractSeries || S <: AbstractSeries)) where
        {T<:Number,S<:Number,R<:Number}
    _evaluate!(x, δx, dest, Val(sorting))
    return nothing
end

function evaluate!(x::AbstractArray{TaylorN{T}}, δx::AbstractVector{TaylorN{T}},
        dest::AbstractArray{TaylorN{T}}; sorting::Bool=false) where {T<:NumberNotSeriesN}
    if sorting
        # Sorted evaluation intentionally follows the scalar allocating path
        # to preserve its summation order and numerical behavior.
        @inbounds for i in eachindex(x, dest)
            dest[i] = evaluate(x[i], δx; sorting=true)
        end
    else
        evaluate!(x, (δx...,), dest; sorting=false)
    end
    return nothing
end

function _evaluate!(a::TaylorN{T}, vals::NTuple{N,TaylorN{T}}, dest::TaylorN{T},
        valscache::Vector{TaylorN{T}}, aux::TaylorN{T}) where {N,T<:Number}
    @inbounds for homPol in eachindex(a)
        _evaluate!(dest, a[homPol], vals, valscache, aux)
    end
    return nothing
end

"""
    _check_taylorN_evaluation(a, vals, dest, valscache, aux)

Validate dimensions, `JetSpace` and order compatibility, and non-aliasing
requirements for `TaylorN` evaluation with explicit workspace. Cache entries
are destructive scratch and must be distinct from every input and from one
another. All checks happen before the destination is mutated.
"""
function _check_taylorN_evaluation(a::TaylorN{T}, vals::NTuple{N,TaylorN{T}},
        dest::TaylorN{T}, valscache::Vector{TaylorN{T}},
        aux::TaylorN{T}) where {N,T<:Number}
    # Validate dimensions, common series space, and truncation order first.
    get_numvars(a) == N || throw(DimensionMismatch(
        "number of evaluation values must match the number of variables"))
    length(valscache) == N || throw(DimensionMismatch(
        "evaluation cache length must match the number of evaluation values"))
    _check_same_space(a, dest, aux)
    order(dest) == order(aux) || throw(DimensionMismatch(
        "destination and scratch value must have the same order"))
    # Destination and auxiliary storage are mutated and may not alias each
    # other or the source polynomial.
    a === dest && throw(ArgumentError(
        "destination must not alias the polynomial being evaluated"))
    a === aux && throw(ArgumentError(
        "scratch value must not alias the polynomial being evaluated"))
    dest === aux && throw(ArgumentError("destination and scratch value must not alias"))
    for i in eachindex(vals)
        _check_same_space(dest, vals[i], valscache[i])
        order(vals[i]) == order(dest) == order(valscache[i]) ||
            throw(DimensionMismatch(
                "evaluation values, destination and cache must have the same order"))
        vals[i] === dest && throw(ArgumentError(
            "destination must not alias an evaluation value"))
        vals[i] === aux && throw(ArgumentError(
            "scratch value must not alias an evaluation value"))
        valscache[i] === dest && throw(ArgumentError(
            "destination must not alias the evaluation cache"))
        valscache[i] === aux && throw(ArgumentError(
            "scratch value must not alias the evaluation cache"))
        valscache[i] === a && throw(ArgumentError(
            "evaluation cache must not alias the polynomial being evaluated"))
        # Each cache entry is repeatedly overwritten with powers of one input.
        for val in vals
            valscache[i] === val && throw(ArgumentError(
                "evaluation cache must not alias an evaluation value"))
        end
        for j in firstindex(valscache):(i-1)
            valscache[i] === valscache[j] && throw(ArgumentError(
                "evaluation cache entries must not alias one another"))
        end
    end
    return nothing
end

"""
    evaluate!(a, vals, dest, valscache, aux; sorting=false)
    evaluate!(a_array, vals, dest_array, valscache, aux; sorting=false)

Evaluate a `TaylorN` polynomial, or an array of them, at series-valued `vals`
and write into `dest`. `valscache` provides one destructive scratch series per
evaluation value and `aux` provides shared arithmetic scratch. Every evaluation
value and workspace series must have the same `JetSpace` and order as `dest`;
the source polynomial must share that `JetSpace` but may have a different
order. Workspace must not alias inputs, outputs, or other cache entries.

With `sorting=false`, these overloads reuse the supplied workspace without
creating evaluation scratch. `sorting=true` preserves magnitude-sorted scalar
evaluation and may allocate an intermediate result.
"""
function evaluate!(a::TaylorN{T}, vals::NTuple{N,TaylorN{T}},
        dest::TaylorN{T}, valscache::Vector{TaylorN{T}},
        aux::TaylorN{T}; sorting::Bool=false) where {N,T<:Number}
    _check_taylorN_evaluation(a, vals, dest, valscache, aux)
    if sorting
        result = evaluate(a, vals; sorting=true)
        zero!(dest)
        for ord in eachindex(dest)
            identity!(dest, result, ord)
        end
    else
        zero!(dest)
        _evaluate!(a, vals, dest, valscache, aux)
    end
    return nothing
end

"""
    evaluate!(a, vals, dest, valscache, aux; sorting=false)

Array form of explicit-workspace `TaylorN` evaluation. The cache and auxiliary
series are reused sequentially for every element of `a`; their requirements
and sorting behavior are the same as for the scalar method above.
"""
function evaluate!(a::AbstractArray{TaylorN{T}}, vals::NTuple{N,TaylorN{T}},
        dest::AbstractArray{TaylorN{T}}, valscache::Vector{TaylorN{T}},
        aux::TaylorN{T}; sorting::Bool=false) where {N,T<:Number}
    for i in eachindex(a, dest)
        evaluate!(a[i], vals, dest[i], valscache, aux; sorting)
    end
    return nothing
end

function evaluate!(a::AbstractArray{TaylorN{T}}, vals::NTuple{N,TaylorN{T}},
        dest::AbstractArray{TaylorN{T}}; sorting::Bool=false) where {N,T<:Number}
    if isempty(dest)
        isempty(a) && return nothing
        throw(DimensionMismatch("source and destination arrays must have matching indices"))
    end
    # Compatibility method: construct one cache set and auxiliary series for
    # this call. The explicit-workspace method 
    # `evaluate!(a, vals, dest, valscache, aux; sorting)` is the allocation-free method.
    valscache = [zero(val) for val in vals]
    aux = zero(dest[firstindex(dest)])
    evaluate!(a, vals, dest, valscache, aux; sorting)
    return nothing
end

# In-place Horner kernels for series-valued evaluation and composition.
function _horner!(suma::Taylor1{T}, a::Taylor1{T}, x::Taylor1{T},
        aux::Taylor1{T}) where {T<:Number}
    @inbounds for k in reverse(eachindex(a))
        for ord in eachindex(suma)
            mul!(aux, suma, x, ord)
        end
        for ord in eachindex(suma)
            add!(suma, aux, a[k], ord)
        end
    end
    return nothing
end

function _horner!(suma::Taylor1{T}, a::Taylor1{Taylor1{T}}, x::Taylor1{T},
        aux::Taylor1{T}) where {T<:Number}
    @inbounds for k in reverse(eachindex(a))
        for ord in eachindex(suma)
            mul!(aux, suma, x, ord)
        end
        for ord in eachindex(suma)
            # An order-zero coefficient is exact; at higher result orders only
            # the Horner product already stored in `aux` contributes.
            if ord <= order(a[k])
                add!(suma, aux, a[k], ord)
            else
                identity!(suma, aux, ord)
            end
        end
    end
    return nothing
end

function _horner!(suma::Taylor1{Taylor1{T}}, a::Taylor1{T}, x::Taylor1{Taylor1{T}},
        aux::Taylor1{Taylor1{T}}) where {T<:Number}
    @inbounds for k in reverse(eachindex(a))
        for ord in eachindex(suma)
            mul!(aux, suma, x, ord)
        end
        for ord in eachindex(suma)
            identity!(suma, aux, ord)
            add!(suma, suma, a[k], ord)
        end
    end
    return nothing
end

function _horner!(suma::TaylorN{T}, a::Taylor1{T}, dx::TaylorN{T},
        aux::TaylorN{T})  where {T<:NumberNotSeries}
    @inbounds for k in reverse(eachindex(a))
        for ordQ in eachindex(suma)
            zero!(aux, ordQ)
            mul!(aux, suma, dx, ordQ)
        end
        for ordQ in eachindex(suma)
            identity!(suma, aux, ordQ)
        end
        add!(suma, suma, a[k], 0)
    end
    return nothing
end

function _horner!(suma::Taylor1{TaylorN{T}}, a::Taylor1{T}, dx::Taylor1{TaylorN{T}},
        aux::Taylor1{TaylorN{T}})  where {T<:NumberNotSeries}
    @inbounds for k in reverse(eachindex(a))
        for ord in eachindex(suma)
            zero!(aux, ord)
            mul!(aux, suma, dx, ord)
        end
        for ord in eachindex(aux)
            add!(suma, aux, a[k], ord)
        end
    end
    return suma
end

function _horner!(suma::TaylorN{T}, a::Taylor1{TaylorN{T}}, dx::S,
        aux::TaylorN{T}) where {T<:Number, S<:Number}
    @inbounds for k in reverse(eachindex(a))
        for ord in eachindex(suma)
            zero!(aux, ord)
            mul!(aux, suma, dx, ord)
        end
        for ord in eachindex(aux)
            # Order-zero coefficients are exact constants; at higher orders
            # only the Horner product contributes.
            if ord <= order(a[k])
                add!(suma, aux, a[k], ord)
            else
                identity!(suma, aux, ord)
            end
        end
    end
    return nothing
end
