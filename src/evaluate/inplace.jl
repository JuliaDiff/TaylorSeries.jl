# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#


# Inplace methods and its auxiliaries
"""
    evaluate!(x, δt, dest)
    evaluate!(x, vals, dest; sorting=...)
    evaluate!(x, δt, dest, aux)
    evaluate!(x, vals, dest, valscache, aux; sorting=false)

Evaluate a polynomial, or an array of polynomials, and write the result into
`dest`, a pre-allocated destination. All these methods return `nothing`.

For `TaylorN` evaluation, `sorting=true` sorts the contributions by magnitude
before adding them, to reduce rounding error. `sorting=false` skips sorting.
The default, given by `_defaultsorting`, is `true` when both the coefficients
and evaluation values are usual numbers (other than Taylor series or intervals),
and `false` otherwise. Passing pre-allocated auxiliaries does not disable
sorting; choose `sorting=false` to avoid the temporary result used by sorted
evaluation.

The methods for a single polynomial are also used by the corresponding array
methods. The allocating `evaluate` methods do not all call `evaluate!`; some
call `evaluate` separately for each array element.

# Examples

Evaluate two polynomials at a number:

```jldoctest
julia> using TaylorSeries

julia> t = Taylor1(3);

julia> polys = [1 + 2t, 3 + 4t];

julia> dest = zeros(2);

julia> evaluate!(polys, 0.5, dest);

julia> dest == [2.0, 5.0]
true
```

Evaluate a multivariable polynomial at Taylor series. Allocate the destination
and auxiliaries once; they can be reused in subsequent calls. The docstrings
below explain which orders and `JetSpace`s are allowed, and which objects
must be kept separate.

```jldoctest
julia> using TaylorSeries

julia> space = JetSpace(order=2, variables=[:x, :y]);

julia> x, y = variables(space);

julia> vals = (x + 1, y + 2);

julia> dest = zero(x); valscache = [zero(x), zero(y)]; aux = zero(x);

julia> evaluate!(x + y^2, vals, dest, valscache, aux);

julia> dest == (x + 1) + (y + 2)^2
true
```
"""
function evaluate!(x::AbstractArray{Taylor1{T}}, δt::S,
        dest::AbstractArray{R}) where
        {T<:NumberNotSeries, S<:NumberNotSeries, R<:NumberNotSeries}
    @inbounds for i in eachindex(x, dest)
        dest[i] = _evaluate(x[i], δt)
    end
    return nothing
end

function evaluate!(x::AbstractArray{Taylor1{T}}, δt::S,
        dest::AbstractArray{R}) where
        {T<:Number, S<:NumberNotSeries, R<:Number}
    dest .= evaluate.(x, δt)
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
    # directly, so this path needs no separate outer TaylorN buffer.
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
    evaluate!(a, δt, dest, aux)
    evaluate!(x, δt, dest, aux)

Evaluate a `Taylor1` polynomial, or an array of them, at a Taylor series `δt`.
Write the result into `dest` and reuse the pre-allocated auxiliary `aux` for
the calculations. The method does not create `aux` on each call.

The coefficients of `a`, `δt`, `dest`, and `aux` must have the same series
type, including their coefficient type, and the same `JetSpace` for `TaylorN`.
Keep `δt`, `dest`, and `aux` separate from one another. Neither `dest` nor
`aux` may share coefficient storage with a coefficient of `a`.
`dest` and `aux` must have the same order, no greater than the order of `δt`.
Each coefficient of `a` must have at least that order, except that order-zero
coefficients are treated as exact constants.
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

@inline function evaluate!(x::AbstractArray{TaylorN{T}}, δx::AbstractVector{S},
        dest::AbstractArray{R};
        sorting::Bool=_defaultsorting(T,S)) where
        {T<:Number,S<:Number,R<:Number}
    _evaluate!(x, δx, dest, Val(sorting))
    return nothing
end

function evaluate!(x::AbstractArray{TaylorN{T}}, δx::AbstractVector{TaylorN{T}},
        dest::AbstractArray{TaylorN{T}};
        sorting::Bool=false) where {T<:NumberNotSeriesN}
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


"""
    evaluate!(a, vals, dest, valscache, aux; sorting=false)
    evaluate!(a_array, vals, dest_array, valscache, aux; sorting=false)

Evaluate a `TaylorN` polynomial, or an array of them, at the Taylor series in
`vals` and write into `dest`. `valscache` provides one pre-allocated auxiliary per
evaluation value, and `aux` is another pre-allocated auxiliary used for the
arithmetic. Both are overwritten during evaluation. Every evaluation value
and auxiliary must have the same `JetSpace` and order as `dest`; the source
polynomial must share that `JetSpace` but may have a different order.
Keep each entry of `valscache` and `aux` separate from the inputs, destination,
and other auxiliaries. The destination must also be separate from the inputs.

With `sorting=false`, these methods reuse the supplied `valscache` and `aux`.
With `sorting=true`, each polynomial is evaluated with its contributions
sorted by magnitude before addition, and the result is copied into `dest`.
This sorted evaluation may allocate an intermediate result.
"""
function evaluate!(a::TaylorN{T}, vals::NTuple{N,TaylorN{T}},
        dest::TaylorN{T}, valscache::Vector{TaylorN{T}},
        aux::TaylorN{T}; sorting::Bool=false) where {N,T<:Number}
    _check_taylorN_evaluation(a, vals, dest, valscache, aux)
    _evaluate!(a, vals, dest, valscache, aux, Val(sorting))
    return nothing
end

"""
    evaluate!(a, vals, dest, valscache, aux; sorting=false)

Evaluate each `TaylorN` polynomial in `a`, writing the results into `dest`.
Reuse the same pre-allocated auxiliaries `valscache` and `aux` for each
polynomial in turn. Their requirements and the sorting behavior are the same
as for evaluating a single polynomial.

Check the shared evaluation values and pre-allocated auxiliaries once per call.
Then check every source polynomial and destination before evaluating any of them.
"""
function evaluate!(a::AbstractArray{TaylorN{T}}, vals::NTuple{N,TaylorN{T}},
        dest::AbstractArray{TaylorN{T}}, valscache::Vector{TaylorN{T}},
        aux::TaylorN{T}; sorting::Bool=false) where {N,T<:Number}
    indices = eachindex(a, dest)
    isempty(indices) && return nothing
    _check_taylorN_auxiliaries(vals, valscache, aux)
    for i in indices
        _check_taylorN_destination(a[i], vals, dest[i], valscache, aux)
    end
    for i in indices
        _evaluate!(a[i], vals, dest[i], valscache, aux, Val(sorting))
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
    # this call. The method with explicit auxiliary (aux) and cache (valcache)
    # `evaluate!(a, vals, dest, valscache, aux; sorting)` is the allocation-free method.
    valscache = [zero(val) for val in vals]
    aux = zero(dest[firstindex(dest)])
    evaluate!(a, vals, dest, valscache, aux; sorting)
    return nothing
end


"""
    _check_taylorN_evaluation(a, vals, dest, valscache, aux)

Check the number of evaluation values and the orders and `JetSpace`s of
the polynomial, destination, values, and pre-allocated auxiliaries.
Also check that the destination and auxiliaries are separate from the inputs
and from one another, since evaluation overwrites them. These checks run
before changing `dest`.
"""
function _check_taylorN_evaluation(a::TaylorN{T}, vals::NTuple{N,TaylorN{T}},
        dest::TaylorN{T}, valscache::Vector{TaylorN{T}},
        aux::TaylorN{T}) where {N,T<:Number}
    _check_taylorN_auxiliaries(vals, valscache, aux)
    _check_taylorN_destination(a, vals, dest, valscache, aux)
    return nothing
end


"""
    _check_series_evaluation(a, δt, dest, aux)

Check that the evaluation value, destination, and pre-allocated auxiliary
`aux` have compatible orders and, for `TaylorN`, the same `JetSpace`.
Also check that updating the destination or auxiliary will not overwrite
an input or the other output. These checks run before changing `dest`.
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

# These checks depend only on the values and pre-allocated auxiliaries, so
# an array evaluation can perform them once for all its polynomials.
function _check_taylorN_auxiliaries(vals::NTuple{N,TaylorN{T}},
        valscache::Vector{TaylorN{T}}, aux::TaylorN{T}) where {N,T<:Number}
    length(valscache) == N || throw(DimensionMismatch(
        "evaluation cache length must match the number of evaluation values"))
    for i in eachindex(vals)
        _check_same_space(aux, vals[i], valscache[i])
        order(vals[i]) == order(aux) == order(valscache[i]) ||
            throw(DimensionMismatch(
                "evaluation values and pre-allocated auxiliaries must have the same order"))
        vals[i] === aux && throw(ArgumentError(
            "auxiliary must not alias an evaluation value"))
        valscache[i] === aux && throw(ArgumentError(
            "auxiliary must not alias the evaluation cache"))
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

# Check the source and destination against the already checked shared values
# and pre-allocated auxiliaries. The source may have a different order.
function _check_taylorN_destination(a::TaylorN{T}, vals::NTuple{N,TaylorN{T}},
        dest::TaylorN{T}, valscache::Vector{TaylorN{T}},
        aux::TaylorN{T}) where {N,T<:Number}
    get_numvars(a) == N || throw(DimensionMismatch(
        "number of evaluation values must match the number of variables"))
    _check_same_space(a, dest, aux)
    order(dest) == order(aux) || throw(DimensionMismatch(
        "destination and auxiliary must have the same order"))
    a === dest && throw(ArgumentError(
        "destination must not alias the polynomial being evaluated"))
    a === aux && throw(ArgumentError(
        "auxiliary must not alias the polynomial being evaluated"))
    dest === aux && throw(ArgumentError("destination and auxiliary must not alias"))
    for i in eachindex(vals)
        vals[i] === dest && throw(ArgumentError(
            "destination must not alias an evaluation value"))
        valscache[i] === dest && throw(ArgumentError(
            "destination must not alias the evaluation cache"))
        valscache[i] === a && throw(ArgumentError(
            "evaluation cache must not alias the polynomial being evaluated"))
    end
    return nothing
end
