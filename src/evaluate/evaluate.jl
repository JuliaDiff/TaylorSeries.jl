# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#


# Kernels
include("kernels.jl")

# Inplace methods
include("inplace.jl")

## Evaluating ##
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
function evaluate(a::TaylorN{T}, vals::Tuple{S,Vararg{S}};
        sorting::Bool=_defaultsorting(T,S)) where {T<:Number,S<:Number}
    @assert get_numvars(a) == length(vals)
    return _evaluate(a, vals, Val(sorting))
end

function evaluate(a::TaylorN, vals::NTuple{N,<:AbstractSeries};
        sorting::Bool=false) where {N}
    @assert get_numvars(a) == N
    return _evaluate(a, vals, Val(sorting))
end

evaluate(a::TaylorN{T}, vals::AbstractVector{S}; sorting::Bool=_defaultsorting(T,S)) where
    {T<:NumberNotSeries, S<:Number} = evaluate(a, (vals...,); sorting)

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

# High-dimensional array evaluation. Keep this allocating interface as a
# broadcast of scalar evaluations so it retains scalar promotion, sorting,
# shape, and StaticArray container behavior while accepting views as inputs.
# TODO: Preserve concrete result element types for empty arrays once this can
# be done without duplicating the scalar promotion rules.
function evaluate(A::AbstractArray{TaylorN{T}}, vals::AbstractVector{S};
        sorting::Bool=_defaultsorting(T,S)) where {T<:Number,S<:Number}
    return evaluate.(A, Ref(vals); sorting)
end

function evaluate(A::AbstractArray{TaylorN{T}}, vals::Tuple{S,Vararg{S}};
        sorting::Bool=_defaultsorting(T,S)) where {T<:Number,S<:Number}
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
