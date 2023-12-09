# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

## Evaluating ##
"""
    evaluate(a, [dx])

Evaluate a `Taylor1` polynomial using Horner's rule (hand coded). If `dx` is
omitted, its value is considered as zero. Note that the syntax `a(dx)` is
equivalent to `evaluate(a,dx)`, and `a()` is equivalent to `evaluate(a)`.
"""
function evaluate(a::Taylor1{T}, dx::T) where {T<:Number}
    @inbounds suma = zero(a[end])
    @inbounds for k in reverse(eachindex(a))
        suma = suma*dx + a[k]
    end
    return suma
end
function evaluate(a::Taylor1{T}, dx::S) where {T<:Number, S<:Number}
    suma = a[end]*zero(dx)
    @inbounds for k in reverse(eachindex(a))
        suma = suma*dx + a[k]
    end
    return suma
end
evaluate(a::Taylor1{T}) where {T<:Number} = a[0]


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

function evaluate(a::Taylor1{T}, x::Taylor1{T}) where {T<:Number}
    if a.order != x.order
        a, x = fixorder(a, x)
    end
    @inbounds suma = a[end]*zero(x)
    @inbounds for k in reverse(eachindex(a))
        suma = suma*x + a[k]
    end
    return suma
end

function evaluate(a::Taylor1{Taylor1{T}}, x::Taylor1{T}) where {T<:NumberNotSeriesN}
    @inbounds suma = a[end]*zero(x)
    @inbounds for k in reverse(eachindex(a))
        suma = suma*x + a[k]
    end
    return suma
end
function evaluate(a::Taylor1{T}, x::Taylor1{Taylor1{T}}) where {T<:NumberNotSeriesN}
    @inbounds suma = a[end]*zero(x)
    @inbounds for k in reverse(eachindex(a))
        suma = suma*x + a[k]
    end
    return suma
end

evaluate(p::Taylor1{T}, x::AbstractArray{S}) where {T<:Number, S<:Number} =
    evaluate.(Ref(p), x)

function evaluate(a::Taylor1{TaylorN{T}}, dx::NumberNotSeries) where
        {T<:NumberNotSeries}
    @inbounds suma = TaylorN( zero(T)*constant_term(dx), a[0].order )
    @inbounds for k in reverse(eachindex(a))
        for ordQ in eachindex(a[k])
            mul!(suma, suma, dx, ordQ)
            add!(suma, suma, a[k], ordQ)
        end
    end
    return suma
end

#function-like behavior for Taylor1
(p::Taylor1)(x) = evaluate(p, x)
(p::Taylor1)() = evaluate(p)

#function-like behavior for Vector{Taylor1} (asumes Julia version >= 1.6)
(p::Array{Taylor1{T}})(x) where {T<:Number} = evaluate.(p, x)
(p::SubArray{Taylor1{T}})(x) where {T<:Number} = evaluate.(p, x)
(p::Array{Taylor1{T}})() where {T<:Number} = evaluate.(p)
(p::SubArray{Taylor1{T}})() where {T<:Number} = evaluate.(p)


"""
    evaluate(a, [vals])

Evaluate a `HomogeneousPolynomial` polynomial at `vals`. If `vals` is omitted,
it's evaluated at zero. Note that the syntax `a(vals)` is equivalent to
`evaluate(a, vals)`; and `a()` is equivalent to `evaluate(a)`.
"""
function evaluate(a::HomogeneousPolynomial, vals::NTuple)
    @assert length(vals) == get_numvars()

    return _evaluate(a, vals)
end

function _evaluate(a::HomogeneousPolynomial, vals::NTuple)
    ct = coeff_table[a.order+1]
    suma = zero(a[1])*vals[1]

    for (i, a_coeff) in enumerate(a.coeffs)
        iszero(a_coeff) && continue
        @inbounds tmp = prod( vals .^ ct[i] )
        suma += a_coeff * tmp
    end

    return suma
end


evaluate(a::HomogeneousPolynomial{T}, vals::AbstractArray{S,1} ) where
        {T<:Number,S<:NumberNotSeriesN} = _evaluate(a, (vals...,))

evaluate(a::HomogeneousPolynomial, v, vals::Vararg) = evaluate(a, (v, vals...,))

evaluate(a::HomogeneousPolynomial, v) = evaluate(a, [v...])

function evaluate(a::HomogeneousPolynomial)
    a.order == 0 && return a[1]
    return zero(a[1])
end

#function-like behavior for HomogeneousPolynomial
(p::HomogeneousPolynomial)(x) = evaluate(p, x)

(p::HomogeneousPolynomial)(x, v::Vararg{T, N}) where {T,N} = evaluate(p, (x, v...,))

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
`evaluate(a)`. No extension exists that incorporates
`sorting`.
"""
evaluate(a::TaylorN, vals; sorting::Bool=true) = _evaluate(a, (vals...,), Val(sorting))

evaluate(a::TaylorN, vals::NTuple; sorting::Bool=true) = _evaluate(a, vals, Val(sorting))

evaluate(a::TaylorN, v, vals::Vararg; sorting::Bool=true) =
    _evaluate(a, (v, vals...,), Val(sorting))

evaluate(a::TaylorN, vals::NTuple{N,<:AbstractSeries}; sorting::Bool=false) where
    {N} = _evaluate(a, vals, Val(sorting))

evaluate(a::TaylorN{Taylor1}, vals::NTuple{N,<:AbstractSeries};
    sorting::Bool=false) where {N} = _evaluate(a, vals, Val(sorting))

evaluate(a::TaylorN{Taylor1}, vals::NTuple; sorting::Bool=false) =
    _evaluate(a, vals, Val(sorting))

evaluate(a::TaylorN{Taylor1}, v::Number, vals::Vararg; sorting::Bool=false) =
    _evaluate(a, (v, vals...,), Val(sorting))

evaluate(a::TaylorN, vals::AbstractVector{T}) where {S, T<:AbstractSeries{S}} =
    _evaluate(a, (vals...,), Val(false))

evaluate(a::TaylorN{Taylor1}, vals::AbstractVector{T}) where {T} =
    _evaluate(a, (vals...,), Val(false))

function evaluate(a::TaylorN{T}, s::Symbol, val::S) where {T<:Number, S<:NumberNotSeriesN}
    vars = get_variables(T)
    ind = lookupvar(s)
    vars[ind] = val + zero(vars[ind])
    return evaluate(a, vars)
end

evaluate(a::TaylorN{T}, x::Pair{Symbol,S}) where {T<:NumberNotSeries,S} =
    evaluate(a, first(x), last(x))

evaluate(a::TaylorN{T}, x::Pair{Symbol,S}) where {T<:AbstractSeries,S} =
    evaluate(a, first(x), last(x))

evaluate(a::TaylorN{T}) where {T<:Number} = a[0][1]


# _evaluate
# Returns a vector with the evaluation of the HomogeneousPolynomials
function _evaluate(a::TaylorN, vals)
    @assert get_numvars() == length(vals)
    R = promote_type(numtype(a), typeof.(vals)...)
    a_length = length(a)
    suma = zeros(R, a_length)
    @inbounds for homPol in eachindex(a)
        suma[homPol+1] = evaluate(a[homPol], vals)
    end
    return suma
end

function _evaluate(a::TaylorN{T}, vals::NTuple, ::Val{true}) where {T<:NumberNotSeries}
    suma = _evaluate(a, vals)
    return sum( sort!(suma, by=abs2) )
end
function _evaluate(a::TaylorN{T}, vals::NTuple, ::Val{false}) where {T<:Number}
    suma = _evaluate(a, vals)
    return sum( suma )
end


#High-dim array evaluation
function evaluate(A::AbstractArray{TaylorN{T},N}, δx::Vector{S}) where {T<:Number,S<:Number,N}
    R = promote_type(T,S)
    return evaluate(convert(Array{TaylorN{R},N},A), convert(Vector{R},δx))
end
function evaluate(A::Array{TaylorN{T}}, δx::Vector{T}) where {T<:Number}
    Anew = Array{T}(undef, size(A)...)
    evaluate!(A, δx, Anew)
    return Anew
end
evaluate(A::AbstractArray{TaylorN{T}}) where {T<:Number} = evaluate.(A)

#function-like behavior for TaylorN
(p::TaylorN)(x) = evaluate(p, x)
(p::TaylorN)() = evaluate(p)
(p::TaylorN)(s::Symbol, x) = evaluate(p, s, x)
(p::TaylorN)(x::Pair) = evaluate(p, first(x), last(x))
(p::TaylorN)(x, v::Vararg{T, N}) where {T,N} = evaluate(p, (x, v...,))
(p::TaylorN)(b::Bool, x) = evaluate(p, x, sorting=b)
(p::TaylorN)(b::Bool, x, v::Vararg{T, N}) where {T,N} = evaluate(p, (x, v...,), sorting=b)

#function-like behavior for AbstractArray{TaylorN{T}}
(p::Array{TaylorN{T}})(x) where {T<:Number} = evaluate(p, x)
(p::SubArray{TaylorN{T}})(x) where {T<:Number} = evaluate(p, x)
(p::Array{TaylorN{T}})() where {T<:Number} = evaluate(p)
(p::SubArray{TaylorN{T}})() where {T<:Number} = evaluate(p)


"""
    evaluate!(x, δt, x0)

Evaluates each element of `x::AbstractArray{Taylor1{T}}`,
representing the Taylor expansion for the dependent variables
of an ODE at *time* `δt`. It updates the vector `x0` with the
computed values.
"""
function evaluate!(x::AbstractArray{Taylor1{T}}, δt::S,
        x0::AbstractArray{T}) where {T<:Number, S<:Number}
    x0 .= evaluate.( x, δt )
    return nothing
end


## In place evaluation of multivariable arrays
function evaluate!(x::AbstractArray{TaylorN{T}}, δx::Array{T,1},
        x0::AbstractArray{T}) where {T<:Number}
    @inbounds for i in eachindex(x, x0)
        x0[i] = evaluate( x[i], δx )
    end
    return nothing
end

function evaluate!(x::AbstractArray{TaylorN{T}}, δx::Array{Taylor1{T},1},
        x0::AbstractArray{Taylor1{T}}) where {T<:NumberNotSeriesN}
    @inbounds for i in eachindex(x, x0)
        x0[i] = evaluate( x[i], δx )
    end
    return nothing
end

function evaluate!(x::AbstractArray{TaylorN{T}}, δx::Array{TaylorN{T},1},
        x0::AbstractArray{TaylorN{T}}; sorting::Bool=true) where {T<:NumberNotSeriesN}
    @inbounds for i in eachindex(x, x0)
        x0[i] = _evaluate( x[i], δx, Val(sorting) )
    end
    return nothing
end

function evaluate!(x::AbstractArray{TaylorN{T}}, δt::T,
        x0::AbstractArray{TaylorN{T}}) where {T<:Number}
    @inbounds for i in eachindex(x, x0)
        x0[i] = evaluate( x[i], δt )
    end
    return nothing
end
