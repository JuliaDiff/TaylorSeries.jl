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
        suma = suma * dx + a[k]
    end
    return suma
end

function evaluate(a::Taylor1{T}, dx::S) where {T<:Number, S<:Number}
    suma = a[end]*zero(dx)
    @inbounds for k in reverse(eachindex(a))
        suma = suma * dx + a[k]
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
    aux = zero(suma)
    _horner!(suma, a, x, aux)
    return suma
end

function evaluate(a::Taylor1{Taylor1{T}}, x::Taylor1{T}) where {T<:NumberNotSeriesN}
    @inbounds suma = a[end]*zero(x)
    aux = zero(suma)
    _horner!(suma, a, x, aux)
    return suma
end

function evaluate(a::Taylor1{T}, x::Taylor1{Taylor1{T}}) where {T<:NumberNotSeriesN}
    @inbounds suma = a[end]*zero(x)
    aux = zero(suma)
    _horner!(suma, a, x, aux)
    return suma
end

evaluate(p::Taylor1{T}, x::AbstractArray{S}) where {T<:Number, S<:Number} =
    evaluate.(Ref(p), x)

# Substitute a TaylorN into a Taylor1
function evaluate(a::Taylor1{T}, dx::TaylorN{T}) where {T<:NumberNotSeries}
    suma = TaylorN( zero(T), dx.order)
    aux = TaylorN( zero(T), dx.order)
    _horner!(suma, a, dx, aux)
    return suma
end

function evaluate(a::Taylor1{T}, dx::Taylor1{TaylorN{T}}) where {T<:NumberNotSeries}
    if a.order != dx.order
        a, dx = fixorder(a, dx)
    end
    suma = Taylor1( zero(dx[0]), a.order)
    aux  = Taylor1( zero(dx[0]), a.order)
    _horner!(suma, a, dx, aux)
    return suma
end


# Evaluate a Taylor1{TaylorN{T}} on Vector{TaylorN} is interpreted
# as a substitution on the TaylorN vars
function evaluate(a::Taylor1{TaylorN{T}}, dx::Vector{TaylorN{T}}) where {T<:NumberNotSeries}
    @assert length(dx) == get_numvars()
    suma = Taylor1( zero(a[0]), a.order)
    suma.coeffs .= evaluate.(a[:], Ref(dx))
    return suma
end

function evaluate(a::Taylor1{TaylorN{T}}, ind::Int, dx::T) where {T<:NumberNotSeries}
    @assert (1 ≤ ind ≤ get_numvars()) "Invalid `ind`; it must be between 1 and `get_numvars()`"
    suma = Taylor1( zero(a[0]), a.order)
    for ord in eachindex(suma)
        for ordQ in eachindex(a[0])
            _evaluate!(suma[ord], a[ord][ordQ], ind, dx)
        end
    end
    return suma
end

function evaluate(a::Taylor1{TaylorN{T}}, ind::Int, dx::TaylorN{T}) where {T<:NumberNotSeries}
    @assert (1 ≤ ind ≤ get_numvars()) "Invalid `ind`; it must be between 1 and `get_numvars()`"
    suma = Taylor1( zero(a[0]), a.order)
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
function evaluate(a::HomogeneousPolynomial, vals::NTuple{N,<:Number}) where {N}
    @assert length(vals) == get_numvars()
    return _evaluate(a, vals)
end

evaluate(a::HomogeneousPolynomial{T}, vals::AbstractArray{S,1} ) where
    {T<:Number,S<:NumberNotSeriesN} = evaluate(a, (vals...,))

evaluate(a::HomogeneousPolynomial, v, vals::Vararg{Number,N}) where {N} =
    evaluate(a, promote(v, vals...,))

evaluate(a::HomogeneousPolynomial, v) = evaluate(a, promote(v...,))

function evaluate(a::HomogeneousPolynomial{T}) where {T}
    a.order == 0 && return a[1]
    return zero(a[1])
end

# Internal method that avoids checking that the length of `vals` is the appropriate
function _evaluate(a::HomogeneousPolynomial{T}, vals::NTuple) where {T}
    # @assert length(vals) == get_numvars()
    a.order == 0 && return a[1]*one(vals[1])
    ct = coeff_table[a.order+1]
    suma = zero(a[1])*vals[1]
    vv = vals .^ ct[1]
    for (i, a_coeff) in enumerate(a.coeffs)
        iszero(a_coeff) && continue
        @inbounds vv .= vals .^ ct[i]
        tmp = prod( vv )
        suma += a_coeff * tmp
    end
    return suma
end

function _evaluate!(res::TaylorN{T}, a::HomogeneousPolynomial{T}, vals::NTuple{N,<:TaylorN{T}},
            valscache::Vector{TaylorN{T}}, aux::TaylorN{T}) where {N,T<:NumberNotSeries}
    ct = coeff_table[a.order+1]
    for el in eachindex(valscache)
        power_by_squaring!(valscache[el], vals[el], aux, ct[1][el])
    end
    for (i, a_coeff) in enumerate(a.coeffs)
        iszero(a_coeff) && continue
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

function _evaluate(a::HomogeneousPolynomial{T}, vals::NTuple{N,<:TaylorN{T}}) where
        {N,T<:NumberNotSeries}
    # @assert length(vals) == get_numvars()
    a.order == 0 && return a[1]*one(vals[1])
    suma = TaylorN(zero(T), vals[1].order)
    valscache = [zero(val) for val in vals]
    aux = zero(suma)
    _evaluate!(suma, a, vals, valscache, aux)
    return suma
end

function _evaluate(a::HomogeneousPolynomial{T}, ind::Int, val::T) where {T<:NumberNotSeries}
    suma = TaylorN(zero(T), get_order())
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
    @assert get_numvars() == N
    return _evaluate(a, vals, Val(sorting))
end

function evaluate(a::TaylorN, vals::NTuple{N,<:AbstractSeries};
        sorting::Bool=false) where {N}
    @assert get_numvars() == N
    return _evaluate(a, vals, Val(sorting))
end

evaluate(a::TaylorN{T}, vals::AbstractVector{<:Number};
        sorting::Bool=true) where {T<:NumberNotSeries} =
    evaluate(a, (vals...,); sorting=sorting)

evaluate(a::TaylorN{T}, vals::AbstractVector{<:AbstractSeries};
        sorting::Bool=false) where {T<:NumberNotSeries} =
    evaluate(a, (vals...,); sorting=sorting)

evaluate(a::TaylorN{Taylor1{T}}, vals::AbstractVector{S};
        sorting::Bool=false) where {T, S} =
    evaluate(a, (vals...,); sorting=sorting)

function evaluate(a::TaylorN{T}, s::Symbol, val::S) where
        {T<:Number, S<:NumberNotSeriesN}
    ind = lookupvar(s)
    @assert (1 ≤ ind ≤ get_numvars()) "Symbol is not a TaylorN variable; see `get_variable_names()`"
    return evaluate(a, ind, val)
end

function evaluate(a::TaylorN{T}, ind::Int, val::S) where
        {T<:Number, S<:NumberNotSeriesN}
    @assert (1 ≤ ind ≤ get_numvars()) "Invalid `ind`; it must be between 1 and `get_numvars()`"
    R = promote_type(T,S)
    return _evaluate(convert(TaylorN{R}, a), ind, convert(R, val))
end

function evaluate(a::TaylorN{T}, s::Symbol, val::TaylorN) where {T<:Number}
    ind = lookupvar(s)
    @assert (1 ≤ ind ≤ get_numvars()) "Symbol is not a TaylorN variable; see `get_variable_names()`"
    return evaluate(a, ind, val)
end

function evaluate(a::TaylorN{T}, ind::Int, val::TaylorN) where {T<:Number}
    @assert (1 ≤ ind ≤ get_numvars()) "Invalid `ind`; it must be between 1 and `get_numvars()`"
    a, val = fixorder(a, val)
    a, val = promote(a, val)
    return _evaluate(a, ind, val)
end

evaluate(a::TaylorN{T}, x::Pair{Symbol,S}) where {T, S} =
    evaluate(a, first(x), last(x))

evaluate(a::TaylorN{T}) where {T<:Number} = constant_term(a)


# _evaluate
_evaluate(a::TaylorN{T}, vals::NTuple, ::Val{true}) where {T<:NumberNotSeries} =
    sum( sort!(_evaluate(a, vals), by=abs2) )

_evaluate(a::TaylorN{T}, vals::NTuple, ::Val{false}) where {T<:Number} =
    sum( _evaluate(a, vals) )

function _evaluate(a::TaylorN{T}, vals::NTuple{N,<:TaylorN}, ::Val{false}) where {N,T<:Number}
    R = promote_type(T, TS.numtype(vals[1]))
    res = TaylorN(zero(R), vals[1].order)
    valscache = [zero(val) for val in vals]
    aux = zero(res)
    @inbounds for homPol in eachindex(a)
        _evaluate!(res, a[homPol], vals, valscache, aux)
    end
    return res
end

function _evaluate(a::TaylorN{T}, vals::NTuple{N,<:Number}) where {N,T<:Number}
    R = promote_type(T, typeof(vals[1]))
    suma = zeros(R, length(a))
    @inbounds for homPol in eachindex(a)
        suma[homPol+1] = _evaluate(a[homPol], vals)
    end
    return suma
end

function _evaluate!(res::Vector{TaylorN{T}}, a::TaylorN{T}, vals::NTuple{N,<:TaylorN},
        valscache::Vector{TaylorN{T}}, aux::TaylorN{T}) where {N,T<:Number}
    @inbounds for homPol in eachindex(a)
        _evaluate!(res[homPol+1], a[homPol], vals, valscache, aux)
    end
    return nothing
end

function _evaluate(a::TaylorN{T}, vals::NTuple{N,<:TaylorN}) where {N,T<:Number}
    R = promote_type(T, TS.numtype(vals[1]))
    suma = [TaylorN(zero(R), vals[1].order) for _ in eachindex(a)]
    valscache = [zero(val) for val in vals]
    aux = zero(suma[1])
    _evaluate!(suma, a, vals, valscache, aux)
    return suma
end


function _evaluate(a::TaylorN{T}, ind::Int, val::T) where {T<:NumberNotSeriesN}
    suma = TaylorN(zero(a[0]*val), a.order)
    vval = convert(numtype(suma), val)
    suma, a = promote(suma, a)
    @inbounds for ordQ in eachindex(a)
        _evaluate!(suma, a[ordQ], ind, vval)
    end
    return suma
end

function _evaluate(a::TaylorN{T}, ind::Int, val::TaylorN{T}) where {T<:NumberNotSeriesN}
    suma = TaylorN(zero(a[0]), a.order)
    aux = zero(suma)
    @inbounds for ordQ in eachindex(a)
        _evaluate!(suma, a[ordQ], ind, val, aux)
    end
    return suma
end

function _evaluate!(suma::TaylorN{T}, a::HomogeneousPolynomial{T}, ind::Int, val::T) where
        {T<:NumberNotSeriesN}
    order = a.order
    if order == 0
        suma[0] = a[1]*one(val)
        return nothing
    end
    vv = val .^ (0:order)
    # ct = @isonethread coeff_table[order+1]
    ct = deepcopy(coeff_table[order+1])
    for (i, a_coeff) in enumerate(a.coeffs)
        iszero(a_coeff) && continue
        if ct[i][ind] == 0
            suma[order][i] += a_coeff
            continue
        end
        vpow = ct[i][ind]
        red_order = order - vpow
        ct[i][ind] -= vpow
        kdic = in_base(get_order(), ct[i])
        ct[i][ind] += vpow
        pos = pos_table[red_order+1][kdic]
        suma[red_order][pos] += a_coeff * vv[vpow+1]
    end
    return nothing
end

function _evaluate!(suma::TaylorN{T}, a::HomogeneousPolynomial{T}, ind::Int,
        val::TaylorN{T}, aux::TaylorN{T}) where {T<:NumberNotSeriesN}
    order = a.order
    if order == 0
        suma[0] = a[1]
        return nothing
    end
    vv = zero(suma)
    ct = coeff_table[order+1]
    za = zero(a)
    for (i, a_coeff) in enumerate(a.coeffs)
        iszero(a_coeff) && continue
        if ct[i][ind] == 0
            suma[order][i] += a_coeff
            continue
        end
        za[i] = a_coeff
        zero!(aux)
        _evaluate!(aux, za, ind, one(T))
        za[i] = zero(T)
        vpow = ct[i][ind]
        # vv = val ^ vpow
        if constant_term(val) == 0
            vv = val ^ vpow
        else
            for ordQ in eachindex(val)
                zero!(vv, ordQ)
                pow!(vv, val, vv, vpow, ordQ)
            end
        end
        for ordQ in eachindex(suma)
            mul!(suma, vv, aux, ordQ)
        end
    end
    return nothing
end


#High-dim array evaluation
function evaluate(A::AbstractArray{TaylorN{T},N}, δx::Vector{S}) where
        {T<:Number, S<:Number, N}
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
(p::TaylorN)(s::S, x) where {S<:Union{Symbol, Int}}= evaluate(p, s, x)
(p::TaylorN)(x::Pair) = evaluate(p, first(x), last(x))
(p::TaylorN)(x, v::Vararg{T}) where {T} = evaluate(p, (x, v...,))
(p::TaylorN)(b::Bool, x) = evaluate(p, x, sorting=b)
(p::TaylorN)(b::Bool, x, v::Vararg{T}) where {T} = evaluate(p, (x, v...,), sorting=b)

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
# function evaluate!(x::AbstractArray{Taylor1{Taylor1{T}}}, δt::Taylor1{T},
#         x0::AbstractArray{Taylor1{T}}) where {T<:Number}
#     x0 .= evaluate.( x, Ref(δt) )
#     # x0 .= evaluate.( x, δt )
#     return nothing
# end

## In place evaluation of multivariable arrays
function evaluate!(x::AbstractArray{TaylorN{T}}, δx::Array{T,1},
        x0::AbstractArray{T}) where {T<:Number}
    x0 .= evaluate.( x, Ref(δx) )
    return nothing
end

function evaluate!(x::AbstractArray{TaylorN{T}}, δx::Array{TaylorN{T},1},
        x0::AbstractArray{TaylorN{T}}; sorting::Bool=true) where {T<:NumberNotSeriesN}
    x0 .= evaluate.( x, Ref(δx), sorting = sorting)
    return nothing
end

function evaluate!(a::TaylorN{T}, vals::NTuple{N,TaylorN{T}}, dest::TaylorN{T},
        valscache::Vector{TaylorN{T}}, aux::TaylorN{T}) where {N,T<:Number}
    @inbounds for homPol in eachindex(a)
        _evaluate!(dest, a[homPol], vals, valscache, aux)
    end
    return nothing
end

function evaluate!(a::AbstractArray{TaylorN{T}}, vals::NTuple{N,TaylorN{T}},
        dest::AbstractArray{TaylorN{T}}) where {N,T<:Number}
    # initialize evaluation cache
    valscache = [zero(val) for val in vals]
    aux = zero(dest[1])
    # loop over elements of `a`
    for i in eachindex(a)
        (!iszero(dest[i])) && zero!(dest[i])
        evaluate!(a[i], vals, dest[i], valscache, aux)
    end
    return nothing
end

# In-place Horner methods, used when the result of an evaluation (substitution)
# is Taylor1{}
function _horner!(suma::Taylor1{T}, a::Taylor1{T}, x::Taylor1{T},
        aux::Taylor1{T}) where {T<:Number}
    @inbounds for k in reverse(eachindex(a))
        for ord in eachindex(aux)
            mul!(aux, suma, x, ord)
        end
        for ord in eachindex(aux)
            identity!(suma, aux, ord)
        end
        add!(suma, suma, a[k], 0)
    end
    return nothing
end

function _horner!(suma::Taylor1{T}, a::Taylor1{Taylor1{T}}, x::Taylor1{T},
        aux::Taylor1{T}) where {T<:Number}
    @inbounds for k in reverse(eachindex(a))
        for ord in eachindex(aux)
            mul!(aux, suma, x, ord)
        end
        for ord in eachindex(aux)
            identity!(suma, aux, ord)
            add!(suma, suma, a[k], ord)
        end
    end
    return nothing
end

function _horner!(suma::Taylor1{Taylor1{T}}, a::Taylor1{T}, x::Taylor1{Taylor1{T}},
        aux::Taylor1{Taylor1{T}}) where {T<:Number}
    @inbounds for k in reverse(eachindex(a))
        for ord in eachindex(aux)
            mul!(aux, suma, x, ord)
        end
        for ord in eachindex(aux)
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
        # aux = suma * dx
        for ord in eachindex(aux)
            zero!(aux, ord)
            mul!(aux, suma, dx, ord)
        end
        for ord in eachindex(aux)
            identity!(suma, aux, ord)
        end
        add!(suma, suma, a[k], 0)
    end
    return suma
end