# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#


"""
    _defaultsorting(T, S)

Function to set defaults for sorting, useful for evaluating with `Interval`s
"""
_defaultsorting(::Type{T}, ::Type{S}) where {T,S} =
    !(T <: AbstractSeries || S <: AbstractSeries)


# evaluate_kernels.jl: the internal machinery. That's _evaluate, _evaluate!, _horner!, _defaultsorting and _valtype. Both public files call into these, and the extension overrides several of them.
"""
    _evaluate(a::Taylor1, dx::NumberNotSeries)

Evaluate `a` at a number `dx` using Horner's rule. Neither the coefficients
nor `dx` are Taylor series. Both `evaluate` and the array form of `evaluate!`
use this method; it does not create temporary arrays.
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
    _evaluate(a::HomogeneousPolynomial, vals)

Evaluate one homogeneous polynomial at `vals`, computing each term directly
without first collecting powers or terms in arrays. Supply one value per
variable. If the values are `TaylorN` series, their `JetSpace`s must match
the polynomial's.
"""
function _evaluate(a::HomogeneousPolynomial{T},
        vals) where {T}
    order(a) == 0 && return a[1]*one(vals[1])
    ct = a.space.coeff_table[order(a)+1]
    suma = zero(a[1]*vals[1])
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

# For ordinary numbers, the coefficient and value types determine the result type.
# For series, the method above uses the values themselves to preserve their JetSpace.
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


# TODO: avoid allocating a new array every time we evaluate with `sorting=true`.
# This could be achieved by passing a reusable buffer for the evaluated contributions,
# as an additional input argument. This buffer would be filled and sorted on each call,
# then summed.
"""
    _evaluate(a::TaylorN, vals::NTuple, ::Val{true})
    _evaluate(a::TaylorN, vals::NTuple, ::Val{false})

Evaluate `a` at `vals` and return the sum of its contributions.
`Val(true)` sorts the contributions from smallest to largest magnitude before
adding them, to reduce rounding error. `Val(false)` skips this sorting step.

The third argument is needed because `_evaluate(a::TaylorN, vals::NTuple)`
returns an array with one contribution per polynomial degree, rather than
their sum. The sorted method needs those separate values so it can reorder
them before adding them. Using `Val(true)` or `Val(false)` selects the method
for the requested summation without changing the two-argument method's result.

The sorted method currently allocates an array to hold the contributions
before sorting and adding them.
"""
_evaluate(a::TaylorN{T}, vals::NTuple, ::Val{true}) where
    {T<:NumberNotSeries} = sum( sort!(_evaluate(a, vals), by=abs2) )

_evaluate(a::TaylorN{T}, vals::NTuple, ::Val{false}) where {T<:Number} =
    sum( _evaluate(a, vals) )

function _evaluate(a::TaylorN{T}, vals::NTuple{N,<:TaylorN}, ::Val{false}) where
        {N,T<:Number}
    R = promote_type(T, TS.numtype(vals[1]))
    _check_same_space(a, vals[1])
    a = convert(TaylorN{R}, a)
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
    _evaluate(a::TaylorN, vals::NTuple)

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
`_evaluate(a::TaylorN, vals::NTuple)`, which returns the separate contributions.
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


# _evaluate!
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

function _evaluate!(a::TaylorN{T}, vals::NTuple{N,TaylorN{T}}, dest::TaylorN{T},
        valscache::Vector{TaylorN{T}}, aux::TaylorN{T}) where {N,T<:Number}
    @inbounds for homPol in eachindex(a)
        _evaluate!(dest, a[homPol], vals, valscache, aux)
    end
    return nothing
end

# The public method checks the inputs before either method changes dest.
function _evaluate!(a::TaylorN{T}, vals::NTuple{N,TaylorN{T}},
        dest::TaylorN{T}, valscache::Vector{TaylorN{T}},
        aux::TaylorN{T}, ::Val{false}) where {N,T<:Number}
    zero!(dest)
    _evaluate!(a, vals, dest, valscache, aux)
    return nothing
end

function _evaluate!(a::TaylorN{T}, vals::NTuple{N,TaylorN{T}},
        dest::TaylorN{T}, valscache::Vector{TaylorN{T}},
        aux::TaylorN{T}, ::Val{true}) where {N,T<:Number}
    result = evaluate(a, vals; sorting=true)
    zero!(dest)
    for ord in eachindex(dest)
        identity!(dest, result, ord)
    end
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
