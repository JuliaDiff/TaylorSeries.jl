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
ommitted, its value is considered as zero. Note that the syntax `a(dx)` is
equivalent to `evaluate(a,dx)`, and `a()` is equivalent to `evaluate(a)`.
"""
function evaluate(a::Taylor1{T}, dx::T) where {T<:Number}
    @inbounds suma = a[end]
    @inbounds for k in a.order-1:-1:0
        suma = suma*dx + a[k]
    end
    suma
end

function evaluate(a::Taylor1{T}, dx::S) where {T<:Number, S<:Number}
    R = promote_type(T, S)
    return evaluate(convert(Taylor1{R}, a), convert(R, dx))
end

evaluate(a::Taylor1{T}) where {T<:Number} = a[0]

"""
    evaluate(x, δt)

Evaluates each element of `x::Union{ Vector{Taylor1{T}}, Matrix{Taylor1{T}} }`,
representing the dependent variables of an ODE, at *time* δt. Note that the
syntax `x(δt)` is equivalent to `evaluate(x, δt)`, and `x()`
is equivalent to `evaluate(x)`.
"""
function evaluate(x::Union{Array{Taylor1{T},1}, SubArray{Taylor1{T},1}}, δt::S) where {T<:Number, S<:Number}
    R = promote_type(T,S)
    return evaluate(convert(Array{Taylor1{R},1},x), convert(R,δt))
end
function evaluate(x::Array{Taylor1{T},1}, δt::T) where {T<:Number}
    xnew = Array{T}(undef, length(x) )
    evaluate!(x, δt, xnew)
    return xnew
end

evaluate(a::Array{Taylor1{T},1}) where {T<:Number} = evaluate(a, zero(T))
evaluate(a::SubArray{Taylor1{T},1}) where {T<:Number} = evaluate(a, zero(T))

function evaluate(A::Union{Array{Taylor1{T},2}, SubArray{Taylor1{T},2}}, δt::S) where {T<:Number, S<:Number}
    R = promote_type(T,S)
    return evaluate(convert(Array{Taylor1{R},2},A), convert(R,δt))
end
function evaluate(A::Array{Taylor1{T},2}, δt::T) where {T<:Number}
    n,m = size(A)
    Anew = Array{T}(undef, n, m )
    xnew = Array{T}(undef, n )

    for i in 1:m
        evaluate!(A[:,i], δt, xnew)
        Anew[:,i] = xnew
    end

    return Anew
end
evaluate(A::Array{Taylor1{T},2}) where {T<:Number} = evaluate.(A)
evaluate(A::SubArray{Taylor1{T},2}) where {T<:Number} = evaluate.(A)

"""
    evaluate!(x, δt, x0)

Evaluates each element of `x::Array{Taylor1{T},1}`,
representing the Taylor expansion for the dependent variables
of an ODE at *time* `δt`. It updates the vector `x0` with the
computed values.
"""
function evaluate!(x::Array{Taylor1{T},1}, δt::T,
        x0::Union{Array{T,1},SubArray{T,1}}) where {T<:Number}

    # @assert length(x) == length(x0)
    @inbounds for i in eachindex(x, x0)
        x0[i] = evaluate( x[i], δt )
    end
    nothing
end
function evaluate!(x::Array{Taylor1{T},1}, δt::S,
        x0::Union{Array{T,1},SubArray{T,1}}) where {T<:Number, S<:Number}

    # @assert length(x) == length(x0)
    @inbounds for i in eachindex(x, x0)
        x0[i] = evaluate( x[i], δt )
    end
    nothing
end

"""
    evaluate(a, x)

Substitute `x::Taylor1` as independent variable in a `a::Taylor1` polynomial.
Note that the syntax `a(x)` is equivalent to `evaluate(a, x)`.
"""
evaluate(a::Taylor1{T}, x::Taylor1{S}) where {T<:Number, S<:Number} =
    evaluate(promote(a,x)...)

function evaluate(a::Taylor1{T}, x::Taylor1{T}) where {T<:Number}
    if a.order != x.order
        a, x = fixorder(a, x)
    end
    @inbounds suma = a[end]
    @inbounds for k = a.order-1:-1:0
        suma = suma*x + a[k]
    end
    suma
end

function evaluate(a::Taylor1{Taylor1{T}}, x::Taylor1{T}) where {T<:Number}
    @inbounds suma = a[end]
    @inbounds for k = a.order-1:-1:0
        suma = suma*x + a[k]
    end
    suma
end

evaluate(p::Taylor1{T}, x::Array{S}) where {T<:Number, S<:Number} =
    evaluate.(p, x)

#function-like behavior for Taylor1
(p::Taylor1)(x) = evaluate(p, x)

(p::Taylor1)() = evaluate(p)

#function-like behavior for Vector{Taylor1}
(p::Array{Taylor1{T},1})(x) where {T<:Number} = evaluate(p, x)
(p::SubArray{Taylor1{T},1})(x) where {T<:Number} = evaluate(p, x)
(p::Array{Taylor1{T},1})() where {T<:Number} = evaluate.(p)
(p::SubArray{Taylor1{T},1})() where {T<:Number} = evaluate.(p)

#function-like behavior for Matrix{Taylor1}
(p::Array{Taylor1{T},2})(x) where {T<:Number} = evaluate(p, x)
(p::SubArray{Taylor1{T},2})(x) where {T<:Number} = evaluate(p, x)
(p::Array{Taylor1{T},2})() where {T<:Number} = evaluate.(p)
(p::SubArray{Taylor1{T},2})() where {T<:Number} = evaluate.(p)

## Evaluation of multivariable
function evaluate!(x::Array{TaylorN{T},1}, δx::Array{T,1},
        x0::Array{T,1}) where {T<:Number}

    # @assert length(x) == length(x0)
    @inbounds for i in eachindex(x, x0)
        x0[i] = evaluate( x[i], δx )
    end
    nothing
end

function evaluate!(x::Array{TaylorN{T},1}, δx::Array{Taylor1{T},1},
        x0::Array{Taylor1{T},1}) where {T<:NumberNotSeriesN}

    # @assert length(x) == length(x0)
    @inbounds for i in eachindex(x, x0)
        x0[i] = evaluate( x[i], δx )
    end
    nothing
end

function evaluate!(x::Array{TaylorN{T},1}, δx::Array{TaylorN{T},1},
        x0::Array{TaylorN{T},1}) where {T<:NumberNotSeriesN}

    # @assert length(x) == length(x0)
    @inbounds for i in eachindex(x, x0)
        x0[i] = evaluate( x[i], δx )
    end
    nothing
end

function evaluate!(x::Array{Taylor1{TaylorN{T}},1}, δt::T,
        x0::Array{TaylorN{T},1}) where {T<:Number}

    # @assert length(x) == length(x0)
    @inbounds for i in eachindex(x, x0)
        x0[i] = evaluate( x[i], δt )
    end
    nothing
end

"""
    evaluate(a, [vals])

Evaluate a `HomogeneousPolynomial` polynomial at `vals`. If `vals` is ommitted,
it's evaluated at zero. Note that the syntax `a(vals)` is equivalent to
`evaluate(a, vals)`; and `a()` is equivalent to `evaluate(a)`.
"""
function evaluate(a::HomogeneousPolynomial{T}, vals::NTuple{N,S} ) where
        {T<:Number, S<:NumberNotSeriesN, N}

    @assert N == get_numvars()

    ct = coeff_table[a.order+1]
    R = promote_type(T,S)
    suma = zero(R)

    for (i,a_coeff) in enumerate(a.coeffs)
        iszero(a_coeff) && continue
        tmp = vals[1]^(ct[i][1])
        for n in 2:N
            tmp *= vals[n]^(ct[i][n])
        end
        suma += a_coeff * tmp
    end

    return suma
end

evaluate(a::HomogeneousPolynomial{T}, vals::Array{S,1} ) where
        {T<:Number, S<:NumberNotSeriesN} = evaluate(a, (vals...,))

evaluate(a::HomogeneousPolynomial, v, vals...) = evaluate(a, (v, vals...,))

evaluate(a::HomogeneousPolynomial, v) = evaluate(a, v...)

function evaluate(a::HomogeneousPolynomial)
    a.order == 0 && return a[1]
    zero(a[1])
end

#function-like behavior for HomogeneousPolynomial
(p::HomogeneousPolynomial)(x) = evaluate(p, x)

(p::HomogeneousPolynomial)(x, v...) = evaluate(p, (x, v...,))

(p::HomogeneousPolynomial)() = evaluate(p)

"""
    evaluate(a, [vals])

Evaluate the `TaylorN` polynomial `a` at `vals`.
If `vals` is ommitted, it's evaluated at zero.
Note that the syntax `a(vals)` is equivalent to `evaluate(a, vals)`; and `a()`
is equivalent to `evaluate(a)`.
"""
function evaluate(a::TaylorN{T}, vals::NTuple{N,S}) where
        {T<:Number,S<:NumberNotSeries, N}

    @assert N == get_numvars()

    R = promote_type(T,S)
    a_length = length(a)
    suma = zeros(R, a_length)
    for homPol in 1:length(a)
        sun = zero(R)
        for (i, a_coeff) in enumerate(a.coeffs[homPol].coeffs)
            iszero(a_coeff) && continue
            tmp = vals[1]^(coeff_table[homPol][i][1])
            for n in 2:N
                tmp *= vals[n]^(coeff_table[homPol][i][n])
            end
            sun += a_coeff * tmp
        end
        suma[homPol] = sun
    end

    return sum( sort!(suma, by=abs2) )
end

evaluate(a::TaylorN, vals) = evaluate(a, (vals...,))

evaluate(a::TaylorN, v, vals...) = evaluate(a, (v, vals...,))

function evaluate(a::TaylorN{T}, vals::NTuple{N,Taylor1{S}}) where
        {T<:Number, S<:NumberNotSeries, N}

    @assert N == get_numvars()

    R = promote_type(T,S)
    a_length = length(a)
    ord = maximum( get_order.(vals) )
    suma = Taylor1(zeros(R, ord))

    for homPol in 1:length(a)
        for (i, a_coeff) in enumerate(a.coeffs[homPol].coeffs)
            iszero(a_coeff) && continue
            tmp = vals[1]^(coeff_table[homPol][i][1])
            for n in 2:N
                tmp *= vals[n]^(coeff_table[homPol][i][n])
            end
            suma += a_coeff * tmp
        end
    end

    return suma
end

evaluate(a::TaylorN{T}, vals::Array{Taylor1{S},1}) where
    {T<:Number, S<:NumberNotSeriesN} = evaluate(a, (vals...,))

function evaluate(a::TaylorN{Taylor1{T}}, vals::NTuple{N, Taylor1{T}}) where
        {T<:NumberNotSeries, N}

    @assert N == get_numvars()

    a_length = length(a)
    ord = maximum( get_order.(vals) )
    suma = Taylor1(zeros(T, ord))

    for homPol in 1:length(a)
        for (i, a_coeff) in enumerate(a.coeffs[homPol].coeffs)
            iszero(a_coeff) && continue
            tmp = vals[1]^(coeff_table[homPol][i][1])
            for n in 2:N
                tmp *= vals[n]^(coeff_table[homPol][i][n])
            end
            suma += a_coeff * tmp
        end
    end

    return suma
end

evaluate(a::TaylorN{Taylor1{T}}, vals::Array{Taylor1{T},1}) where
    {T<:NumberNotSeries} = evaluate(a, (vals...,))

function evaluate(a::TaylorN{T}, vals::NTuple{N, TaylorN{S}}) where
        {T<:Number, S<:NumberNotSeries, N}

    @assert length(vals) == get_numvars()

    num_vars = get_numvars()
    R = promote_type(T,eltype(S))
    a_length = length(a)
    ord = maximum( get_order.(vals) )
    suma = zero(TaylorN{R})

    for homPol in 1:length(a)
        for (i, a_coeff) in enumerate(a.coeffs[homPol].coeffs)
            iszero(a_coeff) && continue
            tmp = vals[1]^(coeff_table[homPol][i][1])
            for n in 2:num_vars
                tmp *= vals[n]^(coeff_table[homPol][i][n])
            end
            suma += a_coeff * tmp
        end
    end

    return suma
end

evaluate(a::TaylorN{T}, vals::Array{TaylorN{S},1}) where
    {T<:Number, S<:NumberNotSeries} = evaluate(a, (vals...,))

function evaluate(a::TaylorN{T}, s::Symbol, val::S) where
        {T<:Number, S<:NumberNotSeriesN}
    vars = get_variables(T)
    ind = lookupvar(s)
    vars[ind] = val
    evaluate(a, vars)
end

evaluate(a::TaylorN{T}, x::Pair{Symbol,S}) where {T<:Number, S<:NumberNotSeriesN} =
    evaluate(a, first(x), last(x))

evaluate(a::TaylorN{T}) where {T<:Number} = a[0][1]

#Vector evaluation
function evaluate(x::Union{Array{TaylorN{T},1},SubArray{TaylorN{T},1}}, δx::Vector{S}) where {T<:Number, S<:Number}
    R = promote_type(T,S)
    return evaluate(convert(Array{TaylorN{R},1},x), convert(Vector{R},δx))
end

function evaluate(x::Array{TaylorN{T},1}, δx::Array{T,1}) where {T<:Number}
    x0 = Array{T}(undef, length(x) )
    evaluate!( x, δx, x0 )
    return x0
end

evaluate(x::Array{TaylorN{T},1}) where {T<:Number} = evaluate.(x)
evaluate(x::SubArray{TaylorN{T},1}) where {T<:Number} = evaluate.(x)

#Matrix evaluation
function evaluate(A::Union{Array{TaylorN{T},2}, SubArray{TaylorN{T},2}}, δx::Vector{S}) where {T<:Number, S<:Number}
    R = promote_type(T,S)
    return evaluate(convert(Array{TaylorN{R},2},A), convert(Vector{R},δx))
end
function evaluate(A::Array{TaylorN{T},2}, δx::Vector{T}) where {T<:Number}
    n,m = size(A)
    Anew = Array{T}(undef, n, m )
    xnew = Array{T}(undef, n )

    for i in 1:m
        evaluate!(A[:,i], δx, xnew)
        Anew[:,i] = xnew
    end

    return Anew
end
evaluate(A::Array{TaylorN{T},2}) where {T<:Number} = evaluate.(A)
evaluate(A::SubArray{TaylorN{T},2}) where {T<:Number} = evaluate.(A)

#function-like behavior for TaylorN
(p::TaylorN)(x) = evaluate(p, x)
(p::TaylorN)() = evaluate(p)
(p::TaylorN)(s::Symbol, x) = evaluate(p, s, x)
(p::TaylorN)(x::Pair) = evaluate(p, first(x), last(x))
(p::TaylorN)(x, v...) = evaluate(p, (x, v...,))

#function-like behavior for Vector{TaylorN}
(p::Array{TaylorN{T},1})(x) where {T<:Number} = evaluate(p, x)
(p::SubArray{TaylorN{T},1})(x) where {T<:Number} = evaluate(p, x)
(p::Array{TaylorN{T},1})() where {T<:Number} = evaluate(p)
(p::SubArray{TaylorN{T},1})() where {T<:Number} = evaluate(p)

#function-like behavior for Matrix{TaylorN}
(p::Array{TaylorN{T},2})(x) where {T<:Number} = evaluate(p, x)
(p::SubArray{TaylorN{T},2})(x) where {T<:Number} = evaluate(p, x)
(p::Array{TaylorN{T},2})() where {T<:Number} = evaluate.(p)
(p::SubArray{TaylorN{T},2})() where {T<:Number} = evaluate.(p)
