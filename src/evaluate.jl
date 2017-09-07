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
ommitted, its value is considered as zero.
"""
function evaluate{T<:Number}(a::Taylor1{T}, dx::T)
    @inbounds suma = a[end]
    @inbounds for k = a.order:-1:1
        suma = suma*dx + a[k]
    end
    suma
end

function evaluate{T<:Number,S<:Number}(a::Taylor1{T}, dx::S)
    R = promote_type(T,S)
    @inbounds suma = convert(R, a[end])
    @inbounds for k = a.order:-1:1
        suma = suma*dx + a[k]
    end
    suma
end

evaluate{T<:Number}(a::Taylor1{T}) = a[1]

doc"""
    evaluate(x, δt)

Evaluates each element of `x::Array{Taylor1{T},1}`, representing
the dependent variables of an ODE, at *time* δt.
"""
function evaluate{T<:Number, S<:Number}(x::Array{Taylor1{T},1}, δt::S)
    R = promote_type(T,S)
    return evaluate(convert(Array{Taylor1{R},1},x), convert(R,δt))
end

function evaluate{T<:Number}(x::Array{Taylor1{T},1}, δt::T)
    xnew = Array{T}( length(x) )
    evaluate!(x, δt, xnew)
    return xnew
end
evaluate{T<:Number}(a::Array{Taylor1{T},1}) = evaluate.(a)

doc"""
    evaluate!(x, δt, x0)

Evaluates each element of `x::Array{Taylor1{T},1}`,
representing the Taylor expansion for the dependent variables
of an ODE at *time* `δt`. It updates the vector `x0` with the
computed values.
"""
function evaluate!{T<:Number}(x::Array{Taylor1{T},1}, δt::T, x0::Union{Array{T,1},SubArray{T,1}})
    @assert length(x) == length(x0)
    @inbounds for i in eachindex(x)
        x0[i] = evaluate( x[i], δt )
    end
    nothing
end

function evaluate!{T<:Number, S<:Number}(x::Array{Taylor1{T},1}, δt::S, x0::Union{Array{T,1},SubArray{T,1}})
    @assert length(x) == length(x0)
    @inbounds for i in eachindex(x)
        x0[i] = evaluate( x[i], δt )
    end
    nothing
end

"""
    evaluate(a, x)

Substitute `x::Taylor1` as independent variable in a `a::Taylor1` polynomial.
"""
evaluate{T<:Number,S<:Number}(a::Taylor1{T}, x::Taylor1{S}) =
    evaluate(promote(a,x)...)
function evaluate{T<:Number}(a::Taylor1{T}, x::Taylor1{T})
    if a.order != x.order
        a, x = fixorder(a, x)
    end
    @inbounds suma = a[end]
    @inbounds for k = a.order:-1:1
        suma = suma*x + a[k]
    end
    suma
end
function evaluate{T<:NumberNotSeries}(a::Taylor1{Taylor1{T}}, x::Taylor1{T})
    @inbounds suma = a[end]
    @inbounds for k = a.order:-1:1
        suma = suma*x + a[k]
    end
    suma
end

evaluate{T<:Number,S<:Number}(p::Taylor1{T},x::Array{S}) = evaluate.(p,x)

#function-like behavior for Taylor1
(p::Taylor1)(x) = evaluate(p, x)

(p::Taylor1)() = evaluate(p)

#function-like behavior for Array{Taylor1,1}
(p::Array{Taylor1{T},1}){T<:Number}(x) = evaluate(p, x)

(p::Array{Taylor1{T},1}){T<:Number}() = evaluate.(p)

## Evaluation of multivariable
function evaluate!{T<:Number}(x::Array{TaylorN{T},1}, δx::Array{T,1},
        x0::Array{T,1})
    @assert length(x) == length(x0)
    @inbounds for i in eachindex(x)
        x0[i] = evaluate( x[i], δx )
    end
    nothing
end

function evaluate!{T<:NumberNotSeries}(x::Array{TaylorN{T},1}, δx::Array{Taylor1{T},1},
        x0::Array{Taylor1{T},1})
    @assert length(x) == length(x0)
    @inbounds for i in eachindex(x)
        x0[i] = evaluate( x[i], δx )
    end
    nothing
end

function evaluate!{T<:NumberNotSeries}(x::Array{TaylorN{Taylor1{T}},1}, δx::Array{Taylor1{T},1},
    x0::Array{Taylor1{T},1})
@assert length(x) == length(x0)
@inbounds for i in eachindex(x)
    x0[i] = evaluate( x[i], δx )
end
nothing
end

function evaluate!{T<:Number}(x::Array{Taylor1{TaylorN{T}},1}, δt::T,
        x0::Array{TaylorN{T},1})
    @assert length(x) == length(x0)
    @inbounds for i in eachindex(x)
        x0[i] = evaluate( x[i], δt )
    end
    nothing
end

"""
    evaluate(a, vals)

Evaluate a `HomogeneousPolynomial` polynomial at `vals`.
"""
function evaluate{T<:Number,S<:NumberNotSeriesN}(a::HomogeneousPolynomial{T},
        vals::Array{S,1} )
    @assert length(vals) == get_numvars()

    num_vars = get_numvars()
    ct = coeff_table[a.order+1]
    R = promote_type(T,S)
    suma = zero(R)

    for (i,a_coeffs) in enumerate(a.coeffs)
        tmp = vals[1]^(ct[i][1])
        for n in 2:num_vars
            tmp *= vals[n]^(ct[i][n])
        end
        suma += a_coeffs * tmp
    end

    return suma
end

evaluate(a::HomogeneousPolynomial) = zero(a[1])

#function-like behavior for HomogeneousPolynomial
(p::HomogeneousPolynomial)(x) = evaluate(p, x)

(p::HomogeneousPolynomial)() = evaluate(p)

"""
    evaluate(a, [vals])

Evaluate the `TaylorN` polynomial `a` at `vals`.
If `vals` is ommitted, it's evaluated at zero.
"""
function evaluate{T<:Number,S<:NumberNotSeries}(a::TaylorN{T},
        vals::Array{S,1})
    @assert length(vals) == get_numvars()

    num_vars = get_numvars()
    ct = coeff_table
    R = promote_type(T,S)
    a_length = length(a)
    suma = zeros(R,a_length)
    for homPol in 1:length(a)
        sun = zero(R)
        for (i,a_coeff) in enumerate(a.coeffs[homPol].coeffs)
            tmp = vals[1]^(ct[homPol][i][1])
            for n in 2:num_vars
                tmp *= vals[n]^(ct[homPol][i][n])
            end
            sun += a_coeff * tmp
        end
        suma[homPol] = sun
    end

    return sum( sort!(suma, by=abs2) )
end

function evaluate{T<:Number,S<:NumberNotSeries}(a::TaylorN{T},
        vals::Array{Taylor1{S},1})
    @assert length(vals) == get_numvars()

    num_vars = get_numvars()
    ct = coeff_table
    R = promote_type(T,S)
    a_length = length(a)
    ord = maximum( get_order.(vals) )
    suma = Taylor1(zeros(R, ord))

    for homPol in 1:length(a)
        sun = zero(R)
        for (i,a_coeff) in enumerate(a.coeffs[homPol].coeffs)
            tmp = vals[1]^(ct[homPol][i][1])
            for n in 2:num_vars
                tmp *= vals[n]^(ct[homPol][i][n])
            end
            suma += a_coeff * tmp
        end
    end

    return suma
end

function evaluate{T<:NumberNotSeries}(a::TaylorN{Taylor1{T}},
        vals::Array{Taylor1{T},1})
    @assert length(vals) == get_numvars()

    num_vars = get_numvars()
    ct = coeff_table
    a_length = length(a)
    ord = maximum( get_order.(vals) )
    suma = Taylor1(zeros(T, ord))

    for homPol in 1:length(a)
        sun = zero(Taylor1{T})
        for (i,a_coeff) in enumerate(a.coeffs[homPol].coeffs)
            tmp = vals[1]^(ct[homPol][i][1])
            for n in 2:num_vars
                tmp *= vals[n]^(ct[homPol][i][n])
            end
            suma += a_coeff * tmp
        end
    end

    return suma
end

evaluate{T<:Number}(a::TaylorN{T}) = a[1][1]

function evaluate{T<:Number}(x::Array{TaylorN{T},1}, δx::Array{T,1})
    x0 = Array{T}( length(x) )
    evaluate!( x, δx, x0 )
    return x0
end

function evaluate{T<:NumberNotSeries}(x::Array{TaylorN{T},1}, δx::Array{Taylor1{T},1})
    x0 = Array{Taylor1{T}}( length(x) )
    evaluate!( x, δx, x0 )
    return x0
end

function evaluate{T<:NumberNotSeries}(x::Array{TaylorN{Taylor1{T}},1}, δx::Array{Taylor1{T},1})
    x0 = Array{Taylor1{T}}( length(x) )
    evaluate!( x, δx, x0 )
    return x0
end

evaluate{T<:Number}(x::Array{TaylorN{T},1}) = evaluate.(x)

#function-like behavior for TaylorN
(p::TaylorN)(x) = evaluate(p, x)

(p::TaylorN)() = evaluate(p)

#function-like behavior for Array{TaylorN,1}
(p::Array{TaylorN{T},1}){T<:Number}(x) = evaluate(p, x)

(p::Array{TaylorN{T},1}){T<:Number}() = evaluate(p)

