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
    @inbounds suma = a.coeffs[end]
    @inbounds for k = a.order:-1:1
        suma = suma*dx + a.coeffs[k]
    end
    suma
end
function evaluate{T<:Number,S<:Number}(a::Taylor1{T}, dx::S)
    R = promote_type(T,S)
    @inbounds suma = convert(R, a.coeffs[end])
    @inbounds for k = a.order:-1:1
        suma = suma*dx + a.coeffs[k]
    end
    suma
end
evaluate{T<:Number}(a::Taylor1{T}) = a.coeffs[1]

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
evaluate{T<:Number}(a::Array{Taylor1{T},1}) = evaluate(a, zero(T))

doc"""
    evaluate!(x, δt, x0)

Evaluates each element of `x::Array{Taylor1{T},1}`,
representing the Taylor expansion for the dependent variables
of an ODE at *time* δt. It updates the vector `x0` with the
computed values.
"""
function evaluate!{T<:Number}(x::Array{Taylor1{T},1}, δt::T, x0::Array{T,1})
    @assert length(x) == length(x0)
    @inbounds for i in eachindex(x)
        x0[i] = evaluate( x[i], δt )
    end
    nothing
end
function evaluate!{T<:Number, S<:Number}(x::Array{Taylor1{T},1}, δt::S, x0::Array{T,1})
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
    a, x = fixorder(a, x)
    @inbounds suma = a.coeffs[end]
    @inbounds for k = a.order:-1:1
        suma = suma*x + a.coeffs[k]
    end
    suma
end



## Evaluation
"""
    evaluate(a, vals)

Evaluate a `HomogeneousPolynomial` polynomial using Horner's rule (hand coded)
at `vals`.
"""
# function evaluate{T<:Number}(a::HomogeneousPolynomial{T}, vals::Array{T,1} )
function evaluate{T<:Number,S<:NumberNotSeriesN}(a::HomogeneousPolynomial{T},
        vals::Array{S,1} )
    @assert length(vals) == get_numvars()
    R = promote_type(T,S)
    numVars = get_numvars()
    suma = convert(TaylorN{R}, a)

    @inbounds for nv = 1:numVars
        suma = horner(suma, (nv, vals[nv]))
    end

    return suma.coeffs[1].coeffs[1]
end
evaluate(a::HomogeneousPolynomial) = zero(a.coeffs[1])

"""
    evaluate(a, [vals])

Evaluate the `TaylorN` polynomial `a` using Horner's rule (hand coded) at `vals`.
If `vals` is ommitted, it is evaluated at zero.
"""
function evaluate{T<:Number,S<:NumberNotSeriesN}(a::TaylorN{T}, vals::Array{S,1} )
    @assert length(vals) == get_numvars()
    numVars = get_numvars()
    R = promote_type(T,S)
    suma = convert(TaylorN{R}, a)

    @inbounds for nv = 1:numVars
        suma = horner(suma, (nv, vals[nv]))
    end

    return suma.coeffs[1].coeffs[1]
end

evaluate{T<:Number}(a::TaylorN{T}) = a.coeffs[1].coeffs[1]

function evaluate!{T<:Number}(x::Array{TaylorN{T},1}, δx::Array{T,1},
        x0::Array{T,1})
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

## Evaluates HomogeneousPolynomials and TaylorN on a val of the nv variable
## using Horner's rule on the nv variable
function horner{T<:Number,S<:NumberNotSeriesN}(a::HomogeneousPolynomial{T},
        b::Tuple{Int,S} )
    nv, val = b
    numVars = get_numvars()
    @assert 1 <= nv <= numVars
    R = promote_type(T,S)
    @inbounds indTb = coeff_table[a.order+1]
    suma = TaylorN(zero(R), a.order)

    # Horner's rule on the nv variable
    for ord = a.order : -1 : 0
        suma_ord = TaylorN(zero(R), a.order)
        posOrd = order_posTb(a.order, nv, ord)
        neworder = a.order-ord
        for pos in posOrd
            c = a.coeffs[pos]
            iIndices = copy(indTb[pos])
            iIndices[nv] = 0
            kdic = in_base(get_order(), iIndices)
            newpos = pos_table[neworder+1][kdic]
            zhp = HomogeneousPolynomial([zero(R)], neworder)
            zhp.coeffs[newpos] = a.coeffs[pos]
            suma_ord += TaylorN(zhp, a.order)
        end
        if ord == a.order
            suma += suma_ord
        else
            suma = suma*val + suma_ord
        end
    end

    return suma
end

function horner{T<:Number,S<:NumberNotSeriesN}(a::TaylorN{T}, b::Tuple{Int,S} )

    nv, val = b
    @assert 1 <= nv <= get_numvars()
    R = promote_type(T,S)

    suma = TaylorN(zero(R), a.order)
    for ord = a.order:-1:0
        @inbounds polH = a.coeffs[ord+1]
        suma += horner( polH, b)
    end
    suma
end
