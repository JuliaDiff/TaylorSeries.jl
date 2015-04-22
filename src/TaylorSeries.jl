# TaylorSeries.jl
#
# Julia module for handling Taylor series of arbitrary but finite order
#
# - utils_Taylor1.jl contains the constructors and methods for 1-variable expansions
#
# - utils_TaylorN.jl contains the constructors and methods for N-variable expansions
#
# Last modification: 2015.04.22
#
# Luis Benet & David P. Sanders
# UNAM
#

module TaylorSeries

## Documentation
if VERSION < v"0.4.0-dev"
    using Docile
end

## Compatibility v0.3 -> 0.4
using Compat
@compat sizehint!
@compat trunc
@compat eachindex

import Base: zero, one, zeros, ones,
    convert, promote_rule, promote, eltype, length, show,
    real, imag, conj, ctranspose,
    rem, mod, mod2pi,
    sqrt, exp, log, sin, cos, tan#, square


## Exported types
export AbstractSeries, Taylor, TaylorN, HomogPol
## Exported methods
export diffTaylor, integTaylor, evalTaylor, deriv, pretty_print,
    set_ParamsTaylorN, show_ParamsTaylorN,
        set_maxOrder, get_maxOrder, set_numVars, get_numVars,
    taylorvar, ∇, jacobian, hessian


@doc "The main overall abstract type in TaylorSeries" ->
abstract AbstractSeries{T<:Number,N} <: Number

include("utils_Taylor1.jl")
include("hashtables.jl")
include("utils_TaylorN.jl")


## The following routines combine Taylor and TaylorN, so they must appear defining
##   Taylor and TaylorN and some of its functionalities

# infostr
infostr{T<:Number}(a::Taylor{T}) =
    string(a.order, "-order Taylor{", T, "}:\n")
infostr{T<:Number}(a::HomogPol{T}) =
    string(a.order, "-order HomogPol{", T, "} in ", _params.numVars, " variables:\n")
infostr{T<:Number}(a::TaylorN{T}) =
    string(a.order, "-order TaylorN{", T, "} in ", _params.numVars, " variables:\n")

# pretty_print
function pretty_print{T<:Number}(a::Taylor{T})
    print( infostr(a) )
    z = zero(T)
    space = utf8(" ")
    a == zero(a) && (println(string( space, z)); return)
    strout::UTF8String = space
    ifirst = true
    for i in eachindex(a.coeffs)
        monom::UTF8String = i==1 ? string("") : i==2 ?
            string("⋅x_{0}") : string("⋅x_{0}^", i-1)
        @inbounds c = a.coeffs[i]
        c == z && continue
        cadena = numbr2str(c, ifirst)
        strout = string(strout, cadena, monom, space)
        ifirst = false
    end
    println(strout)
    return
end
function pretty_print{T<:Number}(a::HomogPol{T})
    print( infostr(a) )
    z = zero(T)
    a == zero(a) && (println(string( " ", z)); return)
    strout = homogPol2str(a)
    println(strout)
    return
end
function pretty_print{T<:Number}(a::TaylorN{T})
    print( infostr(a) )
    a == zero(a) && (println(string( " ", zero(T))); return)
    strout::UTF8String = string("")
    ifirst = true
    for ord in eachindex(a.coeffs)
        pol = a.coeffs[ord]
        pol == zero(a.coeffs[ord]) && continue
        cadena::UTF8String = homogPol2str( pol )
        strsgn = (ifirst || ord == 1 || cadena[2] == '-') ? string("") : string(" +")
        strout = string( strout, strsgn, cadena)
        ifirst = false
    end
    println(strout)
    return
end
function pretty_print{T<:Number}(a::Array{Taylor{T},1})
    for i in eachindex(a)
        pretty_print(a[i])
    end
    return
end
function pretty_print{T<:Number}(a::Array{TaylorN{T},1})
    for i in eachindex(a)
        pretty_print(a[i])
    end
    return
end

# Aux functions related to pretty_print
function homogPol2str{T<:Number}(a::HomogPol{T})
    numVars = _params.numVars
    order = a.order
    varstring = UTF8String[]
    z = zero(T)
    space = utf8(" ")
    for ivar = 1:numVars
        push!(varstring, string(" ⋅ x_{", ivar, "}"))
    end
    strout::UTF8String = space
    ifirst = true
    iIndices = zeros(Int, numVars)
    for pos = 1:sizeTable[order+1]
        monom::UTF8String = string("")
        @inbounds iIndices[:] = indicesTable[order+1][pos]
        for ivar = 1:numVars
            powivar = iIndices[ivar]
            if powivar == 1
                monom = string(monom, varstring[ivar])
            elseif powivar > 1
                monom = string(monom, varstring[ivar], "^", powivar)
            end
        end
        @inbounds c = a.coeffs[pos]
        c == z && continue
        cadena = numbr2str(c, ifirst)
        strout = string(strout, cadena, monom, space)
        ifirst = false
    end
    return strout[1:end-1]
end
function numbr2str{T<:Real}(zz::T, ifirst::Bool=false)
    zz == zero(T) && return string( zz )
    plusmin = zz > zero(T) ? string("+ ") : string("- ")
    if ifirst
        plusmin = zz > zero(T) ? string("") : string("- ")
    end
    return string(plusmin, abs(zz))
end
function numbr2str(zz::Complex, ifirst::Bool=false)
    T = typeof(zz.re)
    zT = zero(T)
    zz == zero(Complex{T}) && return zT
    zre, zim = reim(zz)
    cadena = string("")
    if zre > zT
        if ifirst
            cadena = string(" ( ", abs(zre)," ")
        else
            cadena = string(" + ( ", abs(zre)," ")
        end
        if zim > zT
            cadena = string(cadena, "+ ", abs(zim), " im )")
        elseif zim < zT
            cadena = string(cadena, "- ", abs(zim), " im )")
        else
            cadena = string(cadena, ")")
        end
    elseif zre < zT
        cadena = string(" - ( ", abs(zre), " ")
        if zim > zT
            cadena = string(cadena, "- ", abs(zim), " im )")
        elseif zim < zT
            cadena = string(cadena, "+ ", abs(zim), " im )")
        else
            cadena = string(cadena, ")")
        end
    else
        if zim > zT
            if ifirst
                cadena = string("( ", abs(zim), " im )")
            else
                cadena = string("+ ( ", abs(zim), " im )")
            end
        else
            cadena = string("- ( ", abs(zim), " im )")
        end
    end
    return cadena
end

end
