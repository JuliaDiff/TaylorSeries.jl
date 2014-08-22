# TaylorSeries.jl
#
# Julia module for handling Taylor series of arbitrary but finite order
#
# utils_Taylor1.jl contains the constructors and methods for 1-variable expansions
#
# utils_TaylorN.jl contains the constructors and methods for N-variable expansions
#
# Last modification: 2014.06.07
#
# Luis Benet & David P. Sanders
# UNAM
#

module TaylorSeries

import Base: zero, one, zeros, ones
import Base: convert, promote_rule, promote, eltype, length, show
import Base: real, imag, conj, ctranspose
import Base: rem, mod, mod2pi
import Base: sqrt, exp, log, sin, cos, tan#, square

abstract AbstractSeries{T<:Number,N} <: Number

include("utils_Taylor1.jl")

include("utils_TaylorN.jl")
gc()

## The following routines combine Taylor and TaylorN, so they must appear defining 
##   Taylor and TaylorN and some of its functionalities

# infostr
infostr{T<:Number}(a::Taylor{T}) = 
    string(a.order, "-order Taylor{", T, "}:\n")
infostr{T<:Number}(a::HomogPol{T}) = 
    string(a.order, "-order HomogPol{", T, "} in ", NUMVARS[end], " variables:\n")
infostr{T<:Number}(a::TaylorN{T}) = 
    string(a.order, "-order TaylorN{", T, "} in ", NUMVARS[end], " variables:\n")

# pretty_print
function pretty_print{T<:Number}(a::Taylor{T})
    print( infostr(a) )
    z = zero(T)
    space = string(" ")
    a == zero(a) && (println(string( space, z)); return)
    strout::ASCIIString = space
    ifirst = true
    for i = 0:a.order
        monom::ASCIIString = i==0 ? string("") : i==1 ? string(" * x_{0}") : string(" * x_{0}^", i)
        @inbounds c = a.coeffs[i+1]
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
    strout::ASCIIString = string("")
    ifirst = true
    for ord = 0:a.order
        pol = a.coeffs[ord+1]
        pol == zero(a.coeffs[ord+1]) && continue
        cadena::ASCIIString = homogPol2str( pol )
        strsgn = (ifirst || ord == 0 || cadena[2] == '-') ? string("") : string(" +")
        strout = string( strout, strsgn, cadena)
        ifirst = false
    end
    println(strout)
    return
end
function pretty_print{T<:Number}(a::Array{Taylor{T},1})
    for i=1:length(a)
        pretty_print(a[i])
    end
    return
end
function pretty_print{T<:Number}(a::Array{TaylorN{T},1})
    for i=1:length(a)
        pretty_print(a[i])
    end
    return
end

# Aux functions related to pretty_print
function homogPol2str{T<:Number}(a::HomogPol{T})
    numVars = NUMVARS[end]
    order = a.order
    varstring = ASCIIString[]
    z = zero(T)
    space = string(" ")
    for ivar = 1:numVars
        push!(varstring, string(" * x_{", ivar, "}"))
    end
    strout::ASCIIString = space
    ifirst = true
    iIndices = zeros(Int, numVars)
    for pos = 1:sizeTable[order+1]
        monom::ASCIIString = string("")
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
function numbr2str{T<:Real}(zz::Complex{T}, ifirst::Bool=false)
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

## Exports to Taylor, TaylorN and HomogPol ##
export Taylor, diffTaylor, integTaylor, evalTaylor, deriv, pretty_print
export TaylorN, HomogPol
export set_maxOrder, get_maxOrder, set_numVars, get_numVars
export taylorvar, ∇, jacobian, hessian

end
