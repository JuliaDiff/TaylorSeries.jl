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

## Compatibility v0.3 -> 0.4
using Compat
@compat trunc
@compat Tuple{Int, Int}


import Base: zero, one
import Base: convert, promote_rule, promote, eltype, length
import Base: real, imag, conj, ctranspose
import Base: rem, mod, mod2pi
import Base: sqrt, exp, log, sin, cos, tan#, square

abstract AbstractSeries{T<:Number,N} <: Number

include("utils_Taylor1.jl")

include("utils_TaylorN.jl")

## The following routines combine Taylor and TaylorN, so they must appear defining 
##   Taylor and TaylorN and some of its functionalities

# infostr
infostr{T<:Number}(a::Taylor{T}) = 
    string(a.order, "-order Taylor{", T, "}:\n")
infostr{T<:Number}(a::TaylorN{T}) = 
    string(a.order, "-order TaylorN{", T, "} in ", a.numVars, " variables:\n")

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
function pretty_print{T<:Number}(a::TaylorN{T})
    print( infostr(a) )
    a == zero(a) && (println(string( " ", zero(T))); return)
    z = zero(T)
    space = " "
    varstring = ASCIIString[]
    for ivar=1:a.numVars
        push!(varstring,string(" * x_{", ivar, "}"))
    end
    strout = string(" ")
    ifirst = true
    iIndices = zeros(Int, a.numVars)
    for pos = 1:length(a.coeffs)
        monom = ""
        iIndices = indicesTable[end][pos]
        if pos>1
            for ivar=1:a.numVars
                powivar = iIndices[ivar]
                if powivar == 1
                    monom = string(monom, varstring[ivar])
                elseif powivar > 1
                    monom = string(monom, varstring[ivar], "^", powivar)
                end
            end
        end
        @inbounds c = a.coeffs[pos]
        c == z && continue
        cadena = numbr2str(c, ifirst)
        strout = string(strout, cadena, monom, space)
        ifirst = false
    end
    println(strout)
    return
end
function pretty_print{T<:Number}(a::Union(Array{Taylor{T},1},Array{TaylorN{T},1}))
    for i=1:length(a)
        pretty_print(a[i])
        println("")
    end
    return
end

# make string from a number; for complex numbers, use 
function numbr2str{T<:Real}(zz::T, ifirst::Bool=false)
    zz == zero(T) && return string( zz )
    plusmin = zz > zero(T) ? "+ " : "- "
    if ifirst
        plusmin = zz > zero(T) ? " " : "-"
    end
    return string(plusmin, abs(zz))
end
function numbr2str{T}(zz::Complex{T}, ifirst::Bool=false)
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
            cadena = string(" - ( ", abs(zim), " im )")
        end
    end
    return cadena
end

# evalTaylor(TaylorN, Array{Taylor,1})
function evalTaylor{T<:Number,S<:Number}(a::TaylorN{T}, vT::Array{Taylor{S},1})
    numVars = NUMVARS[end]
    @assert length(vT) == numVars
    R = promote_type(T,S)
    order = length(vT[1])
    sumaT = Taylor(zero(R), order)
    z = zero(sumaT)
    nCoefTot = sizeTable[end]
    iIndices = zeros(Int, numVars)
    for pos = nCoefTot:-1:1
        @inbounds iIndices[1:end] = indicesTable[end][pos]
        a.coeffs[pos] == zero(T) && continue
        val = Taylor(a.coeffs[pos], order)
        for k = 1:numVars
            @inbounds val = val*(vT[k])^iIndices[k]
        end
        sumaT += val
    end
    return sumaT
end


## Exports to Taylor and TaylorN ##
export Taylor, diffTaylor, integTaylor, evalTaylor, deriv, pretty_print
#
export TaylorN
export set_maxOrder, get_maxOrder, set_numVars, get_numVars

end