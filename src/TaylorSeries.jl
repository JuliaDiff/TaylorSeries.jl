# TaylorSeries.jl
#
# Julia module for handling Taylor series of arbitrary but finite order
#
# utils_Taylor1.jl contains the constructors and methods for 1-variable expansions
#
# utils_TaylorN.jl contains the constructors and methods for N-variable expansions
#
# Last modification: 2014.04.12
#
# Luis Benet & David P. Sanders
# UNAM
#

module TaylorSeries

import Base: zero, one
import Base: convert, promote_rule, promote, eltype, length, pretty_print
import Base: real, imag, conj, ctranspose
import Base: rem, mod, mod2pi
import Base: sqrt, exp, log, sin, cos, tan#, square

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
    space = " "
    #a == z && return string( space, z)
    a == z && (print(string( space, z)); return)
    strout = space
    ifirst = true
    for i = 0:a.order
        monom = i==0 ? "" : i==1 ? " * x0" : string(" * x0^", i)
        @inbounds c = a.coeffs[i+1]
        c == z && continue
        cadena = numbr2str(c, ifirst)
        strout = string(strout, cadena, monom, space)
        ifirst = false
    end
    println(strout)
    #return strout
end
function pretty_print{T<:Number}(a::TaylorN{T})
    print( infostr(a) )
    z = zero(T)
    space = " "
    a == z && (print(string( space, z)); return)
    varstring = {}
    for ivar=1:a.numVars
        push!(varstring,string(" * x",ivar))
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
    #return strout
end

# make string from a number; for complex numbers, use 
function numbr2str{T<:Real}(zz::T, ifirst=false::Bool)
    plusmin = zz > 0 ? "+ " : "- "
    if ifirst
        plusmin = zz > 0 ? "" : "-"
    end
    return string(plusmin, abs(zz))
end
function numbr2str{T}(zz::Complex{T}, ifirst=false::Bool)
    zT = zero(T)
    zz == zero(Complex{T}) && return zT
    zre, zim = reim(zz)
    cadena = ""
    if zre > zT
        if ifirst
            cadena = string("( ", abs(zre)," ")
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
                cadena = string(" + ( ", abs(zim), " im )")
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
export Taylor, diffTaylor, integTaylor, evalTaylor, deriv
#
export TaylorN
export set_maxOrder, get_maxOrder, set_numVars, get_numVars

end