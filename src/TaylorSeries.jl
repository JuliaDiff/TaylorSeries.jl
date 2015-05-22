# TaylorSeries.jl
#
# Julia module for handling Taylor series of arbitrary but finite order
#
# - utils_Taylor1.jl contains the constructors and methods for 1-variable expansions
#
# - utils_TaylorN.jl contains the constructors and methods for N-variable expansions
#
# Last modification: 2015.05.08
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
@compat round

import Base: zero, one, zeros, ones,
    convert, promote_rule, promote, eltype, length, show,
    real, imag, conj, ctranspose,
    rem, mod, mod2pi, gradient,
    sqrt, exp, log, sin, cos, tan


## Exported types and methods
export Taylor1, TaylorN, HomogeneousPolynomial
export taylor1_variable, taylorN_variable, get_coeff,
    diffTaylor, integTaylor, evalTaylor, deriv,
    set_params_TaylorN, show_params_TaylorN,
    set_maxOrder, get_maxOrder, set_numVars, get_numVars,
    ∇, jacobian, hessian


include("utils_Taylor1.jl")
include("hashtables.jl")
include("utils_TaylorN.jl")


# pretty_print
function pretty_print{T<:Number}(a::Taylor1{T})
    z = zero(T)
    space = utf8(" ")
    a == zero(a) && return string(space, z)
    strout::UTF8String = space
    ifirst = true
    for i in eachindex(a.coeffs)
        monom::UTF8String = i==1 ? string("") : i==2 ? string("⋅t") :
            string("⋅t", superscriptify(i-1))
        @inbounds c = a.coeffs[i]
        c == z && continue
        cadena = numbr2str(c, ifirst)
        strout = string(strout, cadena, monom, space)
        ifirst = false
    end
    strout
end
function pretty_print{T<:Number}(a::HomogeneousPolynomial{T})
    z = zero(T)
    space = utf8(" ")
    a == zero(a) && return string(space, z)
    strout::UTF8String = homogPol2str(a)
    strout
end
function pretty_print{T<:Number}(a::TaylorN{T})
    z = zero(T)
    space = utf8(" ")
    a == zero(a) && return string(space, z)
    strout::UTF8String = utf8("")
    ifirst = true
    for ord in eachindex(a.coeffs)
        pol = a.coeffs[ord]
        pol == zero(a.coeffs[ord]) && continue
        cadena::UTF8String = homogPol2str( pol )
        strsgn = (ifirst || ord == 1 || cadena[2] == '-') ? string("") : string(" +")
        strout = string( strout, strsgn, cadena)
        ifirst = false
    end
    strout
end

function homogPol2str{T<:Number}(a::HomogeneousPolynomial{T})
    numVars = _params_taylorN.numVars
    order = a.order
    z = zero(T)
    space = utf8(" ")
    strout::UTF8String = space
    ifirst = true
    iIndices = zeros(Int, numVars)
    for pos = 1:sizeTable[order+1]
        monom::UTF8String = string("")
        @inbounds iIndices[:] = indicesTable[order+1][pos]
        for ivar = 1:numVars
            powivar = iIndices[ivar]
            if powivar == 1
                monom = string(monom, name_taylorNvar(ivar))
            elseif powivar > 1
                monom = string(monom, name_taylorNvar(ivar), superscriptify(powivar))
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
    zz == zero(zz) && return string(zT)
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

name_taylorNvar(n::Int) = string("⋅x", subscriptify(n))

# subscriptify is taken from ValidatedNumerics/src/nterval_definition.jl
# and is licensed under MIT "Expat".
# superscriptify is a small variation
function subscriptify(n::Int)
    subscript_digits = [c for c in "₀₁₂₃₄₅₆₇₈₉"]
    dig = reverse(digits(n))
    join([subscript_digits[i+1] for i in dig])
end
function superscriptify(n::Int)
    superscript_digits = [c for c in "⁰¹²³⁴⁵⁶⁷⁸⁹"]
    dig = reverse(digits(n))
    join([superscript_digits[i+1] for i in dig])
end

# summary
summary{T<:Number}(a::Taylor1{T}) = string(a.order, "-order ", typeof(a), ":")
function summary{T<:Number}(a::Union(HomogeneousPolynomial{T}, TaylorN{T}))
    string(a.order, "-order ", typeof(a), " in ", _params_taylorN.numVars, " variables:")
end

# show
function show(io::IO, a::Union(Taylor1, HomogeneousPolynomial, TaylorN))
    # (isa(a, TaylorN) || isa(a, Taylor1)) && println(io, summary(a))
    print(io, pretty_print(a))
end

end # module
