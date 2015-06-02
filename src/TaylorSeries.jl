# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Handles Taylor series of arbitrary but finite order

module TaylorSeries

if VERSION < v"0.4.0-dev"
    using Docile
end

using Compat

import Base: ==, +, -, *, /, ^
import Base: zero, one, zeros, ones,
    convert, promote_rule, promote, eltype, length, show,
    real, imag, conj, ctranspose,
    rem, mod, mod2pi, gradient,
    sqrt, exp, log, sin, cos, tan


export Taylor1, TaylorN, HomogeneousPolynomial
export taylor1_variable, taylorN_variable, get_coeff,
    diffTaylor, integTaylor, evaluate, deriv,
    show_params_TaylorN,
    get_order, get_numvars,
    set_variables, get_variables,
    ∇, jacobian, hessian
export AbstractSeries

@doc "The main overall abstract type in TaylorSeries" ->
abstract AbstractSeries{T<:Number,N}# <: Number

# one variable:
include("Taylor1.jl")

# several variables:
include("parameters.jl")
include("hash_tables.jl")
include("TaylorN.jl")

include("printing.jl")

#-------
function pretty_print(a::TaylorRec)
    z = zero(eltype(a))
    space = utf8(" ")
    a == zero(a) && return string(space, z)
    strout::UTF8String = ""
    ifirst = true
    for i in eachindex(a.coeffs)
        pol = a.coeffs[i]
        pol == zero(pol) && continue
        cadena::UTF8String = taylorrec2str( pol, i-1, 1 )
        strsgn = (ifirst || i == 1 || cadena[2] == '-') ? string("") : string(" +")
        strout = string( strout, strsgn, cadena)
        ifirst = false
    end
    strout
end

function monomial(ord::Int, nv::Int)
    monom::UTF8String = string("")
    if ord == 0
        monom = string(monom)
    elseif ord == 1
        monom = string(monom, name_taylorNvar(nv))
    else
        monom = string(monom, name_taylorNvar(nv), superscriptify(ord))
    end
    monom
end

function taylorrec2str{T<:Number}(c::T, ord::Int, nv::Int)
    z = zero(T)
    space = utf8(" ")
    c == z && return string( space, z)
    strout::UTF8String = space
    ifirst = true

    monom = monomial(ord,nv)
    cadena = numbr2str(c, ifirst)
    strout = string(strout, cadena, monom)

    return strout
end

function taylorrec2str{T<:Number,N}(a::TaylorRec{T,N}, ord::Int, nv::Int)
    z = zero(T)
    space = utf8(" ")
    a == zero(a) && return string( space, z)
    strout::UTF8String = string("")
    ifirst = true

    monom = monomial(ord,nv)
    for i in eachindex(a.coeffs)
        pol = a.coeffs[i]
        pol == zero(pol) && continue
        cadena::UTF8String = taylorrec2str( pol, i-1, nv+1 )
        strsgn = (ifirst || i == 1 || cadena[2] == '-') ? string("") : string(" +")
        strout = string( strout, strsgn, cadena)
        ifirst = false
    end
    strout = string(space,"(",strout," )", monom)
    strout
end
#-------

end # module
