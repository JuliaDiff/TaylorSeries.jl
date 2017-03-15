# This file is part of the TaylorSeries.jl Julia package, MIT license

# Printing of TaylorSeries objects

# subscriptify is taken from the ValidatedNumerics.jl package, licensed under MIT "Expat".
# superscriptify is a small variation

const subscript_digits = [c for c in "₀₁₂₃₄₅₆₇₈₉"]
const superscript_digits = [c for c in "⁰¹²³⁴⁵⁶⁷⁸⁹"]

function subscriptify(n::Int)
    dig = reverse(digits(n))
    join([subscript_digits[i+1] for i in dig])
end

function superscriptify(n::Int)
    dig = reverse(digits(n))
    join([superscript_digits[i+1] for i in dig])
end


function pretty_print{T<:Number}(a::Taylor1{T})
    z = zero(a[1])
    space = string(" ")
    bigO = string("+ 𝒪(t", superscriptify(a.order+1), ")")
    a == zero(a) && return string(space, z, space, bigO)
    strout::String = space
    ifirst = true
    for i in eachindex(a.coeffs)
        monom::String = i==1 ? string("") : i==2 ? string(" t") :
            string(" t", superscriptify(i-1))
        @inbounds c = a[i]
        c == z && continue
        cadena = numbr2str(c, ifirst)
        strout = string(strout, cadena, monom, space)
        ifirst = false
    end
    strout = strout * bigO
    strout
end
function pretty_print{T<:Number}(a::Taylor1{HomogeneousPolynomial{T}})
    z = zero(a[1])
    space = string(" ")
    bigO = string("+ 𝒪(t", superscriptify(a.order+1), ")")
    a == zero(a) && return string(space, z, space, bigO)
    strout::String = space
    ifirst = true
    for i in eachindex(a.coeffs)
        monom::String = i==1 ? string("") : i==2 ? string(" t") :
            string(" t", superscriptify(i-1))
        @inbounds c = a[i]
        c == z && continue
        cadena = numbr2str(c, ifirst)
        ccad::String = i==1 ? cadena : ifirst ? string("(", cadena, ")") :
            string(cadena[1:2], "(", cadena[3:end], ")")
        strout = string(strout, ccad, monom, space)
        ifirst = false
    end
    strout = strout * bigO
    strout
end
function pretty_print{T<:Number}(a::Taylor1{TaylorN{T}})
    z = zero(a[1])
    space = string(" ")
    bigO = string("+ 𝒪(t", superscriptify(a.order+1), ")")
    a == zero(a) && return string(space, z, space, bigO)
    strout::String = space
    ifirst = true
    for i in eachindex(a.coeffs)
        monom::String = i==1 ? string("") : i==2 ? string(" t") :
            string(" t", superscriptify(i-1))
        @inbounds c = a[i]
        c == z && continue
        cadena = numbr2str(c, ifirst)
        ccad::String = i==1 ? cadena : ifirst ? string("(", cadena, ")") :
            string(cadena[1:2], "(", cadena[3:end], ")")
        strout = string(strout, ccad, monom, space)
        ifirst = false
    end
    strout = strout * bigO
    strout
end
function pretty_print{T<:Number}(a::HomogeneousPolynomial{T})
    z = zero(a[1])
    space = string(" ")
    a == zero(a) && return string(space, z)
    strout::String = homogPol2str(a)
    strout
end
function pretty_print{T<:Number}(a::TaylorN{T})
    z = zero(a[1])
    space = string("")
    bigO::String  = string(" + 𝒪(‖x‖", superscriptify(a.order+1), ")")
    a == zero(a) && return string(space, z, bigO)
    strout::String = space#string("")
    ifirst = true
    for ord in eachindex(a.coeffs)
        pol = a[ord]
        pol == zero(a[ord]) && continue
        cadena::String = homogPol2str( pol )
        strsgn = (ifirst || ord == 1 || cadena[2] == '-') ?
            string("") : string(" +")
        strout = string( strout, strsgn, cadena)
        ifirst = false
    end
    strout = strout * bigO
    strout
end

function homogPol2str{T<:Number}(a::HomogeneousPolynomial{T})
    numVars = get_numvars()
    order = a.order
    z = zero(T)
    space = string(" ")
    strout::String = space
    ifirst = true
    iIndices = zeros(Int, numVars)
    for pos = 1:size_table[order+1]
        monom::String = string("")
        @inbounds iIndices[:] = coeff_table[order+1][pos]
        for ivar = 1:numVars
            powivar = iIndices[ivar]
            if powivar == 1
                monom = string(monom, name_taylorNvar(ivar))
            elseif powivar > 1
                monom = string(monom, name_taylorNvar(ivar), superscriptify(powivar))
            end
        end
        @inbounds c = a[pos]
        c == z && continue
        cadena = numbr2str(c, ifirst)
        strout = string(strout, cadena, monom, space)
        ifirst = false
    end
    return strout[1:end-1]
end
function homogPol2str{T<:Number}(a::HomogeneousPolynomial{Taylor1{T}})
    numVars = get_numvars()
    order = a.order
    z = zero(a[1])
    space = string(" ")
    strout::String = space
    ifirst = true
    iIndices = zeros(Int, numVars)
    for pos = 1:size_table[order+1]
        monom::String = string("")
        @inbounds iIndices[:] = coeff_table[order+1][pos]
        for ivar = 1:numVars
            powivar = iIndices[ivar]
            if powivar == 1
                monom = string(monom, name_taylorNvar(ivar))
            elseif powivar > 1
                monom = string(monom, name_taylorNvar(ivar),
                    superscriptify(powivar))
            end
        end
        @inbounds c = a[pos]
        c == z && continue
        cadena = numbr2str(c, ifirst)
        ccad::String = (pos==1 || ifirst) ? string("(", cadena, ")") :
            string(cadena[1:2], "(", cadena[3:end], ")")
        strout = string(strout, ccad, monom, space)
        ifirst = false
    end
    return strout[1:end-1]
end

function numbr2str(zz, ifirst::Bool=false)
    zz == zero(zz) && return string( zz )
    plusmin = ifelse( ifirst, string(""), string("+ ") )
    return string(plusmin, zz)
end
function numbr2str{T<:Union{AbstractFloat,Integer,Rational}}(zz::T, ifirst::Bool=false)
    zz == zero(T) && return string( zz )
    plusmin = ifelse( zz < zero(T), string("- "),
                ifelse( ifirst, string(""), string("+ ")) )
    return string(plusmin, abs(zz))
end
function numbr2str{T<:Complex}(zz::T, ifirst::Bool=false)
    zT = zero(zz.re)
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

name_taylorNvar(i::Int) = string(" ", get_variable_names()[i])


# summary
summary{T<:Number}(a::Taylor1{T}) = string(a.order, "-order ", typeof(a), ":")
function summary{T<:Number}(a::Union{HomogeneousPolynomial{T}, TaylorN{T}})
    string(a.order, "-order ", typeof(a), " in ", get_numvars(), " variables:")
end

# show
function show(io::IO, a::Union{Taylor1, HomogeneousPolynomial, TaylorN})
    print(io, pretty_print(a))
end
