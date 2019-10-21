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


function pretty_print(a::Taylor1{T}) where {T<:NumberNotSeries}
    z = zero(a[0])
    var = _params_Taylor1_.var_name
    space = string(" ")
    bigO = bigOnotation[end] ?
        string("+ 𝒪(", var, superscriptify(a.order+1), ")") :
        string("")
    iszero(a) && return string(space, z, space, bigO)
    strout::String = space
    ifirst = true
    for i in eachindex(a)
        monom::String = i==0 ? string("") : i==1 ? string(" ", var) :
            string(" ", var, superscriptify(i))
        @inbounds c = a[i]
        c == z && continue
        cadena = numbr2str(c, ifirst)
        strout = string(strout, cadena, monom, space)
        ifirst = false
    end
    strout = strout * bigO
    strout
end

function pretty_print(a::Taylor1{T} where {T <: AbstractSeries{S}}) where {S<:Number}
    z = zero(a[0])
    var = _params_Taylor1_.var_name
    space = string(" ")
    bigO = bigOnotation[end] ?
        string("+ 𝒪(", var, superscriptify(a.order+1), ")") :
        string("")
    iszero(a) && return string(space, z, space, bigO)
    strout::String = space
    ifirst = true
    for i in eachindex(a)
        monom::String = i==0 ? string("") : i==1 ? string(" ", var) :
            string(" ", var, superscriptify(i))
        @inbounds c = a[i]
        c == z && continue
        cadena = numbr2str(c, ifirst)
        ccad::String = i==0 ? cadena : ifirst ? string("(", cadena, ")") :
            string(cadena[1:2], "(", cadena[3:end], ")")
        strout = string(strout, ccad, monom, space)
        ifirst = false
    end
    strout = strout * bigO
    strout
end

function pretty_print(a::HomogeneousPolynomial{T}) where {T<:Number}
    z = zero(a[1])
    space = string(" ")
    iszero(a) && return string(space, z)
    strout::String = homogPol2str(a)
    strout
end

function pretty_print(a::TaylorN{T}) where {T<:Number}
    z = zero(a[0])
    space = string("")
    bigO::String = bigOnotation[end] ?
        string(" + 𝒪(‖x‖", superscriptify(a.order+1), ")") :
        string("")
    iszero(a) && return string(space, z, space, bigO)
    strout::String = space
    ifirst = true
    for ord in eachindex(a)
        pol = a[ord]
        iszero(pol) && continue
        cadena::String = homogPol2str( pol )
        strsgn = (ifirst || ord == 0 || cadena[2] == '-') ?
            string("") : string(" +")
        strout = string( strout, strsgn, cadena)
        ifirst = false
    end
    strout = strout * bigO
    strout
end

function homogPol2str(a::HomogeneousPolynomial{T}) where {T<:Number}
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
    return strout[1:prevind(strout, end)]
end

function homogPol2str(a::HomogeneousPolynomial{Taylor1{T}}) where {T<:Number}
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
    return strout[1:prevind(strout, end)]
end

function numbr2str(zz, ifirst::Bool=false)
    iszero(zz) && return string( zz )
    plusmin = ifelse( ifirst, string(""), string("+ ") )
    return string(plusmin, zz)
end

function numbr2str(zz::T, ifirst::Bool=false) where
        {T<:Union{AbstractFloat,Integer,Rational}}
    iszero(zz) && return string( zz )
    plusmin = ifelse( zz < zero(T), string("- "),
                ifelse( ifirst, string(""), string("+ ")) )
    return string(plusmin, abs(zz))
end

function numbr2str(zz::T, ifirst::Bool=false) where {T<:Complex}
    zT = zero(zz.re)
    iszero(zz) && return string(zT)
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
summary(a::Taylor1{T}) where {T<:Number} =
    string(a.order, "-order ", typeof(a), ":")

function summary(a::Union{HomogeneousPolynomial{T}, TaylorN{T}}) where {T<:Number}
    string(a.order, "-order ", typeof(a), " in ", get_numvars(), " variables:")
end

# show
function show(io::IO, a::Union{Taylor1, HomogeneousPolynomial, TaylorN})
    if _show_default[end]
        return Base.show_default(IOContext(io, :compact => false), a)
    else
        return print(io, pretty_print(a))
    end
end
