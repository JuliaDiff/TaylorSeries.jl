# This file is part of the TaylorSeries.jl Julia package, MIT license

# Printing of TaylorSeries objects

# subscriptify is taken from the ValidatedNumerics.jl package, licensed under MIT "Expat".
# superscriptify is a small variation

const subscript_digits = [c for c in "₀₁₂₃₄₅₆₇₈₉"]
const superscript_digits = [c for c in "⁰¹²³⁴⁵⁶⁷⁸⁹"]

function subscriptify(n::Int)
    dig = reverse(digits(n))
    return join([subscript_digits[i+1] for i in dig])
end

function superscriptify(n::Int)
    dig = reverse(digits(n))
    return join([superscript_digits[i+1] for i in dig])
end

var_subscr(var::S, n::Int) where {S<:AbstractString} = string(var, subscriptify(n))
var_supscr(var::S, n::Int) where {S<:AbstractString} = string(var, superscriptify(n))

function str_bigO(var::String, i::Int, bb::Bool = true) :: String
    if bigOnotation[end]
        str = string("+ 𝒪(", var_supscr(var, i+1), ")")
    else
        return string("")
    end
    bb && return str
    return string(" ", str)
end

function monom_string(i::Int, var::String, bb::Bool = true) :: String
    if i==0
        return string("")
    elseif i==1
        str = var
    else
        str = var_supscr(var, i)
    end
    bb && return string(" ", str)
    return string(str)
end



#  Fallback
function pretty_print(a::Taylor1)
    # z = zero(a[0])
    var = _params_Taylor1_.var_name[1]
    space = string(" ")
    bigO = str_bigO(var, a.order)
    TS._isthinzero(a) && return string(space, 0, space, bigO)
    strout::String = space
    ifirst = true
    for i in eachindex(a)
        monom::String = monom_string(i, var)
        @inbounds c = a[i]
        # c == z && continue
        cadena = numbr2str(c, ifirst)
        strout = string(strout, cadena, monom, space)
        ifirst = false
    end
    strout = strout * bigO
    strout
end

function pretty_print(a::Taylor1{T}) where {T<:NumberNotSeries}
    z = zero(a[0])
    var = _params_Taylor1_.var_name[1]
    space = string(" ")
    bigO = str_bigO(var, a.order)
    TS._isthinzero(a) && return string(space, z, space, bigO)
    strout::String = space
    ifirst = true
    for i in eachindex(a)
        monom::String = monom_string(i, var)
        @inbounds c = a[i]
        TS._isthinzero(c) && continue
        cadena = numbr2str(c, ifirst)
        strout = string(strout, cadena, monom, space)
        ifirst = false
    end
    strout = strout * bigO
    return strout
end

function pretty_print(a::Taylor1{<:AbstractSeries})
    z = zero(a[0])
    if get_numvars(a) > length(_params_Taylor1_.var_name)
        set_taylor1_varname(get_numvars(a), _params_Taylor1_.var_name[1][1:1])
    end
    var = _params_Taylor1_.var_name[get_numvars(a)]
    space = string(" ")
    bigO = str_bigO(var, a.order)
    TS._isthinzero(a) && return string(space, z, space, bigO)
    strout::String = space
    ifirst = true
    for i in eachindex(a)
        monom::String = monom_string(i, var)
        @inbounds c = a[i]
        TS._isthinzero(c) && continue
        cadena = numbr2str(c, ifirst)
        ccad::String = i==0 ? cadena : ifirst ? string("(", cadena, ")") :
            string(cadena[1:2], "(", cadena[3:end], ")")
        strout = string(strout, ccad, monom, space)
        ifirst = false
    end
    strout = strout * bigO
    return strout
end

function pretty_print(a::HomogeneousPolynomial{T}) where {T<:Number}
    z = zero(a[1])
    space = string(" ")
    TS._isthinzero(a) && return string(space, z)
    strout::String = homogPol2str(a)
    strout
end

function pretty_print(a::TaylorN{T}) where {T<:Number}
    z = zero(a[0])
    space = string("")
    bigO :: String = str_bigO("‖x‖", a.order, false)
    TS._isthinzero(a) && return string(space, z, space, bigO)
    strout::String = space
    ifirst = true
    for ord in eachindex(a)
        pol = a[ord]
        TS._isthinzero(pol) && continue
        cadena::String = homogPol2str( pol )
        strsgn = (ifirst || ord == 0 || cadena[2] == '-') ? string("") : string(" +")
        strout = string( strout, strsgn, cadena)
        ifirst = false
    end
    strout = strout * bigO
    strout
end

function homogPol2str(a::HomogeneousPolynomial{T}) where {T<:Number}
    numVars = get_numvars()
    order = a.order
    z = zero(a.coeffs[1])
    space = string(" ")
    strout::String = space
    ifirst = true
    iIndices = zeros(Int, numVars)
    for pos = 1:size_table[order+1]
        monom::String = string("")
        @inbounds iIndices[:] = coeff_table[order+1][pos]
        for ivar = 1:numVars
            powivar = iIndices[ivar]
            monom = string(monom, monom_string(powivar, name_taylorNvar(ivar), false))
        end
        @inbounds c = a[pos]
        TS._isthinzero(c) && continue
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
            monom = string(monom, monom_string(powivar, name_taylorNvar(ivar), false))
        end
        @inbounds c = a[pos]
        TS._isthinzero(c) && continue
        cadena = numbr2str(c, ifirst)
        ccad::String = (pos==1 || ifirst) ? string("(", cadena, ")") :
            string(cadena[1:2], "(", cadena[3:end], ")")
        strout = string(strout, ccad, monom, space)
        ifirst = false
    end
    return strout[1:prevind(strout, end)]
end

function numbr2str(zz, ifirst::Bool=false)
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

function numbr2str(zz::Complex, ifirst::Bool=false)
    zT = zero(zz.re)
    iszero(zz) && return string(zT)
    zre, zim = reim(zz)
    if zre > zT
        if ifirst
            cadena = string("( ", zz, " )")
        else
            cadena = string("+ ( ", zz, " )")
        end
    elseif zre < zT
        cadena = string("- ( ", -zz, " )")
    elseif zre == zT
        if zim > zT
            if ifirst
                cadena = string("( ", zz, " )")
            else
                cadena = string("+ ( ", zz, " )")
            end
        elseif zim < zT
            cadena = string("- ( ", -zz, " )")
        else
            if ifirst
                cadena = string("( ", zz, " )")
            else
                cadena = string("+ ( ", zz, " )")
            end
        end
    else
        if ifirst
            cadena = string("( ", zz, " )")
        else
            cadena = string("+ ( ", zz, " )")
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
function show(io::IO, a::AbstractSeries)
    if _show_default[end]
        return Base.show_default(IOContext(io, :compact => false), a)
    else
        return print(io, pretty_print(a))
    end
end
