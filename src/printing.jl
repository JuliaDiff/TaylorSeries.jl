
# pretty_print
function pretty_print{T<:Number}(a::Taylor1{T})
    z = zero(T)
    space = utf8(" ")
    bigO = string("+ 𝒪(t", superscriptify(a.order+1), ")")
    a == zero(a) && return string(space, z, space, bigO)
    strout::UTF8String = space
    ifirst = true
    for i in eachindex(a.coeffs)
        monom::UTF8String = i==1 ? string("") : i==2 ? string(" t") :
            string(" t", superscriptify(i-1))
        @inbounds c = a.coeffs[i]
        c == z && continue
        cadena = numbr2str(c, ifirst)
        strout = string(strout, cadena, monom, space)
        ifirst = false
    end
    strout = strout * bigO
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
    bigO::UTF8String  = string(" + 𝒪(‖x‖", superscriptify(a.order+1), ")")
    a == zero(a) && return string(space, z, bigO)
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
    strout = strout * bigO
    strout
end

function power_string(i, exponent)
    exponent == 0 && return ""
    exponent == 1 && return name_taylorN_var(i)

    string(name_taylorN_var(i), superscriptify(exponent))
end

# Replace string construction by IOBuffer
function homogPol2str{T<:Number}(a::HomogeneousPolynomial{T})
    num_vars = get_numvars()
    order = a.order
    z = zero(T)
    space = utf8(" ")
    output::UTF8String = space
    ifirst = true

    for (pos, indices) in enumerate(index_table[order+1])
        c = a.coeffs[pos]
        c == zero(T) && continue

        monom = string([power_string(i, indices[i]) for i in 1:num_vars]...)
        number = numbr2str(c, ifirst)
        output = string(output, number, monom, space)
        ifirst = false
    end
    return output[1:end-1]
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

name_taylorN_var(i::Int) = string(" ", _params_TaylorN_.variable_names[i])


# summary
summary{T<:Number}(a::Taylor1{T}) = string(a.order, "-order ", typeof(a), ":")
function summary{T<:Number}(a::Union(HomogeneousPolynomial{T}, TaylorN{T}))
    string(a.order, "-order ", typeof(a), " in ", get_numvars(), " variables:")
end

# show
function show(io::IO, a::Union(Taylor1, HomogeneousPolynomial, TaylorN))
    # (isa(a, TaylorN) || isa(a, Taylor1)) && println(io, summary(a))
    print(io, pretty_print(a))
end

