# This file is part of the Taylor1Series.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

## Conversion
convert{T<:Number}(::Type{Taylor1{T}}, a::Taylor1) =
    Taylor1(convert(Array{T,1}, a.coeffs), a.order)

function convert{T<:Integer, S<:AbstractFloat}(::Type{Taylor1{Rational{T}}},
    a::Taylor1{S})
    v = Array{Rational{T}}(length(a.coeffs))
    for i in eachindex(v)
        v[i] = rationalize(a.coeffs[i])
    end
    Taylor1(v)
end

convert{T<:Number}(::Type{Taylor1{T}}, b::Array{T,1}) = Taylor1(b, length(b)-1)

convert{T<:Number, S<:Number}(::Type{Taylor1{T}}, b::Array{S,1}) =
    Taylor1(convert(Array{T,1},b))

convert{T<:Number, S<:Number}(::Type{Taylor1{T}}, b::S) = Taylor1([convert(T,b)], 0)

convert{T<:Number}(::Type{Taylor1{T}}, b::T) = Taylor1([b], 0)

convert(::Type{Taylor1}, a::Number) = Taylor1(a, 0)


convert{T<:Number}(::Type{HomogeneousPolynomial{T}}, a::HomogeneousPolynomial) =
    HomogeneousPolynomial{T}(convert(Array{T,1}, a.coeffs), a.order)

function convert{T<:Integer, S<:AbstractFloat}(
    ::Type{HomogeneousPolynomial{Rational{T}}}, a::HomogeneousPolynomial{S})
    v = Array{Rational{T}}(length(a.coeffs))
    for i in eachindex(v)
        # v[i] = convert(Rational{T}, a.coeffs[i])
        v[i] = rationalize(a.coeffs[i], tol=eps(one(S)))
    end
    HomogeneousPolynomial(v, a.order)
end
convert{T<:Number, S<:Number}(::Type{HomogeneousPolynomial{T}}, b::Array{S,1}) =
    HomogeneousPolynomial{T}(convert(Array{T,1}, b), orderH(b))

convert{T<:Number, S<:Number}(::Type{HomogeneousPolynomial{T}}, b::S) =
    HomogeneousPolynomial{T}([convert(T,b)], 0)

convert{T<:Number}(::Type{HomogeneousPolynomial{T}}, b::Array{T,1}) =
    HomogeneousPolynomial{T}(b, orderH(b))

convert{T<:Number}(::Type{HomogeneousPolynomial{T}}, b::T) =
    HomogeneousPolynomial{T}([b], 0)

convert(::Type{HomogeneousPolynomial}, a::Number) = HomogeneousPolynomial([a],0)


convert{T<:Number}(::Type{TaylorN{T}}, a::TaylorN) =
    TaylorN{T}( convert(Array{HomogeneousPolynomial{T},1}, a.coeffs), a.order)

convert{T<:Number, S<:Number}(::Type{TaylorN{T}}, b::HomogeneousPolynomial{S}) =
    TaylorN{T}( [convert(HomogeneousPolynomial{T}, b)], b.order)

convert{T<:Number, S<:Number}(::Type{TaylorN{T}},
    b::Array{HomogeneousPolynomial{S},1}) =
    TaylorN{T}( convert(Array{HomogeneousPolynomial{T},1}, b), length(b)-1)

convert{T<:Number, S<:Number}(::Type{TaylorN{T}}, b::S) =
    TaylorN( [HomogeneousPolynomial([convert(T, b)], 0)], 0)

convert{T<:Number}(::Type{TaylorN{T}}, b::HomogeneousPolynomial{T}) =
    TaylorN{T}( [b], b.order)

convert{T<:Number}(::Type{TaylorN{T}}, b::Array{HomogeneousPolynomial{T},1}) =
    TaylorN{T}( b, length(b)-1)

convert{T<:Number}(::Type{TaylorN{T}}, b::T) =
    TaylorN( [HomogeneousPolynomial([b], 0)], 0)

convert(::Type{TaylorN}, b::Number) = TaylorN( [HomogeneousPolynomial([b], 0)], 0)


function convert{T<:Number}(::Type{TaylorN{Taylor1{T}}}, s::Taylor1{TaylorN{T}})
    orderN = get_order()
    r = zeros(HomogeneousPolynomial{Taylor1{T}}, orderN)

    v = zeros(T, s.order+1)
    for ordT in 1:s.order+1
        v[ordT] = one(T)
        for ordHP in eachindex(s.coeffs[ordT].coeffs)
            for ic in eachindex(s.coeffs[ordT].coeffs[ordHP].coeffs)
                coef = s.coeffs[ordT].coeffs[ordHP].coeffs[ic]
                r[ordHP].coeffs[ic] += Taylor1(coef*v)
            end
        end
        v[ordT] = zero(T)
    end
    TaylorN(r)
end

function convert{T<:Number}(::Type{Taylor1{TaylorN{T}}}, s::TaylorN{Taylor1{T}})
    ordert = 0
    for ordHP in eachindex(s.coeffs)
        ordert = max(ordert, s.coeffs[ordHP].coeffs[1].order)
    end
    vT = Array{TaylorN{T}}(ordert+1)
    for ordT in eachindex(vT)
        vT[ordT] = TaylorN(zero(T), s.order)
    end

    for ordN in eachindex(s.coeffs)
        vHP = HomogeneousPolynomial(zeros(T,length(s.coeffs[ordN])))
        for ihp in eachindex(s.coeffs[ordN].coeffs)
            for ind in eachindex(s.coeffs[ordN].coeffs[ihp].coeffs)
                c = s.coeffs[ordN].coeffs[ihp].coeffs[ind]
                vHP.coeffs[ihp] = c
                vT[ind] += TaylorN(vHP)
                vHP.coeffs[ihp] = zero(T)
            end
        end
    end
    Taylor1(vT)
end

function convert{T<:Number,N}(::Type{Array{TaylorN{Taylor1{T}},N}},
        s::Array{Taylor1{TaylorN{T}},N})
    v = Array{TaylorN{Taylor1{T}}}(size(s))
    for ind in eachindex(s)
        v[ind] = convert(TaylorN{Taylor1{T}}, s[ind])
    end
    v
end

function convert{T<:Number,N}(::Type{Array{Taylor1{TaylorN{T}},N}},
        s::Array{TaylorN{Taylor1{T}},N})
    v = Array{Taylor1{TaylorN{T}}}(size(s))
    for ind in eachindex(s)
        v[ind] = convert(Taylor1{TaylorN{T}}, s[ind])
    end
    v
end



# Promotion
promote_rule{T<:Number}(::Type{Taylor1{T}}, ::Type{Taylor1{T}}) = Taylor1{T}

promote_rule{T<:Number, S<:Number}(::Type{Taylor1{T}}, ::Type{Taylor1{S}}) =
    Taylor1{promote_type(T, S)}

promote_rule{T<:Number}(::Type{Taylor1{T}}, ::Type{Array{T,1}}) = Taylor1{T}

promote_rule{T<:Number, S<:Number}(::Type{Taylor1{T}}, ::Type{Array{S,1}}) =
    Taylor1{promote_type(T, S)}

promote_rule{T<:Number}(::Type{Taylor1{T}}, ::Type{T}) = Taylor1{T}

promote_rule{T<:Number, S<:Number}(::Type{Taylor1{T}}, ::Type{S}) =
    Taylor1{promote_type(T, S)}


promote_rule{T<:Number, S<:Number}(::Type{HomogeneousPolynomial{T}},
    ::Type{HomogeneousPolynomial{S}}) = HomogeneousPolynomial{promote_type(T,S)}

promote_rule{T<:Number, S<:Number}(::Type{HomogeneousPolynomial{T}},
    ::Type{Array{S,1}}) = HomogeneousPolynomial{promote_type(T, S)}

promote_rule{T<:Number,S<:Union{Real,Complex}}(::Type{HomogeneousPolynomial{T}},
    ::Type{S}) = HomogeneousPolynomial{promote_type(T,S)}


promote_rule{T<:Number, S<:Number}(::Type{TaylorN{T}}, ::Type{TaylorN{S}}) =
    TaylorN{promote_type(T, S)}

promote_rule{T<:Number, S<:Number}(::Type{TaylorN{T}},
    ::Type{HomogeneousPolynomial{S}}) = TaylorN{promote_type(T, S)}

promote_rule{T<:Number, S<:Number}(::Type{TaylorN{T}},
    ::Type{Array{HomogeneousPolynomial{S},1}}) = TaylorN{promote_type(T, S)}

promote_rule{T<:Number, S<:Number}(::Type{TaylorN{T}}, ::Type{S}) =
    TaylorN{promote_type(T, S)}
