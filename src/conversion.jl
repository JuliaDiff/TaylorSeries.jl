# This file is part of the Taylor1Series.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

## Conversion
convert(::Type{Taylor1{T}}, a::Taylor1) where {T<:Number} =
    Taylor1(convert(Array{T,1}, a.coeffs), a.order)

convert(::Type{Taylor1{T}}, a::Taylor1{T}) where {T<:Number} = a

convert(::Type{Taylor1{Rational{T}}}, a::Taylor1{S}) where
        {T<:Integer,S<:AbstractFloat} = Taylor1(rationalize.(a[:]), a.order)

convert(::Type{Taylor1{T}}, b::Array{T,1}) where {T<:Number} =
    Taylor1(b, length(b)-1)

convert(::Type{Taylor1{T}}, b::Array{S,1}) where {T<:Number,S<:Number} =
    Taylor1(convert(Array{T,1},b), length(b)-1)

convert(::Type{Taylor1{T}}, b::S)  where {T<:Number,S<:Number} =
    Taylor1([convert(T,b)], 0)

convert(::Type{Taylor1{T}}, b::T)  where {T<:Number} = Taylor1([b], 0)

convert(::Type{Taylor1}, a::T) where {T<:Number} = Taylor1(a, 0)


convert(::Type{HomogeneousPolynomial{T}}, a::HomogeneousPolynomial) where {T<:Number} =
    HomogeneousPolynomial(convert(Array{T,1}, a.coeffs), a.order)

convert(::Type{HomogeneousPolynomial{T}}, a::HomogeneousPolynomial{T}) where {T<:Number} =
    a

function convert(::Type{HomogeneousPolynomial{Rational{T}}},
        a::HomogeneousPolynomial{S}) where {T<:Integer,S<:AbstractFloat}

    la = length(a.coeffs)
    v = Array{Rational{T}}(undef, la)
    v .= rationalize.(a[1:la], tol=eps(one(S)))
    return HomogeneousPolynomial(v, a.order)
end

convert(::Type{HomogeneousPolynomial{T}}, b::Array{S,1}) where {T<:Number,S<:Number} =
    HomogeneousPolynomial(convert(Array{T,1}, b), orderH(b))

convert(::Type{HomogeneousPolynomial{T}}, b::S) where {T<:Number,S<:Number}=
    HomogeneousPolynomial([convert(T,b)], 0)

convert(::Type{HomogeneousPolynomial{T}}, b::Array{T,1}) where {T<:Number} =
    HomogeneousPolynomial(b, orderH(b))

convert(::Type{HomogeneousPolynomial{T}}, b::T) where {T<:Number} =
    HomogeneousPolynomial([b], 0)

convert(::Type{HomogeneousPolynomial}, a::Number) = HomogeneousPolynomial([a],0)


convert(::Type{TaylorN{T}}, a::TaylorN) where {T<:Number} =
    TaylorN( convert(Array{HomogeneousPolynomial{T},1}, a.coeffs), a.order)

convert(::Type{TaylorN{T}}, a::TaylorN{T}) where {T<:Number} = a

convert(::Type{TaylorN{T}}, b::HomogeneousPolynomial{S}) where {T<:Number,S<:Number} =
    TaylorN( [convert(HomogeneousPolynomial{T}, b)], b.order)

convert(::Type{TaylorN{T}}, b::Array{HomogeneousPolynomial{S},1}) where {T<:Number,S<:Number} =
    TaylorN( convert(Array{HomogeneousPolynomial{T},1}, b), maxorderH(b))

convert(::Type{TaylorN{T}}, b::S)  where {T<:Number,S<:Number} =
    TaylorN( [HomogeneousPolynomial([convert(T, b)], 0)], 0)

convert(::Type{TaylorN{T}}, b::HomogeneousPolynomial{T}) where {T<:Number} =
    TaylorN( [b], b.order)

convert(::Type{TaylorN{T}}, b::Array{HomogeneousPolynomial{T},1}) where {T<:Number} =
    TaylorN( b, maxorderH(b))

convert(::Type{TaylorN{T}}, b::T) where {T<:Number} =
    TaylorN( [HomogeneousPolynomial([b], 0)], 0)

convert(::Type{TaylorN}, b::Number) = TaylorN( [HomogeneousPolynomial([b], 0)], 0)


function convert(::Type{TaylorN{Taylor1{T}}}, s::Taylor1{TaylorN{T}}) where {T<:NumberNotSeries}
    orderN = maximum(get_order.(s[:]))
    r = zeros(HomogeneousPolynomial{Taylor1{T}}, orderN)

    v = zeros(T, s.order+1)
    @inbounds for ordT in eachindex(s)
        v[ordT+1] = one(T)
        @inbounds for ordHP in 0:s[ordT].order
            @inbounds for ic in eachindex(s[ordT][ordHP].coeffs)
                coef = s[ordT][ordHP][ic]
                r[ordHP+1][ic] += Taylor1( coef.*v )
            end
        end
        v[ordT+1] = zero(T)
    end
    return TaylorN(r)
end

function convert(::Type{Taylor1{TaylorN{T}}}, s::TaylorN{Taylor1{T}}) where {T<:NumberNotSeries}
    ordert = 0
    for ordHP in eachindex(s)
        ordert = max(ordert, s[ordHP][1].order)
    end
    vT = Array{TaylorN{T}}(undef, ordert+1)
    @inbounds for ordT in eachindex(vT)
        vT[ordT] = TaylorN(zero(T), s.order)
    end

    @inbounds for ordN in eachindex(s)
        vHP = HomogeneousPolynomial(zeros(T, length(s[ordN])))
        @inbounds for ihp in eachindex(s[ordN].coeffs)
            @inbounds for ind in eachindex(s[ordN][ihp].coeffs)
                c = s[ordN][ihp][ind-1]
                vHP[ihp] = c
                vT[ind] += TaylorN(vHP, s.order)
                vHP[ihp] = zero(T)
            end
        end
    end
    return Taylor1(vT)
end


# Promotion
promote_rule(::Type{Taylor1{T}}, ::Type{Taylor1{T}}) where {T<:Number} = Taylor1{T}

promote_rule(::Type{Taylor1{T}}, ::Type{Taylor1{S}}) where {T<:Number,S<:Number} =
    Taylor1{promote_type(T,S)}

promote_rule(::Type{Taylor1{T}}, ::Type{Array{T,1}}) where {T<:Number} = Taylor1{T}

promote_rule(::Type{Taylor1{T}}, ::Type{Array{S,1}}) where {T<:Number,S<:Number} =
    Taylor1{promote_type(T,S)}


promote_rule(::Type{Taylor1{T}}, ::Type{T}) where {T<:Number} = Taylor1{T}

promote_rule(::Type{Taylor1{T}}, ::Type{S}) where {T<:Number,S<:Number} =
    Taylor1{promote_type(T,S)}

promote_rule(::Type{Taylor1{Taylor1{T}}}, ::Type{Taylor1{T}}) where {T<:Number} =
    Taylor1{Taylor1{T}}


promote_rule(::Type{HomogeneousPolynomial{T}},
    ::Type{HomogeneousPolynomial{S}}) where {T<:Number,S<:Number} =
        HomogeneousPolynomial{promote_type(T,S)}

promote_rule(::Type{HomogeneousPolynomial{T}},
    ::Type{HomogeneousPolynomial{T}}) where {T<:Number} = HomogeneousPolynomial{T}

promote_rule(::Type{HomogeneousPolynomial{T}},
    ::Type{Array{S,1}}) where {T<:Number,S<:Number} = HomogeneousPolynomial{promote_type(T,S)}

promote_rule(::Type{HomogeneousPolynomial{T}}, ::Type{S}) where
    {T<:Number,S<:NumberNotSeries} = HomogeneousPolynomial{promote_type(T,S)}


promote_rule(::Type{TaylorN{T}}, ::Type{TaylorN{S}}) where {T<:Number,S<:Number} =
    TaylorN{promote_type(T,S)}

promote_rule(::Type{TaylorN{T}}, ::Type{HomogeneousPolynomial{S}}) where
    {T<:Number,S<:Number} = TaylorN{promote_type(T,S)}

promote_rule(::Type{TaylorN{T}}, ::Type{Array{HomogeneousPolynomial{S},1}}) where
    {T<:Number,S<:Number} = TaylorN{promote_type(T,S)}

promote_rule(::Type{TaylorN{T}}, ::Type{S}) where {T<:Number,S<:Number} =
    TaylorN{promote_type(T,S)}


promote_rule(::Type{S}, ::Type{T}) where
    {S<:AbstractIrrational,T<:AbstractSeries} = promote_rule(T,S)

promote_rule(::Type{Taylor1{T}}, ::Type{TaylorN{S}}) where {T<:NumberNotSeries,S<:NumberNotSeries} =
    throw(ArgumentError("There is no reasonable promotion among `Taylor1{$T}` and `TaylorN{$S}` types"))

promote_rule(::Type{Taylor1{TaylorN{T}}}, ::Type{TaylorN{Taylor1{S}}}) where
    {T<:NumberNotSeries, S<:NumberNotSeries} = Taylor1{TaylorN{promote_type(T,S)}}
promote_rule(::Type{TaylorN{Taylor1{T}}}, ::Type{Taylor1{TaylorN{S}}}) where
    {T<:NumberNotSeries, S<:NumberNotSeries} = Taylor1{TaylorN{promote_type(T,S)}}

# Different nested Taylor1's promotion methods; consider to 4 levels of nesting
let
    strucTT = :(Taylor1{T})
    strucTR = :(Taylor1{R})
    strucT = :(T)
    strucS = :(S)
    strucR = :(R)
    expr_elem = :(bnew[0])
    for j = 0:5
        for i = 1+j:6
            @eval function Base._promote(a::$strucTT, b::$strucS) where {T<:NumberNotSeries,
                    S<:NumberNotSeries}
                R = promote_type(T, S)
                anew = convert($strucTR, a)
                bnew = zero(anew)
                $(expr_elem) = convert($strucR, b)
                return anew, bnew
            end
            @eval Base._promote(b::$strucS, a::$strucTT) where {T<:NumberNotSeries,
                S<:NumberNotSeries} = reverse(Base._promote(a, b))
            strucTT = :(Taylor1{$strucTT})
            strucTR = :(Taylor1{$strucTR})
            expr_elem = :(($(expr_elem))[0])
        end
        strucS = :(Taylor1{$strucS})
        strucR = :(Taylor1{$strucR})
        strucT = :(Taylor1{$strucT})
        strucTT = deepcopy(strucT)
        strucTT = :(Taylor1{$strucT})
        strucTR = deepcopy(strucR)
        strucTR = :(Taylor1{$strucTR})
        expr_elem = :(bnew[0])
    end
end


# float
float(::Type{Taylor1{T}}) where T<:Number = Taylor1{float(T)}
float(::Type{HomogeneousPolynomial{T}}) where T<:Number = HomogeneousPolynomial{float(T)}
float(::Type{TaylorN{T}}) where T<:Number = TaylorN{float(T)}

float(x::Taylor1{T}) where T<:Number = convert(Taylor1{float(T)}, x)
float(x::HomogeneousPolynomial{T}) where T<:Number = convert(HomogeneousPolynomial{float(T)}, x)
float(x::TaylorN{T}) where T<:Number = convert(TaylorN{float(T)}, x)
