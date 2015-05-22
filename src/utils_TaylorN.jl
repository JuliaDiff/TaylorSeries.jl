# utils_TaylorN.jl: N-variables Taylor expansions through homogeneous polynomials
#
# Last modification: 2015.05.08
#
# Luis Benet & David P. Sanders
# UNAM
#


## Minimum order of a HomogeneousPolynomial compatible with the vector's length
function orderH{T}(coeffs::Array{T,1})
    ord = 0
    ll = length(coeffs)
    for i = 1:_params_taylorN.maxOrder+1
        @inbounds nCoefH = sizeTable[i]
        ll <= nCoefH && break
        ord += 1
    end
    return ord
end


## HomogeneousPolynomial (homogeneous polynomial) constructors ##
@doc """
DataType for *homogenous* polynomials in many (>1) independent variables

Fieldnames:

- `coeffs`: vector containing the expansion coefficients; the vector components
are related to the monomials by `indicesTables` and `posTable`

- `order` : order (degree) of the homogenous polynomial
""" ->
immutable HomogeneousPolynomial{T<:Number} <: Number
    coeffs  :: Array{T,1}
    order   :: Int

    function HomogeneousPolynomial( coeffs::Array{T,1}, order::Int )
        maxOrder = _params_taylorN.maxOrder
        @assert order <= maxOrder
        lencoef = length( coeffs )
        @inbounds nCoefH = sizeTable[order+1]
        @assert lencoef <= nCoefH
        nCoefH == lencoef && return new(coeffs, order)
        resize!(coeffs, nCoefH)
        @simd for i = lencoef+1:nCoefH
            @inbounds coeffs[i] = zero(T)
        end
        new(coeffs, order)
    end
end

HomogeneousPolynomial{T<:Number}(x::HomogeneousPolynomial{T}, order::Int) =
    HomogeneousPolynomial{T}(x.coeffs, order)
HomogeneousPolynomial{T<:Number}(x::HomogeneousPolynomial{T}) = x
HomogeneousPolynomial{T<:Number}(coeffs::Array{T,1}, order::Int) =
    HomogeneousPolynomial{T}(coeffs, order)
HomogeneousPolynomial{T<:Number}(coeffs::Array{T,1}) =
    HomogeneousPolynomial{T}(coeffs, orderH(coeffs))
HomogeneousPolynomial{T<:Number}(x::T, order::Int) =
    HomogeneousPolynomial{T}([x], order)
HomogeneousPolynomial{T<:Number}(x::T) = HomogeneousPolynomial{T}([x], 0)

eltype{T<:Number}(::HomogeneousPolynomial{T}) = T
length(a::HomogeneousPolynomial) = length( a.coeffs )
get_maxOrder(a::HomogeneousPolynomial) = a.order


## zero and one ##
function zero{T<:Number}(a::HomogeneousPolynomial{T})
    a.order == 0 && return HomogeneousPolynomial(zero(T), 0)
    @inbounds nCoefH = sizeTable[a.order+1]
    return HomogeneousPolynomial{T}( zeros(T,nCoefH), a.order)
end

function zeros{T<:Number}(::HomogeneousPolynomial{T}, order::Int)
    order == 0 && return [HomogeneousPolynomial(zero(T),0)]
    v = Array(HomogeneousPolynomial{T}, order+1)
    @simd for ord in eachindex(v)
        @inbounds nCoefH = sizeTable[ord]
        @inbounds v[ord] = HomogeneousPolynomial(zeros(T,nCoefH),ord-1)
    end
    return v
end

zeros{T<:Number}(::Type{HomogeneousPolynomial{T}}, order::Int) =
    zeros( HomogeneousPolynomial(zero(T), 0), order)

function one{T<:Number}(a::HomogeneousPolynomial{T})
    a.order == 0 && return HomogeneousPolynomial(one(T), 0)
    @inbounds nCoefH = sizeTable[a.order+1]
    return HomogeneousPolynomial{T}( ones(T,nCoefH), a.order)
end

function ones{T<:Number}(::HomogeneousPolynomial{T}, order::Int)
    order == 0 && return [HomogeneousPolynomial(one(T),0)]
    v = Array(HomogeneousPolynomial{T}, order+1)
    @simd for ord in eachindex(v)
        @inbounds nCoefH = sizeTable[ord]
        @inbounds v[ord] = HomogeneousPolynomial(ones(T,nCoefH),ord-1)
    end
    return v
end

ones{T<:Number}(::Type{HomogeneousPolynomial{T}}, order::Int) =
    ones( HomogeneousPolynomial(one(T), 0), order)

## Conversion and promotion rules ##
convert{T<:Number}(::Type{HomogeneousPolynomial{T}}, a::HomogeneousPolynomial) =
    HomogeneousPolynomial{T}(convert(Array{T,1}, a.coeffs), a.order)
convert{T<:Number, S<:Number}(::Type{HomogeneousPolynomial{T}}, b::Array{S,1}) =
    HomogeneousPolynomial{T}(convert(Array{T,1}, b), orderH(b))
convert{T<:Number}(::Type{HomogeneousPolynomial{T}}, b::Number) =
    HomogeneousPolynomial{T}([convert(T,b)], 0)
convert{T<:Number}(::Type{HomogeneousPolynomial{T}},
    a::HomogeneousPolynomial{T}) = a
convert{T<:Number}(::Type{HomogeneousPolynomial{T}}, b::Array{T,1}) =
    HomogeneousPolynomial{T}(b, orderH(b))
convert{T<:Number}(::Type{HomogeneousPolynomial{T}}, b::T) =
    HomogeneousPolynomial{T}([b], 0)

promote_rule{T<:Number, S<:Number}(::Type{HomogeneousPolynomial{T}},
    ::Type{HomogeneousPolynomial{S}}) = HomogeneousPolynomial{promote_type(T,S)}
promote_rule{T<:Number, S<:Number}(::Type{HomogeneousPolynomial{T}},
    ::Type{Array{S,1}}) = HomogeneousPolynomial{promote_type(T, S)}
promote_rule{T<:Number,S<:Union(Real,Complex)}(::Type{HomogeneousPolynomial{T}},
    ::Type{S}) = HomogeneousPolynomial{promote_type(T,S)}

## Maximum order of a HomogeneousPolynomial vector; used by TaylorN constructor
function maxorderH{T<:Number}(v::Array{HomogeneousPolynomial{T},1})
    m = 0
    @inbounds for i in eachindex(v)
        m = max(m, v[i].order)
    end
    return m
end



@doc """
DataType for polynomial expansions in many (>1) independent variables

Fieldnames:

- `coeffs`: vector containing the `HomogeneousPolynomial` entries

- `order` : maximum order of the polynomial expansion

""" ->
immutable TaylorN{T<:Number} <: Number
    coeffs  :: Array{HomogeneousPolynomial{T},1}
    order   :: Int

    function TaylorN( v::Array{HomogeneousPolynomial{T},1}, order::Int )
        m = maxorderH(v)
        order = max( m, order )
        coeffs = zeros(HomogeneousPolynomial{T}, order)
        @simd for i in eachindex(v)
            @inbounds ord = v[i].order
            @inbounds coeffs[ord+1] += v[i]
        end
        new(coeffs, order)
    end
end

TaylorN{T<:Number}(x::TaylorN{T}, order::Int) = TaylorN{T}(x.coeffs, order)
TaylorN{T<:Number}(x::TaylorN{T}) = x
TaylorN{T<:Number}(x::Array{HomogeneousPolynomial{T},1}, order::Int) =
    TaylorN{T}(x, order)
TaylorN{T<:Number}(x::Array{HomogeneousPolynomial{T},1}) = TaylorN{T}(x,0)
TaylorN{T<:Number}(x::HomogeneousPolynomial{T},order::Int) =
    TaylorN{T}([x], order)
TaylorN{T<:Number}(x::HomogeneousPolynomial{T}) = TaylorN{T}([x], x.order)
TaylorN{T<:Number}(x::T,order::Int) =
    TaylorN{T}([HomogeneousPolynomial(x)],order)
TaylorN{T<:Number}(x::T) = TaylorN{T}([HomogeneousPolynomial(x)], 0)

## Shortcut to define TaylorN independent variables
function taylorN_variable(T::Type, nv::Int, order::Int=_params_taylorN.maxOrder)
    @assert 0 < nv <= _params_taylorN.numVars
    v = zeros(T, _params_taylorN.numVars)
    @inbounds v[nv] = one(T)
    return TaylorN( HomogeneousPolynomial(v,1), order )
end
taylorN_variable(nv::Int, order::Int=_params_taylorN.maxOrder) =
    taylorN_variable(Float64, nv, order)

## get_coeff
function get_coeff(a::HomogeneousPolynomial, v::Array{Int,1})
    @assert length(v) == _params_taylorN.numVars
    kdic = hash(v)
    @inbounds n = posTable[a.order+1][kdic]
    a.coeffs[n]
end
function get_coeff(a::TaylorN, v::Array{Int,1})
    order = sum(v)
    get_coeff(a.coeffs[order+1], v)
end

## Type, length ##
eltype{T<:Number}(::TaylorN{T}) = T
length(a::TaylorN) = length( a.coeffs )
get_maxOrder(x::TaylorN) = x.order

## zero and one ##
zero{T<:Number}(a::TaylorN{T}) = TaylorN(zero(T), a.order)
one{T<:Number}(a::TaylorN{T}) = TaylorN(one(T), a.order)

## Conversion and promotion rules ##
convert{T<:Number}(::Type{TaylorN{T}}, a::TaylorN) =
    TaylorN{T}( convert(Array{HomogeneousPolynomial{T},1}, a.coeffs), a.order)
convert{T<:Number, S<:Number}(::Type{TaylorN{T}}, b::HomogeneousPolynomial{S}) =
    TaylorN{T}( [convert(HomogeneousPolynomial{T}, b)], b.order)
convert{T<:Number, S<:Number}(::Type{TaylorN{T}},
    b::Array{HomogeneousPolynomial{S},1}) =
    TaylorN{T}( convert(Array{HomogeneousPolynomial{T},1}, b), length(b)-1)
convert{T<:Number}(::Type{TaylorN{T}}, b::Number) =
    TaylorN( [HomogeneousPolynomial(convert(T, b))], 0)
convert{T<:Number}(::Type{TaylorN{T}}, a::TaylorN{T}) = a
convert{T<:Number}(::Type{TaylorN{T}}, b::HomogeneousPolynomial{T}) =
    TaylorN{T}( [b], b.order)
convert{T<:Number}(::Type{TaylorN{T}}, b::Array{HomogeneousPolynomial{T},1}) =
    TaylorN{T}( b, length(b)-1)
convert{T<:Number}(::Type{TaylorN{T}}, b::T) =
    TaylorN( [HomogeneousPolynomial(b)], 0)

promote_rule{T<:Number, S<:Number}(::Type{TaylorN{T}}, ::Type{TaylorN{S}}) =
    TaylorN{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{TaylorN{T}},
    ::Type{HomogeneousPolynomial{S}}) = TaylorN{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{TaylorN{T}},
    ::Type{Array{HomogeneousPolynomial{S},1}}) = TaylorN{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{TaylorN{T}}, ::Type{S}) =
    TaylorN{promote_type(T, S)}

## Auxiliary function ##
fixorder{T<:Number}(a::TaylorN{T}, order::Int64) = TaylorN{T}(a.coeffs, order)

function fixorder(a::TaylorN, b::TaylorN)
    a.order == b.order && return a, b
    a.order < b.order && return TaylorN(a.coeffs, b.order), b
    return a, TaylorN(b.coeffs, a.order)
end

function fixshape(a::HomogeneousPolynomial, b::HomogeneousPolynomial)
    @assert a.order == b.order
    return promote(a, b)
end

function fixshape(a::TaylorN, b::TaylorN)
    if eltype(a) != eltype(b)
        a, b = promote(a, b)
    end
    return fixorder(a, b)
end

## real, imag, conj and ctranspose ##
for TT in (:HomogeneousPolynomial, :TaylorN), f in (:real, :imag, :conj)
    @eval ($f){T<:Number}(a::($TT){T}) = ($TT)(($f)(a.coeffs), a.order)
end

ctranspose{T<:Number}(a::HomogeneousPolynomial{T}) = conj(a)
ctranspose{T<:Number}(a::TaylorN{T}) = conj(a)


## Equality ##
==(a::HomogeneousPolynomial, b::HomogeneousPolynomial) = a.coeffs == b.coeffs
function ==(a::TaylorN, b::TaylorN)
    a, b = fixshape(a, b)
    test = true
    for i in eachindex(a.coeffs)
        @inbounds test = a.coeffs[i] == b.coeffs[i]
        test || return test
    end
    return test
end

function iszero{T}(a::HomogeneousPolynomial{T})
    test = true
    for i in eachindex(a.coeffs)
        @inbounds test = a.coeffs[i] == zero(T)
        test || return test
    end
    return test
end

## Addition and substraction ##
for T in (:HomogeneousPolynomial, :TaylorN), f in (:+, :-)
    @eval begin
        function ($f)(a::($T), b::($T))
            a, b = fixshape(a, b)
            v = Array(eltype(a.coeffs), length(a.coeffs))
            @simd for i in eachindex(v)
                @inbounds v[i] = ($f)(a.coeffs[i], b.coeffs[i])
            end
            return ($T)(v, a.order)
        end
        function ($f)(a::($T))
            v = Array(eltype(a.coeffs), length(a.coeffs))
            @simd for i in eachindex(v)
                @inbounds v[i] = ($f)(a.coeffs[i])
            end
            return ($T)(v, a.order)
        end
    end
end

## Multiplication ##
function *(a::HomogeneousPolynomial, b::HomogeneousPolynomial)
    T = promote_type( eltype(a), eltype(b) )
    order = a.order + b.order
    if order > _params_taylorN.maxOrder
        return HomogeneousPolynomial(zero(T), _params_taylorN.maxOrder)
    end
    (iszero(a) || iszero(b)) && return HomogeneousPolynomial(zero(T), order)

    @inbounds nCoefHa = sizeTable[a.order+1]
    @inbounds nCoefHb = sizeTable[b.order+1]
    @inbounds nCoefH  = sizeTable[order+1]
    if eltype(a) != eltype(b)
        a, b = promote(a, b)
    end

    coeffs = zeros(T, nCoefH)
    iaux = zeros(Int, _params_taylorN.numVars)
    @inbounds posTb = posTable[order+1]
    @inbounds for na = 1:nCoefHa
        ca = a.coeffs[na]
        ca == zero(T) && continue
        inda = indicesTable[a.order+1][na]
        @inbounds for nb = 1:nCoefHb
            cb = b.coeffs[nb]
            cb == zero(T) && continue
            indb = indicesTable[b.order+1][nb]
            @simd for i = 1:_params_taylorN.numVars
                @inbounds iaux[i] = inda[i]+indb[i]
            end
            kdic = hash(iaux)
            pos = posTb[kdic]
            coeffs[pos] += ca * cb
        end
    end

    return HomogeneousPolynomial{T}(coeffs, order)
end

function *(a::TaylorN, b::TaylorN)
    a, b = fixshape(a, b)
    T = eltype(a)
    coeffs = zeros(HomogeneousPolynomial{T}, a.order)

    for ord in eachindex(coeffs)
        @inbounds for i = 0:ord-1
            (iszero(a.coeffs[i+1]) || iszero(b.coeffs[ord-i])) && continue
            coeffs[ord] += a.coeffs[i+1] * b.coeffs[ord-i]
        end
    end

    return TaylorN{T}(coeffs, a.order)
end

## Division ##
/(a::HomogeneousPolynomial, x::Real) = a*inv(x)
/(a::HomogeneousPolynomial, x::Complex) = a*inv(x)
function /(a::TaylorN, b::TaylorN)
    @inbounds b0 = b.coeffs[1].coeffs[1]
    @assert b0 != zero(b0)
    a, b = fixshape(a, b)
    # order and coefficient of first factorized term
    # orddivfact, cdivfact = divfactorization(a, b)
    b0 = inv(b0)
    @inbounds cdivfact = a.coeffs[1] * b0
    T = eltype(cdivfact)
    coeffs = zeros(HomogeneousPolynomial{T}, a.order)
    @inbounds coeffs[1] = cdivfact

    for ord in eachindex(coeffs)
        ord == a.order+1 && continue
        @inbounds for i = 0:ord-1
            (iszero(coeffs[i+1]) || iszero(b.coeffs[ord-i+1])) && continue
            coeffs[ord+1] += coeffs[i+1] * b.coeffs[ord-i+1]
        end
        @inbounds coeffs[ord+1] = (a.coeffs[ord+1] - coeffs[ord+1]) * b0
    end

    return TaylorN{T}(coeffs, a.order)
end

## TODO: Implement factorization (divfactorization) for TaylorN polynomials

## Division functions: rem and mod
for op in (:mod, :rem)
    @eval begin
        @inbounds function ($op){T<:Real,S<:Real}(a::TaylorN{T}, x::S)
            coeffs = copy(a.coeffs)
            y = ($op)(a.coeffs[1].coeffs[1], x)
            coeffs[1] = HomogeneousPolynomial(y)
            return TaylorN{T}( coeffs, a.order )
        end
    end
end

function mod2pi{T<:Real}(a::TaylorN{T})
    coeffs = copy(a.coeffs)
    @inbounds y = mod2pi(a.coeffs[1].coeffs[1])
    @inbounds coeffs[1] = HomogeneousPolynomial(y)
    return TaylorN{T}( coeffs, a.order )
end

## Int power ##
function ^(a::HomogeneousPolynomial, n::Integer)
    n == 0 && return one(a)
    n == 1 && return a
    n == 2 && return square(a)
    n < 0 && throw(DomainError())
    return power_by_squaring(a, n)
end

function ^{T<:Number}(a::TaylorN{T}, n::Integer)
    n == 0 && return one(a)
    n == 1 && return a
    n == 2 && return square(a)
    n < 0 && return inv( a^(-n) )
    return power_by_squaring(a, n)
end

function ^{T<:Integer}(a::TaylorN{T}, n::Integer)
    n == 0 && return one(a)
    n == 1 && return a
    n == 2 && return square(a)
    n < 0 && throw(DomainError())
    return power_by_squaring(a, n)
end

## power_by_squaring; slightly modified from base/intfuncs.jl
## Licensed under MIT "Expat"
for T in (:HomogeneousPolynomial, :TaylorN)
    @eval begin
        function power_by_squaring(x::($T), p::Integer)
            p == 1 && return x# copy(x)
            p == 0 && return one(x)
            p == 2 && return square(x)
            t = trailing_zeros(p) + 1
            p >>= t

            while (t -= 1) > 0
                x = square(x)
            end

            y = x
            while p > 0
                t = trailing_zeros(p) + 1
                p >>= t
                while (t -= 1) >= 0
                    x = square(x)
                end
                y *= x
            end

            return y
        end
    end
end

## Rational power ##
^(a::TaylorN, x::Rational) = a^(x.num/x.den)

^(a::TaylorN, b::TaylorN) = exp( b*log(a) )

## Real power ##
function ^{S<:Real}(a::TaylorN, x::S)
    x == zero(x) && return TaylorN( one(eltype(a)) )
    x == one(x)/2 && return sqrt(a)
    isinteger(x) && return a^round(Int,x)
    @inbounds a0 = a.coeffs[1].coeffs[1]
    @assert a0 != zero(a0)
    aux = ( a0 )^x
    T = typeof(aux)
    coeffs = zeros(HomogeneousPolynomial{T}, a.order)
    @inbounds coeffs[1] = HomogeneousPolynomial( aux )

    for ord in eachindex(coeffs)
        ord == a.order+1 && continue
        @inbounds for i = 0:ord-1
            tt = x*(ord-i)-i
            cpol = coeffs[i+1]
            apol = a.coeffs[ord-i+1]
            (iszero(cpol) || iszero(apol)) && continue
            coeffs[ord+1] += tt * cpol * apol
        end
        coeffs[ord+1] = coeffs[ord+1] / (ord*a0)
    end

    return TaylorN{T}(coeffs, a.order)
end
^(a::TaylorN, x::Complex) = exp( x*log(a) )

## Square ##
function square(a::HomogeneousPolynomial)
    T = eltype(a)
    order = 2*a.order
    if order > _params_taylorN.maxOrder
        return HomogeneousPolynomial(zero(T), _params_taylorN.maxOrder)
    end
    @inbounds nCoefHa = sizeTable[a.order+1]
    @inbounds nCoefH  = sizeTable[order+1]
    two = convert(T,2)
    coeffs = zeros(T, nCoefH)
    iaux = zeros( _params_taylorN.numVars )
    @inbounds posTb = posTable[order+1]

    @inbounds for na = 1:nCoefHa
        ca = a.coeffs[na]
        ca == zero(T) && continue
        inda = indicesTable[a.order+1][na]
        @inbounds for i = 1:_params_taylorN.numVars
            iaux[i] = 2inda[i]
        end
        kdic = hash(iaux)
        pos = posTb[kdic]
        coeffs[pos] += ca * ca
        @inbounds for nb = na+1:nCoefHa
            cb = a.coeffs[nb]
            cb == zero(T) && continue
            indb = indicesTable[a.order+1][nb]
            @simd for i = 1:_params_taylorN.numVars
                @inbounds iaux[i] = inda[i]+indb[i]
            end
            kdic = hash(iaux)
            pos = posTb[kdic]
            coeffs[pos] += two * ca * cb
        end
    end

    return HomogeneousPolynomial{T}(coeffs, order)
end

function square(a::TaylorN)
    T = eltype(a)
    coeffs = zeros(HomogeneousPolynomial{T}, a.order)
    @inbounds coeffs[1] = square(a.coeffs[1])

    two = convert(T,2)
    for ord in eachindex(coeffs)
        ord == a.order+1 && continue
        kodd = ord%2
        kord = div(ord-2+kodd, 2)
        @inbounds for i = 0 : kord
            (iszero(a.coeffs[i+1]) || iszero(a.coeffs[ord-i+1])) && continue
            coeffs[ord+1] += a.coeffs[i+1] * a.coeffs[ord-i+1]
        end
        @inbounds coeffs[ord+1] = two * coeffs[ord+1]
        kodd == 1 && continue
        kodd = div(ord,2)
        @inbounds coeffs[ord+1] += square( a.coeffs[kodd+1] )
    end

    return TaylorN{T}(coeffs, a.order)
end

## sqrt ##
function sqrt(a::TaylorN)
    @inbounds p0 = sqrt( a.coeffs[1].coeffs[1] )
    T = typeof(p0)
    coeffs = zeros(HomogeneousPolynomial{T}, a.order)
    @inbounds coeffs[1] = HomogeneousPolynomial( p0 )
    two = convert(T,2)

    for ord in eachindex(coeffs)
        ord == a.order+1 && continue
        kodd = ord%2
        kord = div(ord-2+kodd, 2)
        @inbounds for i = 1:kord
            coeffs[ord+1] += coeffs[i+1] * coeffs[ord-i+1]
        end
        @inbounds coeffs[ord+1] = a.coeffs[ord+1] - two * coeffs[ord+1]
        if iseven(ord)
            @inbounds coeffs[ord+1]=coeffs[ord+1]-square( coeffs[div(ord,2)+1])
        end
        @inbounds coeffs[ord+1] = coeffs[ord+1] / (two * p0)
    end

    return TaylorN{T}(coeffs, a.order)
end

## exp ##
function exp(a::TaylorN)
    order = a.order
    @inbounds aux = exp( a.coeffs[1].coeffs[1] )
    T = typeof(aux)
    coeffs = zeros(HomogeneousPolynomial{T}, order)
    @inbounds coeffs[1] = HomogeneousPolynomial(aux, 0)

    @inbounds for ord in eachindex(coeffs)
        ord == order+1 && continue
        @inbounds for j = 0:ord-1
            coeffs[ord+1] += (ord-j) * a.coeffs[ord-j+1] * coeffs[j+1]
        end
        @inbounds coeffs[ord+1] = coeffs[ord+1] / ord
    end

    return TaylorN{T}(coeffs, order)
end

## log ##
function log(a::TaylorN)
    order = a.order
    @inbounds a0 = a.coeffs[1].coeffs[1]
    l0 = log( a0 )
    T = typeof(l0)
    coeffs = zeros(HomogeneousPolynomial{T}, order)
    @inbounds coeffs[1] = HomogeneousPolynomial(l0)

    @inbounds for ord in eachindex(coeffs)
        ord == order+1 && continue
        @inbounds for j = 1:ord-1
            coeffs[ord+1] += j * a.coeffs[ord-j+1] * coeffs[j+1]
        end
        coeffs[ord+1] = (a.coeffs[ord+1] - coeffs[ord+1] / ord ) / a0
    end

    return TaylorN{T}(coeffs, order)
end

## sin and cos ##
sin(a::TaylorN) = imag( exp(im*a) )
cos(a::TaylorN) = real( exp(im*a) )
# sin(a::TaylorN) = sincos(a)[1]
# cos(a::TaylorN) = sincos(a)[2]
# function sincos(a::TaylorN)
#     order = a.order
#     @inbounds a0 = a.coeffs[1].coeffs[1]
#     s0 = sin( a0 )
#     c0 = cos( a0 )
#     T = typeof(s0)
#     coeffsSin = zeros(HomogeneousPolynomial{T}, order)
#     coeffsCos = zeros(HomogeneousPolynomial{T}, order)
#     @inbounds coeffsSin[1] = HomogeneousPolynomial(s0)
#     @inbounds coeffsCos[1] = HomogeneousPolynomial(c0)
#
#     @inbounds for ord in eachindex(coeffsSin)
#         ord == order+1 && continue
#         @inbounds for j = 0:ord-1
#             coeffsSin[ord+1] += (ord-j) * a.coeffs[ord-j+1] * coeffsCos[j+1]
#             coeffsCos[ord+1] += (ord-j) * a.coeffs[ord-j+1] * coeffsSin[j+1]
#         end
#         @inbounds coeffsSin[ord+1] =  coeffsSin[ord+1] / ord
#         @inbounds coeffsCos[ord+1] = -coeffsCos[ord+1] / ord
#     end
#     return TaylorN{T}(coeffsSin, order), TaylorN{T}(coeffsCos, order)
# end

## tan ##
function tan(a::TaylorN)
    order = a.order
    @inbounds a0 = a.coeffs[1].coeffs[1]
    t0 = tan(a0)
    T = typeof(t0)
    coeffsTan = zeros(HomogeneousPolynomial{T}, order)
    @inbounds coeffsTan[1] = HomogeneousPolynomial(t0)

    @inbounds for ord in eachindex(coeffsTan)
        ord == order+1 && continue
        v = coeffsTan[1:ord]
        tAux = (TaylorN(v, ord))^2
        @inbounds for j = 0:ord-1
            coeffsTan[ord+1] += (ord-j) * a.coeffs[ord-j+1] * tAux.coeffs[j+1]
        end
        coeffsTan[ord+1] = a.coeffs[ord+1] + coeffsTan[ord+1] / ord
    end

    return TaylorN{T}(coeffsTan, order)
end

## Differentiation ##
"""Partial differentiation of a HomogeneousPolynomial series with respect
to the r-th variable"""
function diffTaylor(a::HomogeneousPolynomial, r::Int)
    @assert 1 <= r <= _params_taylorN.numVars
    T = eltype(a)
    a.order == 0 && return HomogeneousPolynomial(zero(T))
    @inbounds nCoefH = sizeTable[a.order]
    coeffs = zeros(T,nCoefH)
    @inbounds posTb = posTable[a.order]
    @inbounds nCoefH = sizeTable[a.order+1]

    @inbounds for i = 1:nCoefH
        iind = indicesTable[a.order+1][i]
        n = iind[r]
        n == 0 && continue
        iind[r] -= 1
        kdic = hash(iind)
        pos = posTb[kdic]
        coeffs[pos] = n * a.coeffs[i]
        iind[r] += 1
    end

    return HomogeneousPolynomial{T}(coeffs, a.order-1)
end

"""Partial differentiation of a TaylorN series with respect
to the r-th variable"""
function diffTaylor(a::TaylorN, r::Int)
    T = eltype(a)
    coeffs = Array(HomogeneousPolynomial{T},a.order)

    @inbounds for ord in eachindex(coeffs)
        ord == a.order+1 && continue
        coeffs[ord] = diffTaylor( a.coeffs[ord+1], r)
    end
    return TaylorN{T}( coeffs, a.order )
end

diffTaylor(a::TaylorN) = diffTaylor(a, 1)

## Gradient, jacobian and hessian
function gradient(f::TaylorN)
    T = eltype(f)
    numVars = get_numVars()
    grad = Array(TaylorN{T}, numVars)
    @inbounds for nv = 1:numVars
        grad[nv] = diffTaylor(f, nv)
    end
    return grad
end

∇(f::TaylorN) = gradient(f)

function jacobian{T<:Number}(vf::Array{TaylorN{T},1})
    numVars = get_numVars()
    @assert length(vf) == numVars
    jac = Array(T, (numVars,numVars))

    @inbounds for comp = 1:numVars
        jac[:,comp] = vf[comp].coeffs[2].coeffs[1:end]
    end

    return transpose(jac)
end

function jacobian{T<:Number,S<:Number}(vf::Array{TaylorN{T},1},vals::Array{S,1})
    R = promote_type(T,S)
    numVars = get_numVars()
    @assert length(vf) == numVars == length(vals)
    jac = Array(R, (numVars,numVars))

    for comp = 1:numVars
        @inbounds grad = gradient( vf[comp] )
        @inbounds for nv = 1:numVars
            jac[nv,comp] = evalTaylor(grad[nv], vals)
        end
    end

    return transpose(jac)
end

hessian{T<:Number,S<:Number}(f::TaylorN{T}, vals::Array{S,1}) =
    (R = promote_type(T,S); jacobian( gradient(f), vals::Array{R,1}) )
hessian{T<:Number}(f::TaylorN{T}) = hessian( f, zeros(T, get_numVars()) )

## TODO: Integration...

## Evaluates a Taylor polynomial on a given point ##
# NEEDS REVISION since results are not quite precise
function evalHomog{T<:Number,S<:Number}(a::HomogeneousPolynomial{T},
        vals::Array{S,1} )

    numVars = get_numVars()
    @assert length(vals) == numVars
    R = promote_type(T,S)
    suma = zero(R)
    order = a.order
    @inbounds indTb = indicesTable[order+1]
    nCoefH = sizeTable[order+1]

    for pos = nCoefH:-1:1
        @inbounds iIndices = indTb[pos]
        @inbounds c = a.coeffs[pos]
        c == zero(R) && continue
        @inbounds for k = 1:numVars
            c *= vals[k]^iIndices[k]
        end
        suma += c
    end

    return suma
end

function evalTaylor{T<:Number,S<:Number}(a::TaylorN{T}, vals::Array{S,1} )
    @assert length(vals) == get_numVars()
    R = promote_type(T,S)
    suma = zero(R)

    for ord = a.order:-1:0
        @inbounds polH = a.coeffs[ord+1]
        suma += evalHomog( polH, vals )
    end

    return suma
end

evalTaylor{T<:Number}(a::TaylorN{T}) = evalTaylor(a, zeros(T, get_numVars()))
