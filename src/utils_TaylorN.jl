# utils_TaylorN.jl: N-variables Taylor expansions through homogeneous polynomials
#
# Last modification: 2014.06.22
#
# Luis Benet & David P. Sanders
# UNAM
#


## Default values for the maximum degree of polynomials (MAXORDER <--> nPart) and 
##  number of variables considered (NUMVARS <--> lNiv)
const MAXORDER = [6]
const NUMVARS  = [2]
info( string("\n MAXORDER = ", MAXORDER[end], "\n NUMVARS  = ",NUMVARS[end]) )

## Hash tables
#=
  `generateIndicesTable`: generates the array of dictionaries `indicesTable`. Then, `indicesTable[k+1]`
  contains the dictionary with the table defining the indices from the position for the homogeneous
  polynomial of degree k. 
  The number of coefficients of the homogenous polynomial of degree k is stored in `sizeTable[k+1]`.
  `generatePosTable`: generates the dictionary `posTable`, the "inverse" of `indicesTable`
=#
function generateIndicesTable()
    numVars = NUMVARS[end]
    maxOrd = MAXORDER[end]
    arrayDDic = Dict{Int,Array{Int,1}}[]
    arraySize = Int[]
    # if numVars==1
    #     for k = 0:maxOrd
    #         DDic = Dict{Int, Array{Int,1}}()
    #         DDic[k+1] = [k]
    #         push!(arraySize, 1)
    #         push!(arrayDDic, DDic)
    #     end
    #     return arrayDDic, arraySize
    # end
    for kDeg = 0:maxOrd
        nCoefH = binomial( numVars+kDeg-1, kDeg)
        iindices = zeros(Int, numVars)
        pos = 0
        DDic = Dict{Int, Array{Int,1}}()
        iV = numVars
        for iz = 0:kDeg
            @inbounds iindices[end] = iz
            @inbounds pos, DDic[pos] = pos2indices!( iV, kDeg, iindices, pos, DDic)
        end
        ##@assert length(DDic) == nCoefH
        push!(arrayDDic, DDic)
        push!(arraySize, nCoefH)
    end
    return arrayDDic, arraySize
end
function pos2indices!(iV::Int, kDeg::Int, iIndices::Array{Int,1}, pos::Int, dict::Dict{Int,Array{Int,1}})
    #
    jVar = iV-1
    kDegNext = kDeg - iIndices[iV]
    if jVar > 1
        for jDeg=0:kDegNext
            iIndices[jVar] = jDeg
            @inbounds pos, dict[pos] = pos2indices!( jVar, kDegNext, iIndices, pos, dict)
        end
    else
        iIndices[1] = kDegNext
        pos += 1
        @inbounds dict[pos] = iIndices[1:end]
    end
    return pos, iIndices[1:end]
end

const indicesTable, sizeTable = generateIndicesTable()

function generatePosTable()
    maxOrd  = MAXORDER[end]
    arrayDDic = Dict{Array{Int,1},Int}[]
    # if NUMVARS[end]==1
    #     for k = 0:maxOrd
    #         DDic = Dict{Array{Int,1}, Int}()
    #         DDic[ [k] ] = k+1
    #         push!(arrayDDic, DDic)
    #     end
    #     return arrayDDic
    # end
    for kDeg = 0:maxOrd
        DDic = Dict{Array{Int,1}, Int}()
        @inbounds nCoefH = sizeTable[kDeg+1]
        for i = 1:nCoefH
            @inbounds v = indicesTable[kDeg+1][i]
            DDic[ v ] = i
        end
        push!(arrayDDic,DDic)
    end
    return arrayDDic
end
const posTable = generatePosTable()

## Utilities to get/set MAXORDER and NUMVARS; they reset the indicesTable
get_maxOrder() = MAXORDER[end]
function set_maxOrder(n::Int)
    ##@assert n > 0
    oldOrder = MAXORDER[end]
    MAXORDER[end] = n
    if n < oldOrder
        for i = oldOrder:-1:n+1
            pop!(indicesTable);
            pop!(sizeTable);
            pop!(posTable);
        end
    else
        resize!(indicesTable,n+1)
        resize!(sizeTable,n+1)
        resize!(posTable,n+1)
        indicesTable[:], sizeTable[:] = generateIndicesTable()
        posTable[:] = generatePosTable()
    end
    info(string("MAXORDER is now ", n, "; hash tables regenerated.\n"))
    return n
end
#
get_numVars() = NUMVARS[end]
function set_numVars(n::Int)
    #@assert n > 0
    n==1 && error( string("Use `Taylor` rather than `TaylorN` for one independent variable.") )
    NUMVARS[end] = n
    maxOrd = MAXORDER[end]
    indicesTable[:], sizeTable[:] = generateIndicesTable()
    posTable[:] = generatePosTable()
    info(string("NUMVARS is now ", n, "; hash tables regenerated.\n"))
    return n
end

## Obtains the minimum order of a HomogPol compatible with the length of the vector
function orderH{T}(coeffs::Array{T,1})
    ord = 0
    ll = length(coeffs)
    for i=1:MAXORDER[end]+1
        @inbounds nCoefH = sizeTable[i]
        ll <= nCoefH && break
        ord += 1
    end
    return ord
end

## HomogPol (homogeneous polynomial) constructors ##
immutable HomogPol{T<:Number} <: AbstractSeries{T,NUMVARS[end]}
    coeffs  :: Array{T,1}
    order   :: Int
    function HomogPol(coeffs::Array{T,1}, order::Int)
        #@assert order <= MAXORDER[end]
        lencoef = length( coeffs )
        @inbounds nCoefH = sizeTable[order+1]
        #@assert lencoef <= nCoefH
        v = zeros(T, nCoefH)
        @inbounds v[1:lencoef] = coeffs[1:lencoef]
        new(v, order)
    end
end
HomogPol{T<:Number}(x::HomogPol{T}, order::Int) = HomogPol{T}(x.coeffs, order)
HomogPol{T<:Number}(x::HomogPol{T}) = HomogPol{T}(x.coeffs, x.order)
HomogPol{T<:Number}(coeffs::Array{T,1}, order::Int) = HomogPol{T}(coeffs, order)
HomogPol{T<:Number}(coeffs::Array{T,1}) = HomogPol{T}(coeffs, orderH(coeffs))
HomogPol{T<:Number}(x::T, order::Int) = HomogPol{T}([x], order)
HomogPol{T<:Number}(x::T) = HomogPol{T}([x], 0)

eltype{T<:Number}(::HomogPol{T}) = T
length(a::HomogPol) = length( a.coeffs )
get_numVars(x::HomogPol) = NUMVARS[end]
get_maxOrder(a::HomogPol) = a.order

## zero and one ##
function zero{T<:Number}(a::HomogPol{T})
    a.order == 0 && return HomogPol(zero(T), 0)
    @inbounds nCoefH = sizeTable[a.order+1]
    HomogPol(zeros(T,nCoefH), a.order)
end
function zeros{T<:Number}(a::HomogPol{T}, maxOrd::Int)
    order = max( maxOrd, a.order)
    #@assert 0 <= maxOrd <= MAXORDER[end]
    v = HomogPol{T}[]
    for ord = 0:order
        z = HomogPol(zero(T),ord)
        push!(v,z)
    end
    return v
end
zeros{T<:Number}(::Type{HomogPol{T}}, maxOrd::Int) = zeros( HomogPol(zero(T)), maxOrd)
function one{T<:Number}(a::HomogPol{T})
    a.order == 0 && HomogPol(one(T), 0)
    nCoefH = sizeTable[a.order+1]
    HomogPol(ones(T,nCoefH), a.order)
end
function ones{T<:Number}(::HomogPol{T}, maxOrd::Int)
    #@assert 0 <= maxOrd <= MAXORDER[end]
    v = HomogPol{T}[]
    for ord = 0:maxOrd
        z = HomogPol(one(T),ord)
        push!(v,z)
    end
    return v
end
ones{T<:Number}(::Type{HomogPol{T}}, maxOrd::Int) = ones( HomogPol(one(T)), maxOrd)

## Conversion and promotion rules ##
convert{T<:Number}(::Type{HomogPol{T}}, a::HomogPol) = 
    HomogPol(convert(Array{T,1}, a.coeffs), a.order)
convert{T<:Number, S<:Number}(::Type{HomogPol{T}}, b::Array{S,1}) = 
    HomogPol(convert(Array{T,1},b))
convert{T<:Number}(::Type{HomogPol{T}}, b::Number) = 
    HomogPol([convert(T,b)], 0)
#
promote_rule{T<:Number, S<:Number}(::Type{HomogPol{T}}, ::Type{HomogPol{S}}) = 
    HomogPol{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{HomogPol{T}}, ::Type{Array{S,1}}) = 
    HomogPol{promote_type(T, S)}
# Defined this way to permit promotion of HomogPol to TaylorN; see below.
promote_rule{T<:Number, S<:Union(Real,Complex)}(::Type{HomogPol{T}}, ::Type{S}) = 
    HomogPol{promote_type(T, S)}

immutable TaylorN{T<:Number} <: AbstractSeries{T,NUMVARS[end]}
    coeffs  :: Array{HomogPol{T},1}
    order   :: Int
    function TaylorN( v::Array{HomogPol{T},1}, mOrder::Int )
        ll = length(v)
        @inbounds maxOrd = max( [v[i].order for i=1:ll]..., mOrder )
        #@assert maxOrd <= MAXORDER[end]
        coeffs = zeros(HomogPol{T}, maxOrd)
        for i = 1:ll
            @inbounds ord = v[i].order
            @inbounds coeffs[ord+1] += v[i]
        end
        new(coeffs, maxOrd)
    end
end
TaylorN{T<:Number}(x::TaylorN{T}, order::Int) = TaylorN{T}(x.coeffs, order )
TaylorN{T<:Number}(x::TaylorN{T}) = TaylorN{T}(x.coeffs, x.order )
TaylorN{T<:Number}(v::Array{HomogPol{T},1}, order::Int) = TaylorN{T}(v, order )
TaylorN{T<:Number}(v::Array{HomogPol{T},1}) = TaylorN{T}(v, 0 )
TaylorN{T<:Number}(x::HomogPol{T}, order::Int) = TaylorN{T}([x], order )
TaylorN{T<:Number}(x::HomogPol{T}) = TaylorN{T}([x], x.order )
TaylorN{T<:Number}(x::T, order::Int) = TaylorN{T}([HomogPol(x)], order )
TaylorN{T<:Number}(x::T) = TaylorN{T}([HomogPol(x)], 0 )

## Type, length ##
eltype{T<:Number}(::TaylorN{T}) = T
length(a::TaylorN) = length( a.coeffs )
get_numVars(x::TaylorN) = NUMVARS[end]
get_maxOrder(x::TaylorN) = x.order

## zero and one ##
zero{T<:Number}(a::TaylorN{T}) = TaylorN(zero(T), a.order)
one{T<:Number}(a::TaylorN{T}) = TaylorN(one(T), a.order)

## Conversion and promotion rules ##
convert{T<:Number}(::Type{TaylorN{T}}, a::TaylorN) = 
    TaylorN( convert(Array{HomogPol{T}}, a.coeffs), a.order )
convert{T<:Number, S<:Number}(::Type{TaylorN{T}}, b::HomogPol{S}) = 
    TaylorN( convert(HomogPol{T}, b) )
convert{T<:Number, S<:Number}(::Type{TaylorN{T}}, b::Array{HomogPol{S},1}) = 
    TaylorN( convert(Array{HomogPol{T},1}, b) )
convert{T<:Number}(::Type{TaylorN{T}}, b::Number) = TaylorN( convert(T, b), 0)
#
promote_rule{T<:Number, S<:Number}(::Type{TaylorN{T}}, ::Type{TaylorN{S}}) = 
    TaylorN{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{TaylorN{T}}, ::Type{HomogPol{S}}) = 
    TaylorN{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{TaylorN{T}}, ::Type{Array{HomogPol{S},1}}) = 
    TaylorN{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{TaylorN{T}}, ::Type{S}) = TaylorN{promote_type(T, S)}

## Auxiliary function ##
function fixshape{T<:Number, S<:Number}(a::HomogPol{T}, b::HomogPol{S})
    #@assert a.order == b.order
    a1, b1 = promote(a, b)
    return HomogPol(a1, a.order), HomogPol(b1, a.order), a.order
end
function fixshape{T<:Number, S<:Number}(a::TaylorN{T}, b::TaylorN{S})
    #?#@assert a.numVars == b.numVars
    order = max(a.order, b.order)
    a1, b1 = promote(a, b)
    return TaylorN(a1, order), TaylorN(b1, order), order
end

## real, imag, conj and ctranspose ##
for f in (:real, :imag, :conj)
    @eval ($f){T<:Number}(a::HomogPol{T}) = HomogPol(($f)(a.coeffs), a.order)
    @eval ($f){T<:Number}(a::TaylorN{T}) = TaylorN(($f)(a.coeffs), a.order)
end
ctranspose{T<:Number}(a::HomogPol{T}) = conj(a)
ctranspose{T<:Number}(a::TaylorN{T}) = conj(a)

## Equality ##
==(a::HomogPol, b::HomogPol) = a.coeffs == b.coeffs
function ==(a::TaylorN, b::TaylorN)
    a1, b1, order = fixshape(a, b)
    a1.coeffs == b1.coeffs
end

## Addition and substraction ##
for f in (:+, :-), T in (:HomogPol, :TaylorN)
    @eval begin
        function ($f)(a::($T), b::($T))
            a1, b1, order = fixshape(a,b)
            v = ($f)(a1.coeffs, b1.coeffs)
            return ($T)(v, order)
        end
        ($f)(a::($T)) = ($T)(($f)(a.coeffs), a.order)
    end
end

## Multiplication ##
function *(a::HomogPol, b::HomogPol)
    a1, b1 = promote(a, b)
    order = a1.order + b1.order
    order > MAXORDER[end] && return zero(HomogPol(a1, MAXORDER[end]))
    @inbounds nCoefHa = sizeTable[a1.order+1]
    @inbounds nCoefHb = sizeTable[b1.order+1]
    @inbounds nCoefH  = sizeTable[order+1]
    T = eltype(a1)
    coeffs = zeros(T, nCoefH)
    for na = 1:nCoefHa
        @inbounds inda = indicesTable[a1.order+1][na]
        for nb = 1:nCoefHb
            @inbounds indb = indicesTable[b1.order+1][nb]
            @inbounds pos = posTable[order+1][inda+indb]
            @inbounds coeffs[pos] += a1.coeffs[na] * b1.coeffs[nb]
        end
    end
    HomogPol(coeffs, order)
end
function *(a::TaylorN, b::TaylorN)
    order = min(a.order + b.order, MAXORDER[end])
    T = promote_type( eltype(a), eltype(b) )
    coeffs = zeros(HomogPol{T}, order)
    a = TaylorN(a, order)
    b = TaylorN(b, order)
    @inbounds coeffs[1] = a.coeffs[1] * b.coeffs[1]
    for ord = 1:order
        for i = 0:ord
            j = ord-i
            @inbounds coeffs[ord+1] += a.coeffs[i+1] * b.coeffs[j+1]
        end
    end
    TaylorN(coeffs, order)
end

## Division ##
/(a::HomogPol, x::Real) = a*inv(x)
/(a::HomogPol, x::Complex) = a*inv(x)
function /(a::TaylorN, b::TaylorN)
    b0 = b.coeffs[1].coeffs[1]
    #@assert b0 != zero(b0)
    #!?orddivfact, cdivfact = divfactorization(a1, b1) # order and coefficient of first factorized term
    cdivfact = a.coeffs[1] / b0
    order = max(a.order, b.order)
    a = TaylorN(a, order)
    b = TaylorN(b, order)
    T = eltype(cdivfact)
    coeffs = zeros(HomogPol{T}, order)
    @inbounds coeffs[1] = cdivfact
    for ord = 1:order
        for i = 0:ord-1
            j = ord-i
            @inbounds coeffs[ord+1] = coeffs[i+1] * b.coeffs[j+1]
        end
        @inbounds coeffs[ord+1] = (a.coeffs[ord+1] - coeffs[ord+1]) / b0
    end
    TaylorN(coeffs, order)
end
# function divfactorization(a1::Taylor, b1::Taylor)
#     # order of first factorized term; a1 and b1 are assumed to be of the same order (length)
#     a1nz = firstnonzero(a1)
#     b1nz = firstnonzero(b1)
#     orddivfact = min(a1nz, b1nz)
#     if orddivfact > a1.order
#         orddivfact = a1.order
#     end
#     cdivfact = a1.coeffs[orddivfact+1] / b1.coeffs[orddivfact+1]
#     aux = abs2(cdivfact)
#     # Is the polynomial factorizable?
#     if isinf(aux) || isnan(aux)
#         info("Order k=$(orddivfact) => coeff[$(orddivfact+1)]=$(cdivfact)")
#         error("Division does not define a Taylor polynomial\n",
#             " or its first non-zero coefficient is Inf/NaN.\n")
#     ##else orddivfact>0
#     ##    warn("Factorizing the polynomial.\n",
#     ##        "The last k=$(orddivfact) Taylor coefficients ARE SET to 0.\n")
#     end
#     return orddivfact, cdivfact
# end

## Division functions: rem and mod
for op in (:mod, :rem)
    @eval begin
        function ($op){T<:Real}(a::TaylorN{T}, x::Real)
            coeffs = a.coeffs
            @inbounds coeffs[1] = HomogPol(($op)(a.coeffs[1].coeffs[1], x))
            return TaylorN( coeffs, a.order )
        end
    end
end
function mod2pi{T<:Real}(a::TaylorN{T}) 
    coeffs = a.coeffs
    @inbounds coeffs[1] = HomogPol(mod2pi( a.coeffs[1].coeffs[1] ))
    return TaylorN( coeffs, a.order )
end

## Int power ##
function ^(a::HomogPol, n::Integer)
    #@assert n >= 0
    n == 0 && return HomogPol( one(eltype(a)) )
    if n%2 == 0     # even power
        n == 2 && return square(a)
        pow = div(n, 2)
        return square( a^pow )
    else            # odd power
        n == 1 && return a
        pow = div(n-1, 2)
        return a * square( a^pow )
    end
end
function ^(a::TaylorN, n::Integer)
    uno = one(eltype(a))
    n < 0 && return uno / a^(-n)
    n == 0 && return TaylorN( uno )
    if n%2 == 0     # even power
        n == 2 && return square(a)
        pow = div(n, 2)
        return square( a^pow )
    else            # odd power
        n == 1 && return a
        pow = div(n-1, 2)
        return a * square( a^pow )
    end
end
## Rational power ##
^(a::TaylorN, x::Rational) = a^(x.num/x.den)
## Real power ##
function ^(a::TaylorN, x::Real)
    uno = one(eltype(a))
    x == zero(x) && return TaylorN( uno )
    x == 0.5 && return sqrt(a)
    a0 = a.coeffs[1].coeffs[1]
    #@assert a0 != zero(a0)
    aux = ( a0 )^x
    T = typeof(aux)
    coeffs = zeros(HomogPol{T}, a.order)
    @inbounds coeffs[1] = HomogPol( aux )
    for ord = 1:a.order
        for i = 0:ord-1
            aux = x*(ord-i)-i
            @inbounds coeffs[ord+1] += aux * coeffs[i+1] * a.coeffs[ord-i+1]
        end
        @inbounds coeffs[ord+1] = coeffs[ord+1] / (ord*a0)
    end
    TaylorN(coeffs, a.order)
end
^{T<:Number,S<:Number}(a::TaylorN{T}, x::Complex{S}) = exp( x*log(a) )
^(a::TaylorN, b::TaylorN) = exp( b*log(a) )

## Square ##
square{T<:Number}(a::HomogPol{T}) = a * a
function square{T<:Number}(a::TaylorN{T})
    order = min(2*a.order, MAXORDER[end])
    a = TaylorN(a, order)
    coeffs = zeros(HomogPol{T}, order)
    @inbounds coeffs[1] = a.coeffs[1]^2
    for ord = 1:order
        kodd = ord%2
        kord = div(ord-2+kodd, 2)
        for i = 0: kord
            @inbounds coeffs[ord+1] += a.coeffs[i+1] * a.coeffs[ord-i+1]
        end
        @inbounds coeffs[ord+1] = 2 * coeffs[ord+1]
        ord%2 == 1 && continue
        @inbounds coeffs[ord+1] += square( a.coeffs[div(ord,2)+1] )
    end
    TaylorN(coeffs, order)
end

## sqrt ##
function sqrt(a::TaylorN)
    order = a.order
    p0 = sqrt( a.coeffs[1].coeffs[1] )
    T = typeof(p0)
    coeffs = zeros(HomogPol{T}, a.order)
    @inbounds coeffs[1] = HomogPol( p0 )
    for ord = 1:a.order
        kodd = ord%2
        kord = div(ord-2+kodd, 2)
        for i = 1:kord
            @inbounds coeffs[ord+1] += coeffs[i+1] * coeffs[ord-i+1]
        end
        @inbounds coeffs[ord+1] = a.coeffs[ord+1] - 2 * coeffs[ord+1]
        if ord%2 == 0
            @inbounds coeffs[ord+1] = coeffs[ord+1] - square( coeffs[div(ord,2)+1] )
        end
        @inbounds coeffs[ord+1] = coeffs[ord+1] / (2 * p0)
    end
    TaylorN(coeffs, a.order)
end

## exp ##
function exp(a::TaylorN)
    order = a.order
    aux = exp( a.coeffs[1].coeffs[1] )
    T = typeof(aux)
    coeffs = zeros(HomogPol{T}, order)
    @inbounds coeffs[1] = HomogPol(aux)
    for ord = 1:order
        for j = 0:ord-1
            @inbounds coeffs[ord+1] += (ord-j) * a.coeffs[ord-j+1] * coeffs[j+1]
        end
        @inbounds coeffs[ord+1] = coeffs[ord+1] / ord
    end
    TaylorN(coeffs, order)
end

## log ##
function log(a::TaylorN)
    order = a.order
    a0 = a.coeffs[1].coeffs[1]
    l0 = log( a0 )
    T = typeof(l0)
    coeffs = zeros(HomogPol{T}, order)
    @inbounds coeffs[1] = HomogPol(l0)
    for ord = 1:order
        for j = 1:ord-1
            @inbounds coeffs[ord+1] += j * a.coeffs[ord-j+1] * coeffs[j+1]
        end
        @inbounds coeffs[ord+1] = (a.coeffs[ord+1] - coeffs[ord+1] / ord ) / a0
    end
    TaylorN(coeffs, order)
end

## sin and cos ##
sin(a::TaylorN) = sincos(a)[1]
cos(a::TaylorN) = sincos(a)[2]
function sincos(a::TaylorN)
    order = a.order
    a0 = a.coeffs[1].coeffs[1]
    s0 = sin( a0 )
    c0 = cos( a0 )
    T = typeof(s0)
    coeffsSin = zeros(HomogPol{T}, order)
    coeffsCos = zeros(HomogPol{T}, order)
    @inbounds coeffsSin[1] = HomogPol(s0)
    @inbounds coeffsCos[1] = HomogPol(c0)
    for ord = 1:order
        for j = 0:ord-1
            @inbounds coeffsSin[ord+1] += (ord-j) * a.coeffs[ord-j+1] * coeffsCos[j+1]
            @inbounds coeffsCos[ord+1] += (ord-j) * a.coeffs[ord-j+1] * coeffsSin[j+1]
        end
        @inbounds coeffsSin[ord+1] =  coeffsSin[ord+1] / ord
        @inbounds coeffsCos[ord+1] = -coeffsCos[ord+1] / ord
    end
    TaylorN(coeffsSin, order), TaylorN(coeffsCos, order)
end

## tan ##
# function tan(a::TaylorN)

## Differentiation ##
"""Partial differentiation of a HomogPol series with respect to the r-th variable"""
function diffTaylor{T<:Number}(a::HomogPol{T}, r::Int)
    #@assert 1 <= r <= NUMVARS[end]
    order = a.order
    order == 0 && return HomogPol(zero(T))
    nCoefH = sizeTable[order]
    coeffs = zeros(T,nCoefH)
    jind = zeros(Int, NUMVARS[end])
    jind[r] = 1
    for i = 1:sizeTable[order+1]
        iind = indicesTable[order+1][i]
        n = iind[r]
        n == 0 && continue
        pos = posTable[order][iind-jind]
        coeffs[pos] = n * a.coeffs[i]
    end
    HomogPol(coeffs, order-1)
end
"""Partial differentiation of a TaylorN series with respect to the r-th variable"""
function diffTaylor{T<:Number}(a::TaylorN{T}, r::Int)
    order = a.order
    coeffs = zeros(HomogPol{T}, a.order)
    for ord = 1:a.order
        coeffs[ord] = diffTaylor( a.coeffs[ord+1], r)
    end
    return TaylorN( coeffs, a.order )
end
diffTaylor(a::TaylorN) = diffTaylor(a, 1)

# ## TO BE DONE: Integration...

## Evaluates a Taylor polynomial on a given point ##
# NEEDS REVISION since yields results not so exact
function evalHomog{T<:Number,S<:Number}(a::HomogPol{T}, vals::Array{S,1} )
    numVars = NUMVARS[end]
    #@assert length(vals) == numVars
    R = promote_type(T,S)
    suma = zero(R)
    order = a.order
    nCoefH = sizeTable[order+1]
    for pos = 1:nCoefH
        @inbounds iIndices = indicesTable[order+1][pos]
        @inbounds c = a.coeffs[pos]
        c == zero(R) && continue
        for k = 1:numVars
            @inbounds c *= vals[k]^iIndices[k]
        end
        suma += c
    end
    suma
end
function evalTaylor{T<:Number,S<:Number}(a::TaylorN{T}, vals::Array{S,1} )
    #@assert length(vals) == NUMVARS[end]
    R = promote_type(T,S)
    suma = zero(R)
    for ord = MAXORDER[end]:-1:0
        @inbounds polH = a.coeffs[ord+1]
        suma += evalHomog( polH, vals )
    end
    suma
end
evalTaylor{T<:Number}(a::TaylorN{T}) = evalTaylor(a, zeros(T, NUMVARS[end]) )
