# utils_TaylorNv2.jl: N-variables Taylor expansions
#
# Last modification: 2014.06.17
#
# Luis Benet & David P. Sanders
# UNAM
#


## Default values for the maximum degree of polynomials (MAXORDER <--> nPart) and 
##  number of variables considered (NUMVARS <--> lNiv)
const MAXORDER = [6]
const NUMVARS  = [3]
info( string("\n MAXORDER = ", MAXORDER[end], "\n NUMVARS  = ",NUMVARS[end]) )

## Hash tables
#=
  `generateIndicesTable`: generates the array of dictionaries `indicesTable`. Then, `indicesTable[k+1]`
  contains the dictionary with the table defining the indices from the position for the homogeneous
  polynomial of degree k. 
  The number of coefficients of the homogenous polynomial of degree k is stored in `sizeTable[k+1]`.
  `generatePosTable`: generates the dictionary `posTable`; is the "inverse" of `indicesTable`
=#
function generateIndicesTable()
    numVars = NUMVARS[end]
    maxOrd = MAXORDER[end]
    arrayDDic = Dict{Int,Array{Int,1}}[]
    arraySize = Int[]
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
        @assert length(DDic) == nCoefH
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
    for kDeg = 0:maxOrd
        DDic = Dict{Array{Int,1}, Int}()
        nCoefH = sizeTable[kDeg+1]
        for i = 1:nCoefH
            v = indicesTable[kDeg+1][i]
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
    @assert n >= 0
    MAXORDER[end] = n
    info(string("MAXORDER is now ", n, "; hash tables regenerated.\n"))
    indicesTable[:], sizeTable[:] = generateIndicesTable()
    posTable[:] = generatePosTable()
    return n
end
#
get_numVars() = NUMVARS[end]
function set_numVars(n::Int)
    @assert n > 0
    n==1 && info( string("Using `Taylor` constructor rather than `TaylorN` is\n", 
        "MUCH faster for 1-variable expansions.\n") )
    NUMVARS[end] = n
    info(string("NUMVARS is now ", n, "; hash tables regenerated.\n"))
    indicesTable[:], sizeTable[:] = generateIndicesTable()
    posTable[:] = generatePosTable()
    return n
end

## HomogPol (homogeneous polynomial) constructors ##
abstract AbstractSeries{T<:Number,N} <: Number
immutable HomogPol{T<:Number} <: AbstractSeries{T,NUMVARS[end]}
    coeffs  :: Array{T,1}
    order   :: Int
    function HomogPol(coeffs::Array{T,1}, order::Int)
        @assert order <= MAXORDER[end]
        lencoef = length( coeffs )
        @inbounds nCoefH = sizeTable[order+1]
        @assert lencoef <= nCoefH
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

Base.eltype{T<:Number}(::HomogPol{T}) = T
Base.length(a::HomogPol) = length( a.coeffs )
orderH(a::HomogPol) = a.order
function orderH{T}(coeffs::Array{T,1})
    ord=0
    ll = length(coeffs)
    for i=1:MAXORDER[end]+1
        nCoefH = sizeTable[i]
        ll <= nCoefH && break
        ord += 1
    end
    return ord
end

## zero and one ##
function Base.zero{T<:Number}(a::HomogPol{T})
    a.order == 0 && return HomogPol(zero(T), 0)
    nCoefH = sizeTable[a.order+1]
    HomogPol(zeros(T,nCoefH), a.order)
end
function Base.zeros{T<:Number}(a::HomogPol{T}, maxOrd::Int)
    order = max( maxOrd, a.order)
    @assert 0 <= maxOrd <= MAXORDER[end]
    v = HomogPol{T}[]
    for ord = 0:order
        z = HomogPol(zero(T),ord)
        push!(v,z)
    end
    return v
end
Base.zeros{T<:Number}(::Type{HomogPol{T}}, maxOrd::Int) = zeros( HomogPol(zero(T)), maxOrd)
function Base.one{T<:Number}(a::HomogPol{T})
    a.order == 0 && HomogPol(one(T), 0)
    nCoefH = sizeTable[a.order+1]
    HomogPol(ones(T,nCoefH), a.order)
end
function Base.ones{T<:Number}(::HomogPol{T}, maxOrd::Int)
    @assert 0 <= maxOrd <= MAXORDER[end]
    v = HomogPol{T}[]
    for ord = 0:maxOrd
        z = HomogPol(one(T),ord)
        push!(v,z)
    end
    return v
end
Base.ones{T<:Number}(::Type{HomogPol{T}}, maxOrd::Int) = ones( HomogPol(one(T)), maxOrd)

## Conversion and promotion rules ##
Base.convert{T<:Number}(::Type{HomogPol{T}}, a::HomogPol) = HomogPol(convert(Array{T,1}, a.coeffs), a.order)
Base.convert{T<:Number, S<:Number}(::Type{HomogPol{T}}, b::Array{S,1}) = HomogPol(convert(Array{T,1},b))
Base.convert{T<:Number}(::Type{HomogPol{T}}, b::Number) = HomogPol([convert(T,b)], 0)
Base.promote_rule{T<:Number, S<:Number}(::Type{HomogPol{T}}, ::Type{HomogPol{S}}) = HomogPol{promote_type(T, S)}
Base.promote_rule{T<:Number, S<:Number}(::Type{HomogPol{T}}, ::Type{Array{S,1}}) = HomogPol{promote_type(T, S)}
Base.promote_rule{T<:Number, S<:Number}(::Type{HomogPol{T}}, ::Type{S}) = HomogPol{promote_type(T, S)}

immutable TaylorN{T<:Number} <: AbstractSeries{T,NUMVARS[end]}
    coeffs  :: Array{HomogPol{T},1}
    order   :: Int
    numVars :: Int
    function TaylorN( v::Array{HomogPol{T},1}, mOrder::Int, nVars::Int )
        numVars = NUMVARS[end]
        ll = length(v)
        maxOrd = max( [v[i].order for i=1:ll]..., mOrder )
        @assert maxOrd <= MAXORDER[end]
        println(maxOrd)
        coeffs = zeros(HomogPol{T}, maxOrd)
        for i = 1:ll
            ord = v[i].order
            coeffs[ord+1] += v[i]
        end
        new(coeffs, maxOrd, numVars)
    end
end
TaylorN{T<:Number}(x::TaylorN{T}, order::Int) = TaylorN{T}(x.coeffs, order, NUMVARS[end])
TaylorN{T<:Number}(x::TaylorN{T}) = TaylorN{T}(x.coeffs, x.order, NUMVARS[end])
TaylorN{T<:Number}(v::Array{HomogPol{T},1}, order::Int) = TaylorN{T}(v, order, NUMVARS[end])
TaylorN{T<:Number}(v::Array{HomogPol{T},1}) = TaylorN{T}(v, 0, NUMVARS[end])
TaylorN{T<:Number}(x::HomogPol{T}, order::Int) = TaylorN{T}([x], order, NUMVARS[end])
TaylorN{T<:Number}(x::HomogPol{T}) = TaylorN{T}([x], x.order, NUMVARS[end])
TaylorN{T<:Number}(x::T, order::Int) = TaylorN{T}([HomogPol(x)], order, NUMVARS[end])
TaylorN{T<:Number}(x::T) = TaylorN{T}([HomogPol(x)], 0, NUMVARS[end])

## zero and one ##
Base.zero{T<:Number}(a::TaylorN{T}) = TaylorN(zero(T), a.order)
Base.one{T<:Number}(a::TaylorN{T}) = TaylorN(one(T), a.order)

## Auxiliary function ##
function fixshape{T<:Number, S<:Number}(a::HomogPol{T}, b::HomogPol{S})
    @assert a.order == b.order
    a1, b1 = promote(a, b)
    return HomogPol(a1, a.order), HomogPol(b1, a.order), a.order
end

## real, imag, conj and ctranspose ##
for f in (:real, :imag, :conj)
    @eval ($f){T<:Number}(a::HomogPol{T}) = HomogPol(($f)(a.coeffs), a.order)
end
ctranspose{T<:Number}(a::HomogPol{T}) = conj(a)

## Equality ##
==(a::HomogPol, b::HomogPol) = a.coeffs == b.coeffs

## Addition and substraction ##
for f in (:+, :-)
    @eval begin
        function ($f)(a::HomogPol, b::HomogPol)
            a1, b1, order = fixshape(a,b)
            v = ($f)(a1.coeffs, b1.coeffs)
            return HomogPol(v, order)
        end
        ($f)(a::HomogPol) = HomogPol(($f)(a.coeffs), a.order)
    end
end

## Multiplication ##
function *(a::HomogPol, b::HomogPol)
    a1, b1 = promote(a, b)
    order = a.order + b.order
    T = eltype(a1)
    @assert order <= MAXORDER[end]
    nCoefHa = sizeTable[a.order+1]
    nCoefHb = sizeTable[b.order+1]
    nCoefH = sizeTable[order+1]
    coeffs = zeros(T, nCoefH)
    for na = 1:nCoefHa
        inda = indicesTable[a.order+1][na]
        for nb = 1:nCoefHb
            indb = indicesTable[b.order+1][nb]
            pos = posTable[order+1][inda+indb]
            coeffs[pos] += a.coeffs[na] * b.coeffs[nb]
        end
    end
    HomogPol(coeffs, order)
end

# ## Multiplication ##
# function *(a::TaylorN, b::TaylorN)
#     a1, b1, order = fixshape(a, b)
#     T = eltype(a1)
#     nCoefTot = sizeTable[end]
#     coeffs = zeros(T, nCoefTot)
#     coeffs[1] = a1.coeffs[1] * b1.coeffs[1]   ## 0th-order
#     for k = 1:order
#         @inbounds posI, posF = posAuxTable[end][k]
#         @inbounds coeffs[posI:posF] = mulHomogCoefN( k, a1.coeffs, b1.coeffs )
#     end
#     TaylorN(coeffs, order)
# end
# function mulHomogCoefN{T<:Number}( k::Int, ac::Array{T,1}, bc::Array{T,1} )
#     k==0 && return ac[1] * bc[1]
#     @inbounds posI, posF = posAuxTable[end][k]
#     numCoefk = posF - posI + 1
#     coeffs = zeros(T, numCoefk)
#     for ka = 0:k
#         kb = k-ka
#         @inbounds posaI, posaF = posAuxTable[end][ka]
#         @inbounds posbI, posbF = posAuxTable[end][kb]
#         for ja = posaI:posaF
#             @inbounds inda = indicesTable[end][ja]
#             for jb = posbI:posbF
#                 @inbounds indb = indicesTable[end][jb]
#                 @inbounds pos = posTable[end][ inda + indb ]
#                 pos = pos - posI + 1
#                 @inbounds coeffs[pos] += ac[ja] * bc[jb]
#             end
#         end
#     end
#     return coeffs
# end

# ## Division ##
# function /(a::TaylorN, b::TaylorN)
#     a1, b1, order = fixshape(a, b)
#     #!?orddivfact, cdivfact = divfactorization(a1, b1) # order and coefficient of first factorized term
#     @assert b1.coeffs[1] != zero(b1.coeffs[1])
#     cLHopital = a1.coeffs[1] / b1.coeffs[1]
#     T = typeof(cLHopital)
#     v1 = convert(Array{T,1}, a1.coeffs)
#     v2 = convert(Array{T,1}, b1.coeffs)
#     nCoefTot = sizeTable[end]
#     coeffs = zeros(T, nCoefTot)
#     coeffs[1] = cLHopital
#     for k = 1:order
#         @inbounds posI, posF = posAuxTable[end][k]
#         @inbounds coeffs[posI:posF] = divHomogCoefN( k, v1, v2, coeffs, 0 )
#     end
#     TaylorN(coeffs, order)
# end
# # function divfactorization(a1::Taylor, b1::Taylor)
# #     # order of first factorized term; a1 and b1 are assumed to be of the same order (length)
# #     a1nz = firstnonzero(a1)
# #     b1nz = firstnonzero(b1)
# #     orddivfact = min(a1nz, b1nz)
# #     if orddivfact > a1.order
# #         orddivfact = a1.order
# #     end
# #     cdivfact = a1.coeffs[orddivfact+1] / b1.coeffs[orddivfact+1]
# #     aux = abs2(cdivfact)
# #     # Is the polynomial factorizable?
# #     if isinf(aux) || isnan(aux)
# #         info("Order k=$(orddivfact) => coeff[$(orddivfact+1)]=$(cdivfact)")
# #         error("Division does not define a Taylor polynomial\n",
# #             " or its first non-zero coefficient is Inf/NaN.\n")
# #     ##else orddivfact>0
# #     ##    warn("Factorizing the polynomial.\n",
# #     ##        "The last k=$(orddivfact) Taylor coefficients ARE SET to 0.\n")
# #     end
# #     return orddivfact, cdivfact
# # end
# # Homogeneous coefficient for the division
# function divHomogCoefN{T<:Number}( k::Int, ac::Array{T,1}, bc::Array{T,1}, 
#     coeffs::Array{T,1}, ordLHopital::Int)
#     #
#     k == ordLHopital && return ac[ordLHopital+1] / bc[ordLHopital+1]
#     @inbounds posI, posF = posAuxTable[end][k]
#     numCoefk = posF - posI + 1
#     cc = zeros(T, numCoefk)
#     @inbounds cc = mulHomogCoefN(k, coeffs, bc)
#     @inbounds cc[1:end] = (ac[posI:posF]-cc[1:end]) / bc[ordLHopital+1]
#     cc
# end

# ## Division functions: rem and mod
# for op in (:mod, :rem)
#     @eval begin
#         function ($op){T<:Real}(a::TaylorN{T}, x::Real)
#             coeffs = a.coeffs
#             coeffs[1] = ($op)(a.coeffs[1], x)
#             return TaylorN( coeffs, a.order )
#         end
#     end
# end
# function mod2pi{T<:Real}(a::TaylorN{T}) 
#     coeffs = a.coeffs
#     coeffs[1] = mod2pi( a.coeffs[1] )
#     return TaylorN( coeffs, a.order )
# end

# ## Int power ##
# function ^(a::TaylorN, n::Integer)
#     uno = one(a)
#     n < 0 && return uno / a^(-n)
#     n == 0 && return uno
#     if n%2 == 0     # even power
#         n == 2 && return square(a)
#         pow = div(n, 2)
#         return square( a^pow )
#     else            # odd power
#         n == 1 && return a
#         pow = div(n-1, 2)
#         return a*square( a^pow )
#     end
# end
# ## Rational power ##
# ^(a::TaylorN,x::Rational) = a^(x.num/x.den)
# ## Real power ##
# function ^(a::TaylorN, x::Real)
#     uno = one(a)
#     x == zero(x) && return uno
#     x == 0.5 && return sqrt(a)
#     order = a.order
#     # First non-zero coefficient
#     # ...
#     @assert a.coeffs[1] != zero(a.coeffs[1])
#     aux = ( a.coeffs[1] )^x
#     T = typeof(aux)
#     v = convert(Array{T,1}, a.coeffs)
#     nCoefTot = sizeTable[end]
#     coeffs = zeros(T, nCoefTot)
#     coeffs[1] = aux
#     for k = 1:order
#         @inbounds posI, posF = posAuxTable[end][k]
#         @inbounds coeffs[posI:posF] = powHomogCoefN(k, v, x, coeffs, 0)
#     end
#     TaylorN(coeffs,order)
# end
# function powHomogCoefN{T<:Number}(k::Int, ac::Array{T,1}, x::Real, 
#     coeffs::Array{T,1}, knull::Int)
#     #
#     k == knull && return (ac[ordLHopital+1])^x
#     @inbounds posI, posF = posAuxTable[end][k]
#     numCoefk = posF - posI + 1
#     cc = zeros(T, numCoefk)
#     for ka = 0:k-1
#         kb = k-ka
#         @inbounds posaI, posaF = posAuxTable[end][ka]
#         @inbounds posbI, posbF = posAuxTable[end][kb]
#         for ja = posaI:posaF
#             @inbounds inda = indicesTable[end][ja]
#             for jb = posbI:posbF
#                 @inbounds indb = indicesTable[end][jb]
#                 @inbounds pos = posTable[end][ inda + indb ]
#                 pos = pos - posI + 1
#                 @inbounds cc[pos] += ( x*kb-ka ) * ac[jb] * coeffs[ja]
#             end
#         end
#     end
#     @inbounds cc[1:end] = cc[1:end] / ( k * ac[1] )
#     cc
# end
# ^{T<:Number,S<:Number}(a::TaylorN{T}, x::Complex{S}) = exp( x*log(a) )
# ^(a::TaylorN, b::TaylorN) = exp( b*log(a) )

# ## Square ##
# function square{T<:Number}(a::TaylorN{T})
#     order = a.order
#     nCoefTot = sizeTable[end]
#     coeffs = zeros(T, nCoefTot)
#     coeffs[1] = a.coeffs[1]^2
#     for k = 1:order
#         @inbounds posI, posF = posAuxTable[end][k]
#         @inbounds coeffs[posI:posF] = squareHomogCoefN( k, a.coeffs )
#     end
#     TaylorN(coeffs,order)
# end
# # Homogeneous coefficients for square
# function squareHomogCoefN{T<:Number}(k::Int, ac::Array{T,1})
#     k == 0 && return (ac[1])^2
#     @inbounds posI, posF = posAuxTable[end][k]
#     numCoefk = posF - posI + 1
#     coeffs = zeros(T, numCoefk)
#     kodd = k%2
#     kend = div(k - 2 + kodd, 2)
#     for ka = 0:kend
#         kb = k - ka
#         @inbounds posaI, posaF = posAuxTable[end][ka]
#         @inbounds posbI, posbF = posAuxTable[end][kb]
#         for ja = posaI:posaF
#             @inbounds inda = indicesTable[end][ja]
#             for jb = posbI:posbF
#                 @inbounds indb = indicesTable[end][jb]
#                 @inbounds pos = posTable[end][ inda + indb ]
#                 pos = pos - posI + 1
#                 @inbounds coeffs[pos] += ac[ja] * ac[jb]
#             end
#         end
#     end
#     @inbounds coeffs[1:end] = 2coeffs[1:end]
#     #
#     if kodd == 0
#         ka = div(k,2)
#         @inbounds posaI, posaF = posAuxTable[end][ka]
#         for ja = posaI:posaF
#             @inbounds inda = indicesTable[end][ja]
#             for jb = posaI:posaF
#                 @inbounds indb = indicesTable[end][jb]
#                 @inbounds pos = posTable[end][ inda + indb ]
#                 pos = pos - posI + 1
#                 @inbounds coeffs[pos] += ac[ja] * ac[jb]
#             end
#         end
#     end
#     coeffs
# end

# ## sqrt ##
# function sqrt(a::TaylorN)
#     order = a.order
#     nCoefTot = sizeTable[end]
#     aux = sqrt(a.coeffs[1])
#     T = typeof(aux)
#     v = convert(Array{T,1}, a.coeffs)
#     coeffs = zeros(T, nCoefTot)
#     coeffs[1] = aux
#     for k = 1:order
#         @inbounds posI, posF = posAuxTable[end][k]
#         @inbounds coeffs[posI:posF] = sqrtHomogCoefN( k, v, coeffs, 0 )
#     end
#     TaylorN(coeffs,order)
# end
# # Homogeneous coefficients for sqrt
# function sqrtHomogCoefN{T<:Number}(k::Int, ac::Array{T,1}, coeffs::Array{T,1}, knull::Int)
#     k == 0 && return sqrt( ac[1] )
#     @inbounds posI, posF = posAuxTable[end][k]
#     numCoefk = posF - posI + 1
#     cc = zeros(T, numCoefk)
#     kodd = k%2
#     kend = div(k - 2 + kodd, 2)
#     for ka = 1:kend
#         kb = k - ka
#         @inbounds posaI, posaF = posAuxTable[end][ka]
#         @inbounds posbI, posbF = posAuxTable[end][kb]
#         for ja = posaI:posaF
#             @inbounds inda = indicesTable[end][ja]
#             for jb = posbI:posbF
#                 @inbounds indb = indicesTable[end][jb]
#                 @inbounds pos = posTable[end][ inda + indb ]
#                 pos = pos - posI + 1
#                 @inbounds cc[pos] += coeffs[ja] * coeffs[jb]
#             end
#         end
#     end
#     @inbounds cc[1:end] = 2cc[1:end]
#     #
#     if kodd == 0
#         ka = div(k,2)
#         @inbounds posaI, posaF = posAuxTable[end][ka]
#         for ja = posaI:posaF
#             @inbounds inda = indicesTable[end][ja]
#             for jb = posaI:posaF
#                 @inbounds indb = indicesTable[end][jb]
#                 @inbounds pos = posTable[end][ inda + indb ]
#                 pos = pos - posI + 1
#                 @inbounds cc[pos] += coeffs[ja] * coeffs[jb]
#             end
#         end
#     end
#     @inbounds cc[1:end] = (ac[posI:posF] - cc[1:end]) / (2coeffs[1])
#     cc
# end

# ## exp ##
# function exp(a::TaylorN)
#     order = a.order
#     nCoefTot = sizeTable[end]
#     aux = exp(a.coeffs[1])
#     T = typeof(aux)
#     v = convert(Array{T,1}, a.coeffs)
#     coeffs = zeros(T, nCoefTot)
#     coeffs[1] = aux
#     for k = 1:order
#         @inbounds posI, posF = posAuxTable[end][k]
#         @inbounds coeffs[posI:posF] = expHomogCoefN( k, v, coeffs )
#     end
#     TaylorN(coeffs,order)
# end
# # Homogeneous coefficients for exp
# function expHomogCoefN{T<:Number}(k::Int, ac::Array{T,1}, coeffs::Array{T,1})
#     k == 0 && return exp( ac[1] )
#     @inbounds posI, posF = posAuxTable[end][k]
#     numCoefk = posF - posI + 1
#     cc = zeros(T, numCoefk)
#     for ka = 0:k-1
#         kb = k - ka
#         @inbounds posaI, posaF = posAuxTable[end][ka]
#         @inbounds posbI, posbF = posAuxTable[end][kb]
#         for ja = posaI:posaF
#             @inbounds inda = indicesTable[end][ja]
#             for jb = posbI:posbF
#                 @inbounds indb = indicesTable[end][jb]
#                 @inbounds pos = posTable[end][ inda + indb ]
#                 pos = pos - posI + 1
#                 @inbounds cc[pos] += kb * coeffs[ja] * ac[jb]
#             end
#         end
#     end
#     @inbounds cc[1:end] = cc[1:end] / k
#     cc
# end

# ## log ##
# function log(a::TaylorN)
#     order = a.order
#     nCoefTot = sizeTable[end]
#     aux = log(a.coeffs[1])
#     T = typeof(aux)
#     v = convert(Array{T,1}, a.coeffs)
#     coeffs = zeros(T, nCoefTot)
#     coeffs[1] = aux
#     for k = 1:order
#         @inbounds posI, posF = posAuxTable[end][k]
#         @inbounds coeffs[posI:posF] = logHomogCoefN( k, v, coeffs )
#     end
#     TaylorN(coeffs, order)
# end
# # Homogeneous coefficients for square
# function logHomogCoefN{T<:Number}(k::Int, ac::Array{T,1}, coeffs::Array{T,1})
#     k == 0 && return log( ac[1] )
#     @inbounds posI, posF = posAuxTable[end][k]
#     numCoefk = posF - posI + 1
#     cc = zeros(T, numCoefk)
#     for ka = 1:k-1
#         kb = k - ka
#         @inbounds posaI, posaF = posAuxTable[end][ka]
#         @inbounds posbI, posbF = posAuxTable[end][kb]
#         for ja = posaI:posaF
#             @inbounds inda = indicesTable[end][ja]
#             for jb = posbI:posbF
#                 @inbounds indb = indicesTable[end][jb]
#                 @inbounds pos = posTable[end][ inda + indb ]
#                 pos = pos - posI + 1
#                 @inbounds cc[pos] += ka * coeffs[ja] * ac[jb]
#             end
#         end
#     end
#     @inbounds cc[1:end] = ( ac[posI:posF] - cc[1:end] / k ) / ac[1]
#     cc
# end

# ## sin and cos ##
# sin(a::TaylorN) = sincos(a)[1]
# cos(a::TaylorN) = sincos(a)[2]
# function sincos(a::TaylorN)
#     order = a.order
#     nCoefTot = sizeTable[end]
#     aux = sin(a.coeffs[1])
#     T = typeof(aux)
#     v = convert(Array{T,1}, a.coeffs)
#     sincoeffs = zeros(T,nCoefTot)
#     coscoeffs = zeros(T,nCoefTot)
#     sincoeffs[1] = aux
#     coscoeffs[1] = cos( a.coeffs[1] )
#     for k = 1:order
#         @inbounds posI, posF = posAuxTable[end][k]
#         @inbounds sincoeffs[posI:posF], coscoeffs[posI:posF] = 
#             sincosHomogCoefN(k, v, sincoeffs, coscoeffs)
#     end
#     return TaylorN( sincoeffs, order ), TaylorN( coscoeffs, order )
# end
# # Homogeneous coefficients for sin and cos
# function sincosHomogCoefN{T<:Number}(k::Int, ac::Array{T,1}, scoeffs::Array{T,1}, ccoeffs::Array{T,1})
#     k == 0 && return sin( ac[1] ), cos( ac[1] )
#     @inbounds posI, posF = posAuxTable[end][k]
#     numCoefk = posF - posI + 1
#     sv = zeros(T, numCoefk)
#     cv = zeros(T, numCoefk)
#     for ka = 0:k-1
#         kb = k - ka
#         @inbounds posaI, posaF = posAuxTable[end][ka]
#         @inbounds posbI, posbF = posAuxTable[end][kb]
#         for ja = posaI:posaF
#             @inbounds inda = indicesTable[end][ja]
#             for jb = posbI:posbF
#                 @inbounds indb = indicesTable[end][jb]
#                 @inbounds pos = posTable[end][ inda + indb ]
#                 pos = pos - posI + 1
#                 @inbounds begin
#                     x = kb * ac[jb]
#                     sv[pos] += x * ccoeffs[ja]
#                     cv[pos] -= x * scoeffs[ja]
#                 end
#             end
#         end
#     end
#     @inbounds sv[1:end] = sv[1:end] / k
#     @inbounds cv[1:end] = cv[1:end] / k
#     sv, cv
# end

# ## Differentiation ##
# """Partial differentiation of a TaylorN series with respect to the r-th variable"""
# function diffTaylor{T<:Number}(a::TaylorN{T}, r::Int)
#     @assert 1 <= r <= NUMVARS[end]
#     order = a.order
#     nCoefTot = sizeTable[end]
#     coeffs = zeros(T, nCoefTot)
#     iIndices = zeros(Int, NUMVARS[end])
#     for pos = 1:nCoefTot
#         @inbounds iIndices[1:end] = indicesTable[end][pos]
#         @inbounds n = iIndices[r]
#         n == 0 && continue
#         @inbounds iIndices[r] -= 1
#         @inbounds posI = posTable[end][ iIndices ]
#         @inbounds coeffs[posI] = n * a.coeffs[pos]
#     end
#     return TaylorN( coeffs, order )
# end
# diffTaylor(a::TaylorN) = diffTaylor(a, 1)

# ## TO BE DONE: Integration...

# ## Evaluates a Taylor polynomial on a given point ##
# ## NEEDS REVISION since yields only approx results
# function evalTaylor{T<:Number,S<:Number}(a::TaylorN{T}, vals::Array{S,1} )
#     numVars = NUMVARS[end]
#     @assert length(vals) == numVars
#     R = promote_type(T,S)
#     suma = zero(R)
#     nCoefTot = sizeTable[end]
#     iIndices = zeros(Int, numVars)
#     for pos = nCoefTot:-1:1
#         @inbounds iIndices[1:end] = indicesTable[end][pos]
#         @inbounds val = a.coeffs[pos]
#         val == zero(R) && continue
#         for k = 1:numVars
#             @inbounds val *= vals[k]^iIndices[k]
#         end
#         suma += val
#     end
#     suma
# end
# evalTaylor{T<:Number}(a::TaylorN{T}) = evalTaylor(a, zeros(T, a.numVars) )
