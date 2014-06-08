# utils_TaylorN.jl: N-variables Taylor expansions
#
# Last modification: 2014.06.07
#
# Luis Benet & David P. Sanders
# UNAM
#


## Default values for the maximum degree of polynomials (MAXORDER) and 
##  number of variables considered (NUMVARS)
const MAXORDER = [4]
const NUMVARS = [2]

## Hash tables
"""
`generateIndicesTable`: generates the dictionary `indicesTable`.
NUMVARS: number of variables for the polynomial expansions <--> lNiv
MAXORDER: maximum degree of polynomials <--> nPart
`generatePosTable`: generates the dictionary `posTable`; is the "inverse" of `indicesTable`
"""
function generateIndicesTable()
    numVars::Int = NUMVARS[end]
    maxDeg::Int = MAXORDER[end]
    DDic = Dict{Int, Array{Int,1}}()
    #
    if numVars==1
        info( string("Using `Taylor` constructor rather than `TaylorN` is\n", 
            "MUCH faster for 1-variable expansions.\n") )
        for k = 0:maxDeg
            DDic[k+1] = [k]
        end
        return DDic
    end
    #
    nCoefTot = binomial( numVars+maxDeg, maxDeg)
    iindices = zeros(Int, numVars)
    pos = 0
    for kDeg = 0:maxDeg
        iV = numVars
        for iz = 0:kDeg
            iindices[end] = iz
            @inbounds pos, DDic[pos] = pos2indices!( iV, kDeg, iindices, pos, DDic)
        end
    end
    return DDic
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

const indicesTable = [ generateIndicesTable() ]
const sizeTable = [ length( indicesTable[end]) ]

function generatePosTable()
    DDic = Dict{Array{Int,1}, Int}()
    for i = 1:sizeTable[end]
        v = indicesTable[end][i]
        DDic[ v ] = i
    end
    return DDic
end
const posTable = [ generatePosTable() ]

## Functions to obtain the number of homogenous coefficients of given degree
# """Returns the number of homogeneous coefficients of degree k for NUMVARS:
#     binomial(k+NUMVARS-1,k)"""
# function numHomogCoefK(k::Int)
#     k == 0 && return 1
#     return binomial( k + NUMVARS[end] - 1, k )
# end
"""Returns the position of the first homogeneous coefficient of degree k for NUMVARS:
    binomial(k+NUMVARS-1,k-1)+1"""
function posHomogCoefK(k::Int)
    k == 0 && return 1
    return binomial( k + NUMVARS[end]-1, k-1 ) + 1
end
"""Generates the auxiliary table of the initial and final position for the coefficients of 
homogeneous polynomial of given degree (`posAuxTable`)"""
function generateAuxTable()
    DDic = Dict{Int,(Int, Int)}()
    for k=0:MAXORDER[end]
        posI = posHomogCoefK(k)
        posF = posHomogCoefK(k+1)-1
        DDic[k] = ( posI, posF )
    end
    return DDic
end
const posAuxTable = [ generateAuxTable() ]

## Utilities to get/set MAXORDER and NUMVARS; they reset the indicesTable
get_maxOrder() = MAXORDER[end]
function set_maxOrder(n::Int)
    @assert n >= 0
    MAXORDER[end] = n
    info(string("MAXORDER is now ", n, "; hash tables regenerated.\n"))
    indicesTable[end] = generateIndicesTable()
    sizeTable[end] = length( indicesTable[end] )
    posTable[end] = generatePosTable()
    posAuxTable[end] = generateAuxTable()
    return n
end
#
get_numVars() = NUMVARS[end]
function set_numVars(n::Int)
    @assert n > 0
    NUMVARS[end] = n
    info(string("NUMVARS is now ", n, "; hash tables regenerated.\n"))
    indicesTable[end] = generateIndicesTable()
    sizeTable[end] = length( indicesTable[end] )
    posTable[end] = generatePosTable()
    posAuxTable[end] = generateAuxTable()
    return n
end

## TaylorN Constructors ##
immutable TaylorN{T<:Number} <: AbstractSeries{T,NUMVARS[end]}
   coeffs :: Array{T,1}
   order :: Int
   numVars :: Int
   function TaylorN(coeffs::Array{T,1}, order::Int, numVars::Int)
        @assert order <= MAXORDER[end]
        lencoef = length(coeffs)
        nCoefTot = sizeTable[end]
        @assert lencoef <= nCoefTot
        v = zeros(T, nCoefTot)
        @inbounds v[1:lencoef] = coeffs[1:lencoef]
        new(v, order, NUMVARS[end])
   end
end
TaylorN{T<:Number}(x::TaylorN{T}, order::Int) = TaylorN{T}(x.coeffs, order, NUMVARS[end])
TaylorN{T<:Number}(x::TaylorN{T}) = TaylorN{T}(x.coeffs, x.order, NUMVARS[end])
#
TaylorN{T<:Number}(coeffs::Array{T,1}, order::Int) = TaylorN{T}(coeffs, order, NUMVARS[end])
TaylorN{T<:Number}(coeffs::Array{T,1}) = TaylorN{T}(coeffs, MAXORDER[end], NUMVARS[end])
#
TaylorN{T<:Number}(x::T, order::Int) = TaylorN{T}([x], order, NUMVARS[end])
TaylorN{T<:Number}(x::T) = TaylorN{T}([x], 0, NUMVARS[end])

## Type, length ##
eltype{T<:Number}(::TaylorN{T}) = T
length(a::TaylorN) = length( a.coeffs )
get_numVars(x::TaylorN) = x.numVars
get_maxOrder(x::TaylorN) = x.order

## Conversion and promotion rules ##
convert{T<:Number}(::Type{TaylorN{T}}, a::TaylorN) = TaylorN(convert(Array{T,1}, a.coeffs), a.order)
convert{T<:Number, S<:Number}(::Type{TaylorN{T}}, b::Array{S,1}) = TaylorN(convert(Array{T,1},b))
convert{T<:Number}(::Type{TaylorN{T}}, b::Number) = TaylorN([convert(T,b)], 0)
promote_rule{T<:Number, S<:Number}(::Type{TaylorN{T}}, ::Type{TaylorN{S}}) = TaylorN{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{TaylorN{T}}, ::Type{Array{S,1}}) = TaylorN{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{TaylorN{T}}, ::Type{S}) = TaylorN{promote_type(T, S)}

## Auxiliary functions ##
# function firstnonzero{T<:Number}(a::Taylor{T})
#     order = a.order
#     nonzero::Int = order+1
#     z = zero(T)
#     for i = 1:order+1
#         if a.coeffs[i] != z
#             nonzero = i-1
#             break
#         end
#     end
#     nonzero
# end
function fixshape{T<:Number, S<:Number}(a::TaylorN{T}, b::TaylorN{S})
    @assert a.numVars == b.numVars
    order = max(a.order, b.order)
    a1, b1 = promote(a, b)
    return TaylorN(a1, order), TaylorN(b1, order), order
end

## real, imag, conj and ctranspose ##
for f in (:real, :imag, :conj)
    @eval ($f){T<:Number}(a::TaylorN{T}) = TaylorN(($f)(a.coeffs), a.order)
end
ctranspose{T<:Number}(a::TaylorN{T}) = conj(a)

## zero and one ##
zero{T<:Number}(a::TaylorN{T}) = TaylorN(zero(T), a.order)
one{T<:Number}(a::TaylorN{T}) = TaylorN(one(T), a.order)

## Equality ##
function ==(a::TaylorN, b::TaylorN)
    a1, b1, order = fixshape(a, b)
    return a1.coeffs == b1.coeffs
end

## Addition and substraction ##
for f in (:+, :-)
    @eval begin
        function ($f)(a::TaylorN, b::TaylorN)
            a1, b1, order = fixshape(a, b)
            v = ($f)(a1.coeffs, b1.coeffs)
            return TaylorN(v, order)
        end
        ($f)(a::TaylorN) = TaylorN(($f)(a.coeffs), a.order)
    end
end

## Multiplication ##
function *(a::TaylorN, b::TaylorN)
    a1, b1, order = fixshape(a, b)
    T = eltype(a1)
    nCoefTot = sizeTable[end]
    coeffs = zeros(T, nCoefTot)
    coeffs[1] = a1.coeffs[1] * b1.coeffs[1]   ## 0th-order
    for k = 1:order
        @inbounds posI, posF = posAuxTable[end][k]
        @inbounds coeffs[posI:posF] = mulHomogCoefN( k, a1.coeffs, b1.coeffs )
    end
    TaylorN(coeffs, order)
end
function mulHomogCoefN{T<:Number}( k::Int, ac::Array{T,1}, bc::Array{T,1} )
    k==0 && return ac[1] * bc[1]
    @inbounds posI, posF = posAuxTable[end][k]
    numCoefk = posF - posI + 1
    coeffs = zeros(T, numCoefk)
    for ka = 0:k
        kb = k-ka
        @inbounds posaI, posaF = posAuxTable[end][ka]
        @inbounds posbI, posbF = posAuxTable[end][kb]
        for ja = posaI:posaF
            @inbounds inda = indicesTable[end][ja]
            for jb = posbI:posbF
                @inbounds indb = indicesTable[end][jb]
                @inbounds pos = posTable[end][ inda + indb ]
                pos = pos - posI + 1
                @inbounds coeffs[pos] += ac[ja] * bc[jb]
            end
        end
    end
    return coeffs
end

## Division ##
function /(a::TaylorN, b::TaylorN)
    a1, b1, order = fixshape(a, b)
    #!?ordLHopital, cLHopital = divlhopital(a1, b1) # L'Hôpital order and coefficient
    @assert b1.coeffs[1] != zero(b1.coeffs[1])
    cLHopital = a1.coeffs[1] / b1.coeffs[1]
    T = typeof(cLHopital)
    v1 = convert(Array{T,1}, a1.coeffs)
    v2 = convert(Array{T,1}, b1.coeffs)
    nCoefTot = sizeTable[end]
    coeffs = zeros(T, nCoefTot)
    coeffs[1] = cLHopital
    for k = 1:order
        @inbounds posI, posF = posAuxTable[end][k]
        @inbounds coeffs[posI:posF] = divHomogCoefN( k, v1, v2, coeffs, 0 )
    end
    TaylorN(coeffs, order)
end
# function divlhopital(a1::TaylorN, b1::TaylorN)
#     # L'Hôpital order is calculated; a1 and b1 are assumed to be of the same order (length)
#     a1nz = firstnonzero(a1)
#     b1nz = firstnonzero(b1)
#     ordLHopital = min(a1nz, b1nz)
#     if ordLHopital > a1.order
#         ordLHopital = a1.order
#     end
#     cLHopital = a1.coeffs[ordLHopital+1] / b1.coeffs[ordLHopital+1]
#     aux = abs2(cLHopital)
#     # Can L'Hôpital be applied?
#     if isinf(aux)
#         info("Order k=$(ordLHopital) => coeff[$(ordLHopital+1)]=$(cLHopital)")
#         error("Division does not define a Taylor polynomial or its first coefficient is infinite.\n")
#     elseif isnan(aux)
#         info("Order k=$(ordLHopital) => coeff[$(ordLHopital+1)]=$(cLHopital)")
#         error("Impossible to apply L'Hôpital...\n")
#     elseif ordLHopital>0
#         warn("Applying L'Hôpital. The last k=$(ordLHopital) Taylor coefficients ARE SET to 0.\n")
#     end
#     return ordLHopital, cLHopital
# end
# Homogeneous coefficient for the division
function divHomogCoefN{T<:Number}( k::Int, ac::Array{T,1}, bc::Array{T,1}, 
    coeffs::Array{T,1}, ordLHopital::Int)
    #
    k == ordLHopital && return ac[ordLHopital+1] / bc[ordLHopital+1]
    @inbounds posI, posF = posAuxTable[end][k]
    numCoefk = posF - posI + 1
    cc = zeros(T, numCoefk)
    @inbounds cc = mulHomogCoefN(k, coeffs, bc)
    @inbounds cc[1:end] = (ac[posI:posF]-cc[1:end]) / bc[ordLHopital+1]
    cc
end

## Division functions: rem and mod
for op in (:mod, :rem)
    @eval begin
        function ($op){T<:Real}(a::TaylorN{T}, x::Real)
            coeffs = a.coeffs
            coeffs[1] = ($op)(a.coeffs[1], x)
            return TaylorN( coeffs, a.order )
        end
    end
end
function mod2pi{T<:Real}(a::TaylorN{T}) 
    coeffs = a.coeffs
    coeffs[1] = mod2pi( a.coeffs[1] )
    return TaylorN( coeffs, a.order )
end

## Int power ##
function ^(a::TaylorN, n::Integer)
    uno = one(a)
    n < 0 && return uno / a^(-n)
    n == 0 && return uno
    if n%2 == 0     # even power
        n == 2 && return square(a)
        pow = div(n, 2)
        return square( a^pow )
    else            # odd power
        n == 1 && return a
        pow = div(n-1, 2)
        return a*square( a^pow )
    end
end
## Rational power ##
^(a::TaylorN,x::Rational) = a^(x.num/x.den)
## Real power ##
function ^(a::TaylorN, x::Real)
    uno = one(a)
    x == zero(x) && return uno
    x == 0.5 && return sqrt(a)
    order = a.order
    # First non-zero coefficient
    # ...
    @assert a.coeffs[1] != zero(a.coeffs[1])
    aux = ( a.coeffs[1] )^x
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    nCoefTot = sizeTable[end]
    coeffs = zeros(T, nCoefTot)
    coeffs[1] = aux
    for k = 1:order
        @inbounds posI, posF = posAuxTable[end][k]
        @inbounds coeffs[posI:posF] = powHomogCoefN(k, v, x, coeffs, 0)
    end
    TaylorN(coeffs,order)
end
function powHomogCoefN{T<:Number}(k::Int, ac::Array{T,1}, x::Real, 
    coeffs::Array{T,1}, knull::Int)
    #
    k == knull && return (ac[ordLHopital+1])^x
    @inbounds posI, posF = posAuxTable[end][k]
    numCoefk = posF - posI + 1
    cc = zeros(T, numCoefk)
    for ka = 0:k-1
        kb = k-ka
        @inbounds posaI, posaF = posAuxTable[end][ka]
        @inbounds posbI, posbF = posAuxTable[end][kb]
        for ja = posaI:posaF
            @inbounds inda = indicesTable[end][ja]
            for jb = posbI:posbF
                @inbounds indb = indicesTable[end][jb]
                @inbounds pos = posTable[end][ inda + indb ]
                pos = pos - posI + 1
                @inbounds cc[pos] += ( x*kb-ka ) * ac[jb] * coeffs[ja]
            end
        end
    end
    @inbounds cc[1:end] = cc[1:end] / ( k * ac[1] )
    cc
end
^{T<:Number,S<:Number}(a::TaylorN{T}, x::Complex{S}) = exp( x*log(a) )
^(a::TaylorN, b::TaylorN) = exp( b*log(a) )

## Square ##
function square{T<:Number}(a::TaylorN{T})
    order = a.order
    nCoefTot = sizeTable[end]
    coeffs = zeros(T, nCoefTot)
    coeffs[1] = a.coeffs[1]^2
    for k = 1:order
        @inbounds posI, posF = posAuxTable[end][k]
        @inbounds coeffs[posI:posF] = squareHomogCoefN( k, a.coeffs )
    end
    TaylorN(coeffs,order)
end
# Homogeneous coefficients for square
function squareHomogCoefN{T<:Number}(k::Int, ac::Array{T,1})
    k == 0 && return (ac[1])^2
    @inbounds posI, posF = posAuxTable[end][k]
    numCoefk = posF - posI + 1
    coeffs = zeros(T, numCoefk)
    kodd = k%2
    kend = div(k - 2 + kodd, 2)
    for ka = 0:kend
        kb = k - ka
        @inbounds posaI, posaF = posAuxTable[end][ka]
        @inbounds posbI, posbF = posAuxTable[end][kb]
        for ja = posaI:posaF
            @inbounds inda = indicesTable[end][ja]
            for jb = posbI:posbF
                @inbounds indb = indicesTable[end][jb]
                @inbounds pos = posTable[end][ inda + indb ]
                pos = pos - posI + 1
                @inbounds coeffs[pos] += ac[ja] * ac[jb]
            end
        end
    end
    @inbounds coeffs[1:end] = 2coeffs[1:end]
    #
    if kodd == 0
        ka = div(k,2)
        @inbounds posaI, posaF = posAuxTable[end][ka]
        for ja = posaI:posaF
            @inbounds inda = indicesTable[end][ja]
            for jb = posaI:posaF
                @inbounds indb = indicesTable[end][jb]
                @inbounds pos = posTable[end][ inda + indb ]
                pos = pos - posI + 1
                @inbounds coeffs[pos] += ac[ja] * ac[jb]
            end
        end
    end
    coeffs
end

## sqrt ##
function sqrt(a::TaylorN)
    order = a.order
    nCoefTot = sizeTable[end]
    aux = sqrt(a.coeffs[1])
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    coeffs = zeros(T, nCoefTot)
    coeffs[1] = aux
    for k = 1:order
        @inbounds posI, posF = posAuxTable[end][k]
        @inbounds coeffs[posI:posF] = sqrtHomogCoefN( k, v, coeffs, 0 )
    end
    TaylorN(coeffs,order)
end
# Homogeneous coefficients for sqrt
function sqrtHomogCoefN{T<:Number}(k::Int, ac::Array{T,1}, coeffs::Array{T,1}, knull::Int)
    k == 0 && return sqrt( ac[1] )
    @inbounds posI, posF = posAuxTable[end][k]
    numCoefk = posF - posI + 1
    cc = zeros(T, numCoefk)
    kodd = k%2
    kend = div(k - 2 + kodd, 2)
    for ka = 1:kend
        kb = k - ka
        @inbounds posaI, posaF = posAuxTable[end][ka]
        @inbounds posbI, posbF = posAuxTable[end][kb]
        for ja = posaI:posaF
            @inbounds inda = indicesTable[end][ja]
            for jb = posbI:posbF
                @inbounds indb = indicesTable[end][jb]
                @inbounds pos = posTable[end][ inda + indb ]
                pos = pos - posI + 1
                @inbounds cc[pos] += coeffs[ja] * coeffs[jb]
            end
        end
    end
    @inbounds cc[1:end] = 2cc[1:end]
    #
    if kodd == 0
        ka = div(k,2)
        @inbounds posaI, posaF = posAuxTable[end][ka]
        for ja = posaI:posaF
            @inbounds inda = indicesTable[end][ja]
            for jb = posaI:posaF
                @inbounds indb = indicesTable[end][jb]
                @inbounds pos = posTable[end][ inda + indb ]
                pos = pos - posI + 1
                @inbounds cc[pos] += coeffs[ja] * coeffs[jb]
            end
        end
    end
    @inbounds cc[1:end] = (ac[posI:posF] - cc[1:end]) / (2coeffs[1])
    cc
end

## exp ##
function exp(a::TaylorN)
    order = a.order
    nCoefTot = sizeTable[end]
    aux = exp(a.coeffs[1])
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    coeffs = zeros(T, nCoefTot)
    coeffs[1] = aux
    for k = 1:order
        @inbounds posI, posF = posAuxTable[end][k]
        @inbounds coeffs[posI:posF] = expHomogCoefN( k, v, coeffs )
    end
    TaylorN(coeffs,order)
end
# Homogeneous coefficients for exp
function expHomogCoefN{T<:Number}(k::Int, ac::Array{T,1}, coeffs::Array{T,1})
    k == 0 && return exp( ac[1] )
    @inbounds posI, posF = posAuxTable[end][k]
    numCoefk = posF - posI + 1
    cc = zeros(T, numCoefk)
    for ka = 0:k-1
        kb = k - ka
        @inbounds posaI, posaF = posAuxTable[end][ka]
        @inbounds posbI, posbF = posAuxTable[end][kb]
        for ja = posaI:posaF
            @inbounds inda = indicesTable[end][ja]
            for jb = posbI:posbF
                @inbounds indb = indicesTable[end][jb]
                @inbounds pos = posTable[end][ inda + indb ]
                pos = pos - posI + 1
                @inbounds cc[pos] += kb * coeffs[ja] * ac[jb]
            end
        end
    end
    @inbounds cc[1:end] = cc[1:end] / k
    cc
end

## log ##
function log(a::TaylorN)
    order = a.order
    nCoefTot = sizeTable[end]
    aux = log(a.coeffs[1])
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    coeffs = zeros(T, nCoefTot)
    coeffs[1] = aux
    for k = 1:order
        @inbounds posI, posF = posAuxTable[end][k]
        @inbounds coeffs[posI:posF] = logHomogCoefN( k, v, coeffs )
    end
    TaylorN(coeffs, order)
end
# Homogeneous coefficients for square
function logHomogCoefN{T<:Number}(k::Int, ac::Array{T,1}, coeffs::Array{T,1})
    k == 0 && return log( ac[1] )
    @inbounds posI, posF = posAuxTable[end][k]
    numCoefk = posF - posI + 1
    cc = zeros(T, numCoefk)
    for ka = 1:k-1
        kb = k - ka
        @inbounds posaI, posaF = posAuxTable[end][ka]
        @inbounds posbI, posbF = posAuxTable[end][kb]
        for ja = posaI:posaF
            @inbounds inda = indicesTable[end][ja]
            for jb = posbI:posbF
                @inbounds indb = indicesTable[end][jb]
                @inbounds pos = posTable[end][ inda + indb ]
                pos = pos - posI + 1
                @inbounds cc[pos] += ka * coeffs[ja] * ac[jb]
            end
        end
    end
    @inbounds cc[1:end] = ( ac[posI:posF] - cc[1:end] / k ) / ac[1]
    cc
end

## sin and cos ##
sin(a::TaylorN) = sincos(a)[1]
cos(a::TaylorN) = sincos(a)[2]
function sincos(a::TaylorN)
    order = a.order
    nCoefTot = sizeTable[end]
    aux = sin(a.coeffs[1])
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    sincoeffs = zeros(T,nCoefTot)
    coscoeffs = zeros(T,nCoefTot)
    sincoeffs[1] = aux
    coscoeffs[1] = cos( a.coeffs[1] )
    for k = 1:order
        @inbounds posI, posF = posAuxTable[end][k]
        @inbounds sincoeffs[posI:posF], coscoeffs[posI:posF] = 
            sincosHomogCoefN(k, v, sincoeffs, coscoeffs)
    end
    return TaylorN( sincoeffs, order ), TaylorN( coscoeffs, order )
end
# Homogeneous coefficients for sin and cos
function sincosHomogCoefN{T<:Number}(k::Int, ac::Array{T,1}, scoeffs::Array{T,1}, ccoeffs::Array{T,1})
    k == 0 && return sin( ac[1] ), cos( ac[1] )
    @inbounds posI, posF = posAuxTable[end][k]
    numCoefk = posF - posI + 1
    sv = zeros(T, numCoefk)
    cv = zeros(T, numCoefk)
    for ka = 0:k-1
        kb = k - ka
        @inbounds posaI, posaF = posAuxTable[end][ka]
        @inbounds posbI, posbF = posAuxTable[end][kb]
        for ja = posaI:posaF
            @inbounds inda = indicesTable[end][ja]
            for jb = posbI:posbF
                @inbounds indb = indicesTable[end][jb]
                @inbounds pos = posTable[end][ inda + indb ]
                pos = pos - posI + 1
                @inbounds begin
                    x = kb * ac[jb]
                    sv[pos] += x * ccoeffs[ja]
                    cv[pos] -= x * scoeffs[ja]
                end
            end
        end
    end
    @inbounds sv[1:end] = sv[1:end] / k
    @inbounds cv[1:end] = cv[1:end] / k
    sv, cv
end

## Differentiation ##
"""Partial differentiation of a TaylorN series with respect to the r-th variable"""
function diffTaylor{T<:Number}(a::TaylorN{T}, r::Int)
    @assert 1 <= r <= NUMVARS[end]
    order = a.order
    nCoefTot = sizeTable[end]
    coeffs = zeros(T, nCoefTot)
    iIndices = zeros(Int, NUMVARS[end])
    for pos = 1:nCoefTot
        @inbounds iIndices[1:end] = indicesTable[end][pos]
        @inbounds n = iIndices[r]
        n == 0 && continue
        @inbounds iIndices[r] -= 1
        @inbounds posI = posTable[end][ iIndices ]
        @inbounds coeffs[posI] = n * a.coeffs[pos]
    end
    return TaylorN( coeffs, order )
end
diffTaylor(a::TaylorN) = diffTaylor(a, 1)

## TO BE DONE: Integration...

## Evaluates a Taylor polynomial on a given point ##
## NEEDS REVISION since yields only approx results
function evalTaylor{T<:Number,S<:Number}(a::TaylorN{T}, vals::Array{S,1} )
    numVars = NUMVARS[end]
    @assert length(vals) == numVars
    R = promote_type(T,S)
    suma = zero(R)
    nCoefTot = sizeTable[end]
    iIndices = zeros(Int, numVars)
    for pos = nCoefTot:-1:1
        @inbounds iIndices[1:end] = indicesTable[end][pos]
        @inbounds val = a.coeffs[pos]
        val == zero(R) && continue
        for k = 1:numVars
            @inbounds val *= vals[k]^iIndices[k]
        end
        suma += val
    end
    suma
end
evalTaylor{T<:Number}(a::TaylorN{T}) = evalTaylor(a, zeros(T, a.numVars) )
