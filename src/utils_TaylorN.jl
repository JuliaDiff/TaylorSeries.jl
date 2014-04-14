# utils_TaylorN.jl: N-variables Taylor expansions
#
# Last modification: 2014.04.14
#
# Luis Benet & David P. Sanders
# UNAM
#


## Default values for the maximum degree of polynomials (MAXDEG) and 
##  number of variables considered (NUMVARS)
const MAXDEG = [4]
const NUMVARS = [2]

## Hash table: coeffsTable
"""
`generateCoeffsTable`: generates a dictionary with the HashTable.
numVars: number of (real?) variables for the polynomial expansions <--> lNiv
maxDeg: maximum degree of polynomials <--> nPart
"""
function generateCoeffsTable()
    numVars::Int = NUMVARS[end]
    maxDeg::Int = MAXDEG[end]
    DDic = Dict{Int64, Array{Int64,1}}()

    if numVars==1
        info( string("`TaylorSeries.jl` package is more appropriate\n", 
            "(much MUCH faster!) for 1-variable expansions.\n", 
            "Instead of `TaylorN` use `Taylor` constructor.\n") )
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
            @inbounds pos, DDic[pos] = coef2indices!( iV, kDeg, iindices, pos, DDic)
        end
    end
    return DDic
end
function coef2indices!(iV::Int, kDeg::Int, iIndices::Array{Int,1}, pos::Int, 
    dict::Dict{Int64,Array{Int64,1}})

    jVar = iV-1
    kDegNext = kDeg - iIndices[iV]
    if jVar > 1
        for jDeg=0:kDegNext
            iIndices[jVar] = jDeg
            @inbounds pos, dict[pos] = coef2indices!( jVar, kDegNext, iIndices, pos, dict)
        end
    else
        iIndices[1] = kDegNext
        pos += 1
        @inbounds dict[pos] = iIndices[1:end]
    end
    return pos, iIndices[1:end]
end

const coeffsTable = [ generateCoeffsTable() ]
const sizeCoeffsTable = [ length( coeffsTable[end]) ]

"""Returns dictionary-key from indices"""
function indices2coef(iIndices::Array{Int,1})
    for i = 1:sizeCoeffsTable[end]
        iIndices == coeffsTable[end][i] && return i
    end
    error(string("Provided indices ARE NOT in `coeffsTable`"))
end

## Utilities to get/set MAXDEG and NUMVARS; they reset the coeffsTable
get_maxDeg() = MAXDEG[end]
function set_maxDeg(n::Int)
    @assert n >= 0
    info("MAXDEG changed; `coeffsTable` regenerated.\n")
    MAXDEG[end] = n
    coeffsTable[end] = generateCoeffsTable()
    sizeCoeffsTable[end] = length( coeffsTable[end] )
    n
end
#
get_numVars() = NUMVARS[end]
function set_numVars(n::Int)
    @assert n > 0
    info("NUMVARS changed; `coeffsTable` regenerated.\n")
    NUMVARS[end] = n
    coeffsTable[end] = generateCoeffsTable()
    sizeCoeffsTable[end] = length( coeffsTable[end] )
    n
end

## Constructors ##
immutable TaylorN{T<:Number}
   coeffs :: Array{T,1}
   order :: Int
   numVars :: Int
   function TaylorN(coeffs::Array{T,1}, order::Int, numVars::Int)
        @assert order <= MAXDEG[end]
        lencoef = length(coeffs)
        nCoefTot = sizeCoeffsTable[end]
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
TaylorN{T<:Number}(coeffs::Array{T,1}) = TaylorN{T}(coeffs, MAXDEG[end], NUMVARS[end])
#
TaylorN{T<:Number}(x::T, order::Int) = TaylorN{T}([x], order, NUMVARS[end])
TaylorN{T<:Number}(x::T) = TaylorN{T}([x], 0, NUMVARS[end])

## Functions to obtain the number of homogenous coefficients of given degree
"""Returns the number of homogeneous coefficients of degree k for NUMVARS"""
function numHomogCoefN(k::Int)
    k == 0 && return 1
    return binomial( k + NUMVARS[end] - 1, k )
end
"""Returns the position (key) of the first homogeneous coefficient of degree k for NUMVARS"""
function posHomogCoefN(k::Int)
    k == 0 && return 1
    return binomial( k + NUMVARS[end]-1, k-1 ) + 1
end

## Type, length ##
eltype{T<:Number}(::TaylorN{T}) = T
length(a::TaylorN) = length( a.coeffs )
get_numVars(x::TaylorN) = x.numVars

## Conversion and promotion rules ##
convert{T<:Number}(::Type{TaylorN{T}}, a::TaylorN) = TaylorN(convert(Array{T,1}, a.coeffs), a.order)
convert{T<:Number, S<:Number}(::Type{TaylorN{T}}, b::Array{S,1}) = TaylorN(convert(Array{T,1},b))
convert{T<:Number, S<:Number}(::Type{TaylorN{T}}, b::S) = TaylorN([convert(T,b)], 0)
promote_rule{T<:Number, S<:Number}(::Type{TaylorN{T}}, ::Type{TaylorN{S}}) = TaylorN{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{TaylorN{T}}, ::Type{Array{S,1}}) = TaylorN{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{Array{S,1}}, ::Type{TaylorN{T}}) = TaylorN{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{TaylorN{T}}, ::Type{S}) = TaylorN{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{S}, ::Type{TaylorN{T}}) = TaylorN{promote_type(T, S)}

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
    @eval ($f)(a::TaylorN) = TaylorN(($f)(a.coeffs), a.order)
end
ctranspose(a::TaylorN) = conj(a)

## zero and one ##
zero{T<:Number}(a::TaylorN{T}) = TaylorN(zero(T), a.order)
one{T<:Number}(a::TaylorN{T}) = TaylorN(one(T), a.order)

## Equality ##
function ==(a::TaylorN, b::TaylorN)
    a1, b1, order = fixshape(a, b)
    return a1.coeffs == b1.coeffs
end
==(a::TaylorN, b::Number) = ==(a, TaylorN(b, a.order))
==(a::Number, b::TaylorN) = ==(b, TaylorN(a, b.order))

## Addition and substraction ##
for f in (:+, :-)
    @eval begin
        function ($f)(a::TaylorN, b::TaylorN)
            a1, b1, order = fixshape(a, b)
            v = ($f)(a1.coeffs, b1.coeffs)
            return TaylorN(v, order)
        end
        ($f)(a::TaylorN, b::Number) = ($f)(a, TaylorN(b, a.order))
        ($f)(a::Number, b::TaylorN) = ($f)(TaylorN(a, b.order), b)
        ($f)(a::TaylorN) = TaylorN(($f)(a.coeffs), a.order)
    end
end

## Multiplication ##
function *(a::TaylorN, b::TaylorN)
    a1, b1, order = fixshape(a, b)
    T = eltype(a1)
    nCoefTot = sizeCoeffsTable[end]
    coeffs = zeros(T, nCoefTot)
    coeffs[1] = a1.coeffs[1] * b1.coeffs[1]   ## 0th-order
    for k = 1:order
        posI = posHomogCoefN(k)
        posF = posHomogCoefN(k+1)-1
        @inbounds coeffs[posI:posF] = mulHomogCoefN( k, a1.coeffs, b1.coeffs )
    end
    TaylorN(coeffs, order)
end
function mulHomogCoefN{T<:Number}( k::Int, ac::Array{T,1}, bc::Array{T,1} )
    k==0 && return ac[1] * bc[1]
    posI = posHomogCoefN(k)
    numCoefk = numHomogCoefN(k)
    coeffs = zeros(T, numCoefk)
    for ka = 0:k
        kb = k-ka
        posaI = posHomogCoefN(ka)
        posaF = posHomogCoefN(ka+1)-1
        posbI = posHomogCoefN(kb)
        posbF = posHomogCoefN(kb+1)-1
        for ja = posaI:posaF
            @inbounds inda = coeffsTable[end][ja]
            for jb = posbI:posbF
                @inbounds indb = coeffsTable[end][jb]
                pos = indices2coef( inda + indb )
                pos = pos - posI + 1
                @inbounds coeffs[pos] += ac[ja] * bc[jb]
            end
        end
    end
    return coeffs
end
*(a::TaylorN, b::Number) = TaylorN(b*a.coeffs, a.order)
*(a::Number, b::TaylorN) = TaylorN(a*b.coeffs, b.order)

## Division ##
function /(a::TaylorN, b::TaylorN)
    a1, b1, order = fixshape(a, b)
    #!?ordLHopital, cLHopital = divlhopital(a1, b1) # L'Hôpital order and coefficient
    cLHopital = a1.coeffs[1] / b1.coeffs[1]
    T = typeof(cLHopital)
    v1 = convert(Array{T,1}, a1.coeffs)
    v2 = convert(Array{T,1}, b1.coeffs)
    nCoefTot = sizeCoeffsTable[end]
    coeffs = zeros(T, nCoefTot)
    coeffs[1] = cLHopital
    for k = 1:order
        posI = posHomogCoefN(k)
        posF = posHomogCoefN(k+1)-1
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
    posI = posHomogCoefN(k)
    numCoefk = numHomogCoefN(k)
    posF = posI + numCoefk - 1
    cc = zeros(T, numCoefk)
    @inbounds cc = mulHomogCoefN(k, coeffs, bc)
    @inbounds cc[1:end] = (ac[posI:posF]-cc[1:end]) / bc[ordLHopital+1]
    cc
end
/(a::TaylorN,b::Number) = TaylorN(a.coeffs/b, a.order)
/(a::Number,b::TaylorN) = TaylorN([a], b.order) / b

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
## Real power ##
function ^(a::TaylorN, x::Real)
    uno = one(a)
    x == zero(x) && return uno
    #?x == 0.5 && return sqrt(a)
    order = a.order
    # First non-zero coefficient
    # ...
    @assert a.coeffs[1] != zero(a.coeffs[1])
    aux = ( a.coeffs[1] )^x
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    nCoefTot = sizeCoeffsTable[end]
    coeffs = zeros(T, nCoefTot)
    coeffs[1] = aux
    for k = 1:order
        posI = posHomogCoefN(k)
        posF = posHomogCoefN(k+1)-1
        @inbounds coeffs[posI:posF] = powHomogCoefN(k, v, x, coeffs, 0)
    end
    TaylorN(coeffs,order)
end
function powHomogCoefN{T<:Number}(k::Int, ac::Array{T,1}, x::Real, 
    coeffs::Array{T,1}, knull::Int)
    #
    k == knull && return (ac[ordLHopital+1])^x
    posI = posHomogCoefN(k)
    numCoefk = numHomogCoefN(k)
    posF = posI + numCoefk - 1
    cc = zeros(T, numCoefk)
    for ka = 0:k-1
        kb = k-ka
        posaI = posHomogCoefN(ka)
        posaF = posHomogCoefN(ka+1)-1
        posbI = posHomogCoefN(kb)
        posbF = posHomogCoefN(kb+1)-1
        for ja = posaI:posaF
            @inbounds inda = coeffsTable[end][ja]
            for jb = posbI:posbF
                @inbounds indb = coeffsTable[end][jb]
                pos = indices2coef( inda + indb )
                pos = pos - posI + 1
                @inbounds cc[pos] += ( x*kb-ka ) * ac[jb] * coeffs[ja]
            end
        end
    end
    @inbounds cc[1:end] = cc[1:end] / ( k * ac[1] )
    cc
end

## Square ##
function square{T<:Number}(a::TaylorN{T})
    order = a.order
    nCoefTot = sizeCoeffsTable[end]
    coeffs = zeros(T, nCoefTot)
    coeffs[1] = a.coeffs[1]^2
    for k = 1:order
        posI = posHomogCoefN(k)
        posF = posHomogCoefN(k+1)-1
        @inbounds coeffs[posI:posF] = squareHomogCoefN( k, a.coeffs )
    end
    TaylorN(coeffs,order)
end
# Homogeneous coefficients for square
function squareHomogCoefN{T<:Number}(k::Int, ac::Array{T,1})
    k == 0 && return (ac[1])^2
    posI = posHomogCoefN(k)
    numCoefk = numHomogCoefN(k)
    posF = posI + numCoefk - 1
    coeffs = zeros(T, numCoefk)
    kodd = k%2
    kend = div(k - 2 + kodd, 2)
    for ka = 0:kend
        kb = k - ka
        posaI = posHomogCoefN(ka)
        posaF = posHomogCoefN(ka+1)-1
        posbI = posHomogCoefN(kb)
        posbF = posHomogCoefN(kb+1)-1
        for ja = posaI:posaF
            @inbounds inda = coeffsTable[end][ja]
            for jb = posbI:posbF
                @inbounds indb = coeffsTable[end][jb]
                pos = indices2coef( inda + indb )
                pos = pos - posI + 1
                @inbounds coeffs[pos] += ac[ja] * ac[jb]
            end
        end
    end
    @inbounds coeffs[1:end] = 2coeffs[1:end]
    #
    if kodd == 0
        ka = div(k,2)
        posaI = posHomogCoefN(ka)
        posaF = posHomogCoefN(ka+1)-1
        for ja = posaI:posaF
            @inbounds inda = coeffsTable[end][ja]
            for jb = posaI:posaF
                @inbounds indb = coeffsTable[end][jb]
                pos = indices2coef( inda + indb )
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
    nCoefTot = sizeCoeffsTable[end]
    aux = sqrt(a.coeffs[1])
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    coeffs = zeros(T, nCoefTot)
    coeffs[1] = aux
    for k = 1:order
        posI = posHomogCoefN(k)
        posF = posHomogCoefN(k+1)-1
        @inbounds coeffs[posI:posF] = sqrtHomogCoefN( k, v, coeffs, 0 )
    end
    TaylorN(coeffs,order)
end
# Homogeneous coefficients for sqrt
function sqrtHomogCoefN{T<:Number}(k::Int, ac::Array{T,1}, coeffs::Array{T,1}, knull::Int)
    k == 0 && return sqrt( ac[1] )
    posI = posHomogCoefN(k)
    numCoefk = numHomogCoefN(k)
    posF = posI + numCoefk - 1
    cc = zeros(T, numCoefk)
    kodd = k%2
    kend = div(k - 2 + kodd, 2)
    for ka = 1:kend
        kb = k - ka
        posaI = posHomogCoefN(ka)
        posaF = posHomogCoefN(ka+1)-1
        posbI = posHomogCoefN(kb)
        posbF = posHomogCoefN(kb+1)-1
        for ja = posaI:posaF
            @inbounds inda = coeffsTable[end][ja]
            for jb = posbI:posbF
                @inbounds indb = coeffsTable[end][jb]
                pos = indices2coef( inda + indb )
                pos = pos - posI + 1
                @inbounds cc[pos] += coeffs[ja] * coeffs[jb]
            end
        end
    end
    @inbounds cc[1:end] = 2cc[1:end]
    #
    if kodd == 0
        ka = div(k,2)
        posaI = posHomogCoefN(ka)
        posaF = posHomogCoefN(ka+1)-1
        for ja = posaI:posaF
            @inbounds inda = coeffsTable[end][ja]
            for jb = posaI:posaF
                @inbounds indb = coeffsTable[end][jb]
                pos = indices2coef( inda + indb )
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
    nCoefTot = sizeCoeffsTable[end]
    aux = exp(a.coeffs[1])
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    coeffs = zeros(T, nCoefTot)
    coeffs[1] = aux
    for k = 1:order
        posI = posHomogCoefN(k)
        posF = posHomogCoefN(k+1)-1
        @inbounds coeffs[posI:posF] = expHomogCoefN( k, v, coeffs )
    end
    TaylorN(coeffs,order)
end
# Homogeneous coefficients for exp
function expHomogCoefN{T<:Number}(k::Int, ac::Array{T,1}, coeffs::Array{T,1})
    k == 0 && return exp( ac[1] )
    posI = posHomogCoefN(k)
    numCoefk = numHomogCoefN(k)
    posF = posI + numCoefk - 1
    cc = zeros(T, numCoefk)
    for ka = 0:k-1
        kb = k - ka
        posaI = posHomogCoefN(ka)
        posaF = posHomogCoefN(ka+1)-1
        posbI = posHomogCoefN(kb)
        posbF = posHomogCoefN(kb+1)-1
        for ja = posaI:posaF
            @inbounds inda = coeffsTable[end][ja]
            for jb = posbI:posbF
                @inbounds indb = coeffsTable[end][jb]
                pos = indices2coef( inda + indb )
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
    nCoefTot = sizeCoeffsTable[end]
    aux = log(a.coeffs[1])
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    coeffs = zeros(T, nCoefTot)
    coeffs[1] = aux
    for k = 1:order
        posI = posHomogCoefN(k)
        posF = posHomogCoefN(k+1)-1
        @inbounds coeffs[posI:posF] = logHomogCoefN( k, v, coeffs )
    end
    TaylorN(coeffs, order)
end
# Homogeneous coefficients for square
function logHomogCoefN{T<:Number}(k::Int, ac::Array{T,1}, coeffs::Array{T,1})
    k == 0 && return log( ac[1] )
    posI = posHomogCoefN(k)
    numCoefk = numHomogCoefN(k)
    posF = posI + numCoefk - 1
    cc = zeros(T, numCoefk)
    for ka = 1:k-1
        kb = k - ka
        posaI = posHomogCoefN(ka)
        posaF = posHomogCoefN(ka+1)-1
        posbI = posHomogCoefN(kb)
        posbF = posHomogCoefN(kb+1)-1
        for ja = posaI:posaF
            @inbounds inda = coeffsTable[end][ja]
            for jb = posbI:posbF
                @inbounds indb = coeffsTable[end][jb]
                pos = indices2coef( inda + indb )
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
    nCoefTot = sizeCoeffsTable[end]
    aux = sin(a.coeffs[1])
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    sincoeffs = zeros(T,nCoefTot)
    coscoeffs = zeros(T,nCoefTot)
    sincoeffs[1] = aux
    coscoeffs[1] = cos( a.coeffs[1] )
    for k = 1:order
        posI = posHomogCoefN(k)
        posF = posHomogCoefN(k+1)-1
        @inbounds sincoeffs[posI:posF], coscoeffs[posI:posF] = 
            sincosHomogCoefN(k, v, sincoeffs, coscoeffs)
    end
    return TaylorN( sincoeffs, order ), TaylorN( coscoeffs, order )
end
# Homogeneous coefficients for sin and cos
function sincosHomogCoefN{T<:Number}(k::Int, ac::Array{T,1}, scoeffs::Array{T,1}, ccoeffs::Array{T,1})
    k == 0 && return sin( ac[1] ), cos( ac[1] )
    posI = posHomogCoefN(k)
    numCoefk = numHomogCoefN(k)
    posF = posI + numCoefk - 1
    sv = zeros(T, numCoefk)
    cv = zeros(T, numCoefk)
    for ka = 0:k-1
        kb = k - ka
        posaI = posHomogCoefN(ka)
        posaF = posHomogCoefN(ka+1)-1
        posbI = posHomogCoefN(kb)
        posbF = posHomogCoefN(kb+1)-1
        for ja = posaI:posaF
            @inbounds inda = coeffsTable[end][ja]
            for jb = posbI:posbF
                @inbounds indb = coeffsTable[end][jb]
                pos = indices2coef( inda + indb )
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

