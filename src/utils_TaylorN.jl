# utils_TaylorN.jl: N-variables Taylor expansions
#
# Last modification: 2014.04.12
#
# Luis Benet & David P. Sanders
# UNAM
#



## Default values for the maximum degree of polynomials (MAX_DEG) and 
##  number of variables considered (NUM_VARS)
const MAX_DEG = [10]
const NUM_VARS = [2]

## Hash table: coeffsTable
"""
`generateCoeffsTable`: generates a dictionary with the Hash Table.
numVars: number of (real?) variables for the polynomial expansions <--> lNiv
maxDeg: maximum degree of polynomials <--> nPart
"""
function generateCoeffsTable()
    numVars::Int = NUM_VARS[end]
    maxDeg::Int = MAX_DEG[end]
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

function coef2indices!(iV::Int, kDeg::Int, iIndices::Array{Int,1}, pos::Int, dict::Dict{Int64,Array{Int64,1}})

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

## Utilities to get/set MAX_DEG and NUM_VARS; they reset the coeffsTable
get_maxDeg() = MAX_DEG[end]
function set_maxDeg(n::Int)
    @assert n >= 0
    info("MAX_DEG changed; `coeffsTable` regenerated.\n")
    MAX_DEG[end] = n
    coeffsTable[end] = generateCoeffsTable()
    n
end
#
get_numVars() = NUM_VARS[end]
function set_numVars(n::Int)
    @assert n > 0
    info("NUM_VARS changed; `coeffsTable` regenerated.\n")
    NUM_VARS[end] = n
    coeffsTable[end] = generateCoeffsTable()
    n
end


## Constructors ##
immutable TaylorN{T<:Number}
   coeffs :: Array{T,1}
   order :: Int
   nVars :: Int
   function TaylorN(coeffs::Array{T,1}, order::Int, nVars::Int)
        @assert order <= MAX_DEG[end]
        nCoefTot = binomial( NUM_VARS[end]+order, order)
        lencoef = length(coeffs)
        v = zeros(T, nCoefTot)
        @inbounds v[1:lencoef] = coeffs[1:lencoef]
        new(v, order, NUM_VARS[end])
   end
end
#
TaylorN{T<:Number}(x::TaylorN{T}, order::Int) = TaylorN{T}(x.coeffs, order, NUM_VARS[end])
TaylorN{T<:Number}(x::TaylorN{T}) = TaylorN{T}(x.coeffs, x.order, NUM_VARS[end])
#
TaylorN{T<:Number}(coeffs::Array{T,1}, order::Int) = TaylorN{T}(coeffs, order, NUM_VARS[end])
TaylorN{T<:Number}(coeffs::Array{T,1}) = TaylorN{T}(coeffs, MAX_DEG[end], NUM_VARS[end])
#
TaylorN{T<:Number}(x::T, order::Int) = TaylorN{T}([x], order, NUM_VARS[end])
TaylorN{T<:Number}(x::T) = TaylorN{T}([x], 0, NUM_VARS[end])

## Type, length ##
eltype{T<:Number}(::TaylorN{T}) = T
length{T<:Number}(a::TaylorN{T}) = length( a.coeffs )

