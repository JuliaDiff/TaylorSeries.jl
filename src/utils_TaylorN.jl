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
   numVars :: Int
   function TaylorN(coeffs::Array{T,1}, order::Int, numVars::Int)
        @assert order <= MAX_DEG[end]
        nCoefTot = binomial( NUM_VARS[end]+order, order)
        lencoef = length(coeffs)
        v = zeros(T, nCoefTot)
        @inbounds v[1:lencoef] = coeffs[1:lencoef]
        new(v, order, NUM_VARS[end])
   end
end
TaylorN{T<:Number}(x::TaylorN{T}, order::Int) = TaylorN{T}(x.coeffs, order, NUM_VARS[end])
TaylorN{T<:Number}(x::TaylorN{T}) = TaylorN{T}(x.coeffs, x.order, NUM_VARS[end])
#
TaylorN{T<:Number}(coeffs::Array{T,1}, order::Int) = TaylorN{T}(coeffs, order, NUM_VARS[end])
TaylorN{T<:Number}(coeffs::Array{T,1}) = TaylorN{T}(coeffs, MAX_DEG[end], NUM_VARS[end])
#
TaylorN{T<:Number}(x::T, order::Int) = TaylorN{T}([x], order, NUM_VARS[end])
TaylorN{T<:Number}(x::T) = TaylorN{T}([x], 0, NUM_VARS[end])

## Functions to obtain the number of homogenous coefficients of given degree
"""Returns the number of homogeneous coefficients of degree k for NUM_VARS"""
function numHomogCoefN(k::Int)
    k == 0 && return 1
    binomial( k + NUM_VARS[end] - 1, k )
end
"""Returns the position (key) of the first homogeneous coefficient of degree k for NUM_VARS"""
function posHomogCoefN(k)
    k == 0 && return 1
    binomial( k + NUM_VARS[end]-1, k-1 ) + 1
end
function getHomogCoefN{T<:Number}(a::TaylorN{T}, k::Int)
    posini = posHomogCoefN(k)
    numCoef = numHomogCoefN(k)
    posfin = posini+numCoef-1
    coeffs = zeros(T, numCoef)
    for ipos=1:numCoef
       @inbounds coeffs[ipos] = a.coeffs[posini+ipos]
    end
    coeffs
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
    numCoefs = posHomogCoefN(order+1) - 1
    coeffs = zeros(T, numCoefs)
    for k = 0:order
        posini = posHomogCoefN(k)
        posfin = posHomogCoefN(k+1)
        #!?Call SumProdSHk( s, a, k, b, order-k )
        coeffs[posini:posfin] = 
            mulHomogCoefN(a1.coeffs[1:numCoefs], k, b1.coeffs[1:numCoefs], order-k )
    end
    TaylorN(coeffs, order)
end
function mulHomogCoefN(ac, ka, ab, kb)
    order = ka+kb
    posa = posHomogCoefN(ka)
    numa = numHomogCoefN(ka)
    posb = posHomogCoefN(kb)
    numb = numHomogCoefN(kb)
    for ia = posa:
#   Subroutine ProductSHk(s, a, b, korder)
#     Call ZeroHk( s, korder )
#     Do ik=0, korder
#        Call SumProdSHk( s, a, ik, b, korder-ik )
#     End Do
#     Return
#   End Subroutine ProductSHk
#   !-------------------------
#   !> \brief Implements the product of two homogeneous polynomials
#   !!   (s|(ka+kb))+=(a|ka)*(b|kb)
#   Subroutine SumProdSHk(s, a, ka, b, kb)
#     posa = POSHOMO(ka,NVAR)
#     ncoefa = NCOEFH(ka,NVAR) -1
#     posb = POSHOMO(kb,NVAR)
#     ncoefb = NCOEFH(kb,NVAR) -1
#     Do ja=posa, posa+ncoefa
#        inda = coef2ind(ja)%space
#        Do jb=posb, posb+ncoefb
#           indb = coef2ind(jb)%space
#           pos = ind2coef( inda+indb )
#           s%coef(pos) = s%coef(pos) + a%coef(ja) * b%coef(jb)
#        End Do
#     End Do
#     Return
#   End Subroutine SumProdSHk
#   !-------------------------




*(a::Taylor, b::Number) = Taylor(b*a.coeffs, a.order)
*(a::Number, b::Taylor) = Taylor(a*b.coeffs, b.order)

