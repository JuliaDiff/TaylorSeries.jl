# utils_TaylorN.jl: N-variables Taylor expansions through homogeneous polynomials
#
# Last modification: 2014.08.19
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
    @assert n > 0
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
    info(string("MAXORDER is now ", n, "; hash tables resetted.\n"))
    return n
end
#
get_numVars() = NUMVARS[end]
function set_numVars(n::Int)
    @assert n > 0
    n==1 && error( string("Use `Taylor` rather than `TaylorN` for one independent variable.") )
    NUMVARS[end] = n
    maxOrd = MAXORDER[end]
    indicesTable[:], sizeTable[:] = generateIndicesTable()
    posTable[:] = generatePosTable()
    info(string("NUMVARS is now ", n, "; hash tables resetted.\n"))
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
immutable HomogPol{T<:Number} <: AbstractSeries{T, NUMVARS[end]}
    coeffs  :: Array{T,1}
    order   :: Int
    function HomogPol( coeffs::Array{T,1}, order::Int )
        @assert order <= MAXORDER[end]
        lencoef = length( coeffs )
        @inbounds nCoefH = sizeTable[order+1]
        @assert lencoef <= nCoefH
        z = zero(T)
        for i=lencoef+1:nCoefH
            push!(coeffs, z)
        end
        new(coeffs, order)
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
    return HomogPol(zeros(T,nCoefH), a.order)
end
function zeros{T<:Number}(::HomogPol{T}, order::Int)
    @assert order <= MAXORDER[end]
    order == 0 && return [HomogPol(zero(T),0)]
    v = HomogPol{T}[]
    sizehint(v, order+1)
    v = push!(v, HomogPol(zero(T),0))
    for ord = 1:order
        z = HomogPol(zero(T),ord)
        push!(v,z)
    end
    return v
end
zeros{T<:Number}(::Type{HomogPol{T}}, order::Int) = zeros( HomogPol(zero(T), 0), order)
function one{T<:Number}(a::HomogPol{T})
    a.order == 0 && HomogPol(one(T), 0)
    nCoefH = sizeTable[a.order+1]
    return HomogPol(ones(T,nCoefH), a.order)
end
function ones{T<:Number}(::HomogPol{T}, order::Int)
    @assert order <= MAXORDER[end]
    order == 0 && return [HomogPol(zero(T),0)]
    v = HomogPol{T}[]
    sizehint(v, order+1)
    v = push!(v,HomogPol(one(T),0))
    for ord = 1:order
        z = HomogPol(one(T),ord)
        push!(v,z)
    end
    return v
end
ones{T<:Number}(::Type{HomogPol{T}}, order::Int) = ones( HomogPol(one(T), 0), order)

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
    function TaylorN( v::Array{HomogPol{T},1}, order::Int )
        ll = length(v)
        @inbounds coeffs = [v[k].order for k=1:ll]; push!(coeffs, order)
        order = maximum(coeffs)
        @assert order <= MAXORDER[end]
        coeffs = zeros(HomogPol{T}, order)
        @simd for i = 1:ll
            @inbounds ord = v[i].order
            @inbounds coeffs[ord+1] += v[i]
        end
        new(coeffs, order)
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

# Fast way to define independent variables
function taylorvar(T::Type, nv::Int, order::Int=MAXORDER[end] )
    @assert (0 < nv <= NUMVARS[end] && order <= MAXORDER[end])
    v = zeros(T, NUMVARS[end])
    @inbounds v[nv] = one(T)
    return TaylorN( HomogPol(v,1), order )
end
taylorvar(nv::Int) = taylorvar(Float64, nv)

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
    TaylorN( convert(Array{HomogPol{T},1}, a.coeffs), a.order )
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
function fixorder(a::TaylorN, order::Int64)
    order <= a.order && return a
    @assert MAXORDER[end] >= order
    T = eltype(a)
    for ord = a.order+1:order
        @inbounds nCoefH = sizeTable[ord+1]
        z = zeros(T, nCoefH)
        zH = HomogPol(z, ord)
        push!(a.coeffs, zH)
    end
    return TaylorN(a.coeffs, order)
end
function fixshape(a::HomogPol, b::HomogPol)
    @assert a.order == b.order
    eltype(a) == eltype(b) && return a, b, a.order
    a, b = promote(a, b)
    return a, b, a.order
end
function fixshape(a::TaylorN, b::TaylorN)
    eltype(a) == eltype(b) && a.order == b.order && return a, b, a.order
    order = a.order < b.order ? b.order : a.order
    if eltype(a) != eltype(b)
        a, b = promote(a, b)
    end
    a.order == b.order && return a, b, a.order
    a, b = fixorder(a, order), fixorder(b, order)
    return a, b, order
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
    a, b, order = fixshape(a, b)
    return a.coeffs == b.coeffs
end

iszero(a::HomogPol) = a == zero(a)

## Addition and substraction ##
for f in (:+, :-), T in (:HomogPol, :TaylorN)
    @eval begin
        function ($f)(a::($T), b::($T))
            a, b, order = fixshape(a,b)
            v = ($f)(a.coeffs, b.coeffs)
            return ($T)(v, order)
        end
        ($f)(a::($T)) = ($T)(($f)(a.coeffs), a.order)
    end
end

## Multiplication ##
function *(a::HomogPol, b::HomogPol)
    T = promote_type( eltype(a), eltype(b) )
    order = a.order + b.order
    order > MAXORDER[end] && return HomogPol(zero(T), MAXORDER[end])
    (iszero(a) || iszero(b)) && return HomogPol(zero(T), order)
    @inbounds nCoefHa = sizeTable[a.order+1]
    @inbounds nCoefHb = sizeTable[b.order+1]
    @inbounds nCoefH  = sizeTable[order+1]
    coeffs = zeros(T, nCoefH)
    z = zero(T)
    for na = 1:nCoefHa
        @inbounds ca = a.coeffs[na]
        ca == z && continue
        @inbounds inda = indicesTable[a.order+1][na]
        for nb = 1:nCoefHb
            @inbounds cb = b.coeffs[nb]
            cb == z && continue
            @inbounds indb = indicesTable[b.order+1][nb]
            @inbounds pos = posTable[order+1][inda+indb]
            @inbounds coeffs[pos] += ca * cb
        end
    end
    return HomogPol(coeffs, order)
end
function *(a::TaylorN, b::TaylorN)
    a, b, order = fixshape(a, b)
    T = eltype(a)
    coeffs = zeros(HomogPol{T}, order)
    @inbounds coeffs[1] = a.coeffs[1] * b.coeffs[1]
    for ord = 1:order
        @inbounds for i = 0:ord
            (iszero(a.coeffs[i+1]) || iszero(b.coeffs[ord-i+1])) && continue
            coeffs[ord+1] += a.coeffs[i+1] * b.coeffs[ord-i+1]
        end
    end
    return TaylorN(coeffs, order)
end

## Division ##
/(a::HomogPol, x::Real) = a*inv(x)
/(a::HomogPol, x::Complex) = a*inv(x)
function /(a::TaylorN, b::TaylorN)
    b0 = b.coeffs[1].coeffs[1]
    @assert b0 != zero(b0)
    a, b, order = fixshape(a, b)
    #!?orddivfact, cdivfact = divfactorization(a, b) # order and coefficient of first factorized term
    invb0 = inv(b0)
    cdivfact = a.coeffs[1] * invb0
    T = eltype(cdivfact)
    coeffs = zeros(HomogPol{T}, order)
    @inbounds coeffs[1] = cdivfact
    for ord = 1:order
        @inbounds for i = 0:ord-1
            (iszero(coeffs[i+1]) || iszero(b.coeffs[ord-i+1])) && continue
            coeffs[ord+1] = coeffs[i+1] * b.coeffs[ord-i+1]
        end
        @inbounds coeffs[ord+1] = (a.coeffs[ord+1] - coeffs[ord+1]) * invb0
    end
    return TaylorN(coeffs, order)
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
            @inbounds y = ($op)(a.coeffs[1].coeffs[1], x)
            a.coeffs[1] = HomogPol(y)
            return TaylorN( a.coeffs, a.order )
        end
    end
end
function mod2pi{T<:Real}(a::TaylorN{T}) 
    @inbounds y = mod2pi(a.coeffs[1].coeffs[1])
    a.coeffs[1] = HomogPol(y)
    return TaylorN( a.coeffs, a.order )
end

## Int power ##
function ^(a::HomogPol, n::Integer)
    @assert n >= 0
    T = eltype(a)
    n == 0 && return HomogPol( one(T) )
    n*a.order > MAXORDER[end] && return HomogPol(zero(T), MAXORDER[end])
    n == 1 && return a
    n == 2 && return square(a)
    pow, rest = divrem(n,2)
    b = square(a)
    b = b^pow
    rest == 0 && return b     # even power
    return a * b              # odd power
end
function ^(a::TaylorN, n::Integer)
    uno = one(eltype(a))
    n < 0 && return uno / a^(-n)
    n == 0 && return TaylorN( uno, a.order )
    n == 1 && return a
    n == 2 && return square(a)
    pow, rest = divrem(n,2)
    b = square(a)
    b = b^pow
    rest == 0 && return b     # even power
    return a * b              # odd power
end
## Rational power ##
^(a::TaylorN, x::Rational) = a^(x.num/x.den)
## Real power ##
function ^(a::TaylorN, x::Real)
    uno = one(eltype(a))
    x == zero(x) && return TaylorN( uno )
    x == 0.5 && return sqrt(a)
    x == int(x) && return a^int(x)
    @inbounds a0 = a.coeffs[1].coeffs[1]
    @assert a0 != zero(a0)
    aux = ( a0 )^x
    T = typeof(aux)
    #order = MAXORDER[end]
    #a = TaylorN(a, order)
    coeffs = zeros(HomogPol{T}, a.order)
    @inbounds coeffs[1] = HomogPol( aux )
    for ord = 1:a.order
        for i = 0:ord-1
            tt = x*(ord-i)-i
            @inbounds cpol = coeffs[i+1]
            @inbounds apol = a.coeffs[ord-i+1]
            (iszero(cpol) || iszero(apol)) && continue
            @inbounds coeffs[ord+1] += tt * cpol * apol
        end
        @inbounds coeffs[ord+1] = coeffs[ord+1] / (ord*a0)
    end
    return TaylorN(coeffs, a.order)
end
^(a::TaylorN, x::Complex) = exp( x*log(a) )
^(a::TaylorN, b::TaylorN) = exp( b*log(a) )

## Square ##
function square(a::HomogPol)
    T = eltype(a)
    order = 2*a.order
    @inbounds order > MAXORDER[end] && return HomogPol(zero(T), MAXORDER[end])
    @inbounds nCoefHa = sizeTable[a.order+1]
    @inbounds nCoefH  = sizeTable[order+1]
    two = 2*one(T)
    coeffs = zeros(T, nCoefH)
    for na = 1:nCoefHa
        @inbounds ca = a.coeffs[na]
        ca == zero(T) && continue
        @inbounds inda = indicesTable[a.order+1][na]
        iind = 2*inda
        @inbounds pos = posTable[order+1][iind]
        @inbounds coeffs[pos] += ca * ca
        for nb = na+1:nCoefHa
            @inbounds cb = a.coeffs[nb]
            cb == zero(T) && continue
            @inbounds indb = indicesTable[a.order+1][nb]
            @inbounds pos = posTable[order+1][inda + indb]
            @inbounds coeffs[pos] += two * ca * cb
        end
    end
    return HomogPol(coeffs, order)
end
function square(a::TaylorN)
    T = eltype(a)
    coeffs = zeros(HomogPol{T}, a.order)
    @inbounds coeffs[1] = square(a.coeffs[1])
    for ord = 1:a.order
        kodd = ord%2
        kord = div(ord-2+kodd, 2)
        for i = 0: kord
            @inbounds coeffs[ord+1] += a.coeffs[i+1] * a.coeffs[ord-i+1]
        end
        @inbounds coeffs[ord+1] = 2 * coeffs[ord+1]
        ord%2 == 1 && continue
        kodd = div(ord,2)
        @inbounds coeffs[ord+1] += square( a.coeffs[kodd+1] )
    end
    return TaylorN(coeffs, a.order)
end

## sqrt ##
function sqrt(a::TaylorN)
    @inbounds p0 = sqrt( a.coeffs[1].coeffs[1] )
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
    return TaylorN(coeffs, a.order)
end

## exp ##
function exp(a::TaylorN)
    order = a.order
    @inbounds aux = exp( a.coeffs[1].coeffs[1] )
    T = typeof(aux)
    coeffs = zeros(HomogPol{T}, order)
    @inbounds coeffs[1] = HomogPol(aux, 0)
    for ord = 1:order
        for j = 0:ord-1
            @inbounds coeffs[ord+1] += (ord-j) * a.coeffs[ord-j+1] * coeffs[j+1]
        end
        @inbounds coeffs[ord+1] = coeffs[ord+1] / ord
    end
    return TaylorN(coeffs, order)
end

## log ##
function log(a::TaylorN)
    order = a.order
    @inbounds a0 = a.coeffs[1].coeffs[1]
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
    return TaylorN(coeffs, order)
end

## sin and cos ##
sin(a::TaylorN) = sincos(a)[1]
cos(a::TaylorN) = sincos(a)[2]
function sincos(a::TaylorN)
    order = a.order
    @inbounds a0 = a.coeffs[1].coeffs[1]
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
    return TaylorN(coeffsSin, order), TaylorN(coeffsCos, order)
end

## tan ##
function tan(a::TaylorN)
    order = a.order
    @inbounds a0 = a.coeffs[1].coeffs[1]
    t0 = tan(a0)
    T = typeof(t0)
    coeffsTan = zeros(HomogPol{T}, order)
    coeffsAux = zeros(HomogPol{T}, order)
    @inbounds coeffsTan = HomogPol(t0)
    @inbounds coeffsAux = HomogPol(t0^2)
    for ord = 1:order
        @inbounds coeffsAux[ord] = coeffsTan[ord] * coeffsTan[ord]
        for j = 0:ord-1
            @inbounds coeffsTan[ord+1] += (ord-j) * coeffsTan[ord-j+1] * coeffsAux[j+1]
        end
        @inbounds coeffsTan[ord+1] = a.coeffs[ord+1] + coeffsTan[ord+1] * inv(ord)
    end
    return TaylorN(coeffsTan, order)
end

## Differentiation ##
"""Partial differentiation of a HomogPol series with respect to the r-th variable"""
function diffTaylor(a::HomogPol, r::Int)
    @assert 1 <= r <= NUMVARS[end]
    T = eltype(a)
    a.order == 0 && return HomogPol(zero(T))
    nCoefH = sizeTable[a.order]
    coeffs = zeros(T,nCoefH)
    jind = zeros(Int, NUMVARS[end])
    @inbounds jind[r] = 1
    @inbounds for i = 1:sizeTable[a.order+1]
        iind = indicesTable[a.order+1][i]
        n = iind[r]
        n == 0 && continue
        pos = posTable[a.order][iind-jind]
        coeffs[pos] = n * a.coeffs[i]
    end
    return HomogPol(coeffs, a.order-1)
end
"""Partial differentiation of a TaylorN series with respect to the r-th variable"""
function diffTaylor(a::TaylorN, r::Int)
    T = eltype(a)
    coeffs = zeros(HomogPol{T}, a.order)
    for ord = 1:a.order
        @inbounds coeffs[ord] = diffTaylor( a.coeffs[ord+1], r)
    end
    return TaylorN( coeffs, a.order )
end
diffTaylor(a::TaylorN) = diffTaylor(a, 1)

## Gradient, jacobian and hessian
function gradient(f::TaylorN)
    T = eltype(f)
    numVars = NUMVARS[end]
    grad = zeros(TaylorN{T}, numVars)
    for nv = 1:numVars
        @inbounds grad[nv] = diffTaylor(f, nv)
    end
    return grad
end
∇(f::TaylorN) = gradient(f)
function jacobian{T<:Number}(vf::Array{TaylorN{T},1})
    numVars = NUMVARS[end]
    @assert length(vf) == numVars
    jac = zeros(T, (numVars,numVars))
    for comp = 1:numVars
        @inbounds jac[:,comp] = vf[comp].coeffs[2].coeffs[1:end]
    end
    return transpose(jac)
end
function jacobian{T<:Number,S<:Number}(vf::Array{TaylorN{T},1}, vals::Array{S,1})
    R = promote_type(T,S)
    numVars = NUMVARS[end]
    @assert length(vf) == numVars == length(vals)
    jac = zeros(R, (numVars,numVars))
    for comp = 1:numVars
        grad = gradient( vf[comp] )
        for nv = 1:numVars
            @inbounds jac[nv,comp] = evalTaylor(grad[nv], vals)
        end
    end
    return transpose(jac)
end
hessian{T<:Number,S<:Number}(f::TaylorN{T}, vals::Array{S,1}) = 
    (R = promote_type(T,S); jacobian( gradient(f), vals::Array{R,1}) )
hessian{T<:Number}(f::TaylorN{T}) = hessian( f, zeros(T, NUMVARS[end]))

# ## TO BE DONE: Integration...

## Evaluates a Taylor polynomial on a given point ##
# NEEDS REVISION since results are not quite precise
function evalHomog{T<:Number,S<:Number}(a::HomogPol{T}, vals::Array{S,1} )
    numVars = NUMVARS[end]
    @assert length(vals) == numVars
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
    return suma
end
function evalTaylor{T<:Number,S<:Number}(a::TaylorN{T}, vals::Array{S,1} )
    @assert length(vals) == NUMVARS[end]
    R = promote_type(T,S)
    suma = zero(R)
    for ord = a.order:-1:0
        @inbounds polH = a.coeffs[ord+1]
        suma += evalHomog( polH, vals )
    end
    return suma
end
evalTaylor{T<:Number}(a::TaylorN{T}) = evalTaylor(a, zeros(T, NUMVARS[end]) )

## show
function show(io::IO, a::TaylorN)
    a == zero(a) && (print(io, string(typeof(a), "( [", a.coeffs[1], "], ", a.order, ")")); return)
    strout = string(typeof(a), "( [")
    for ord = 0:a.order
        pol = a.coeffs[ord+1]
        pol == zero(pol) && continue
        strout = string( strout, pol, ", ")
    end
    strout = string(strout[1:end-2], "], ", a.order,")")
    print(io, strout)
end
