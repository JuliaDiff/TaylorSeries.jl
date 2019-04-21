# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

## Auxiliary function ##

"""
    resize_coeffs1!{T<Number}(coeffs::Array{T,1}, order::Int)

If the length of `coeffs` is smaller than `order+1`, it resizes
`coeffs` appropriately filling it with zeros.
"""
function resize_coeffs1!(coeffs::Array{T,1}, order::Int) where {T<:Number}
    lencoef = length(coeffs)
    resize!(coeffs, order+1)
    if order > lencoef-1
        z = zero(coeffs[1])
        coeffs[lencoef+1:order+1] .= z
    end
    return nothing
end

"""
    resize_coeffsHP!{T<Number}(coeffs::Array{T,1}, order::Int)

If the length of `coeffs` is smaller than the number of coefficients
correspondinf to `order` (given by `size_table[order+1]`), it resizes
`coeffs` appropriately filling it with zeros.
"""
function resize_coeffsHP!(coeffs::Array{T,1}, order::Int) where {T<:Number}
    lencoef = length( coeffs )
    @inbounds num_coeffs = size_table[order+1]
    @assert order ≤ get_order() && lencoef ≤ num_coeffs
    num_coeffs == lencoef && return nothing
    resize!(coeffs, num_coeffs)
    z = zero(coeffs[1])
    coeffs[lencoef+1:num_coeffs] .= z
    return nothing
end

## Minimum order of an HomogeneousPolynomial compatible with the vector's length
function orderH(coeffs::Array{T,1}) where {T<:Number}
    ord = 0
    ll = length(coeffs)
    for i = 1:get_order()+1
        @inbounds num_coeffs = size_table[i]
        ll ≤ num_coeffs && break
        ord += 1
    end
    return ord
end

## Maximum order of a HomogeneousPolynomial vector; used by TaylorN constructor
function maxorderH(v::Array{HomogeneousPolynomial{T},1}) where {T<:Number}
    m = 0
    @inbounds for i in eachindex(v)
        m = max(m, v[i].order)
    end
    return m
end



## getcoeff ##
"""
    getcoeff(a, n)

Return the coefficient of order `n::Int` of a `a::Taylor1` polynomial.
"""
getcoeff(a::Taylor1, n::Int) = (@assert 0 ≤ n ≤ a.order; return a[n])

getindex(a::Taylor1, n::Int) = a.coeffs[n+1]
getindex(a::Taylor1, u::UnitRange) = view(a.coeffs, u .+ 1 )
getindex(a::Taylor1, c::Colon) = view(a.coeffs, c)

setindex!(a::Taylor1{T}, x::T, n::Int) where {T<:Number} = a.coeffs[n+1] = x
setindex!(a::Taylor1{T}, x::T, u::UnitRange) where {T<:Number} =
    a.coeffs[u .+ 1] .= x
function setindex!(a::Taylor1{T}, x::Array{T,1}, u::UnitRange) where {T<:Number}
    # a.coeffs[u .+ 1] .= x
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a.coeffs[u[ind]+1] = x[ind]
    end
end
setindex!(a::Taylor1{T}, x::T, c::Colon) where {T<:Number} = a.coeffs[c] .= x
setindex!(a::Taylor1{T}, x::Array{T,1}, c::Colon) where {T<:Number} = a.coeffs[c] .= x


"""
    getcoeff(a, v)

Return the coefficient of `a::HomogeneousPolynomial`, specified by `v`,
which is a tuple (or vector) with the indices of the specific
monomial.
"""
function getcoeff(a::HomogeneousPolynomial, v::NTuple{N,Int}) where {N}
    @assert N == get_numvars() && all(v .>= 0)
    kdic = in_base(get_order(),v)
    @inbounds n = pos_table[a.order+1][kdic]
    a[n]
end
getcoeff(a::HomogeneousPolynomial, v::Array{Int,1}) = getcoeff(a, (v...,))

getindex(a::HomogeneousPolynomial, n::Int) = a.coeffs[n]
getindex(a::HomogeneousPolynomial, n::UnitRange) = view(a.coeffs, n)
getindex(a::HomogeneousPolynomial, c::Colon) = view(a.coeffs, c)

setindex!(a::HomogeneousPolynomial{T}, x::T, n::Int) where {T<:Number} =
    a.coeffs[n] = x
setindex!(a::HomogeneousPolynomial{T}, x::T, n::UnitRange) where {T<:Number} =
    a.coeffs[n] .= x
setindex!(a::HomogeneousPolynomial{T}, x::Array{T,1}, n::UnitRange) where {T<:Number} =
    a.coeffs[n] .= x
setindex!(a::HomogeneousPolynomial{T}, x::T, c::Colon) where {T<:Number} =
    a.coeffs[c] .= x
setindex!(a::HomogeneousPolynomial{T}, x::Array{T,1}, c::Colon) where {T<:Number} =
    a.coeffs[c] = x


"""
    getcoeff(a, v)

Return the coefficient of `a::TaylorN`, specified by `v`,
which is a tuple (or vector) with the indices of the specific
monomial.
"""
function getcoeff(a::TaylorN, v::NTuple{N,Int}) where {N}
    order = sum(v)
    @assert order ≤ a.order
    getcoeff(a[order], v)
end
getcoeff(a::TaylorN, v::Array{Int,1}) = getcoeff(a, (v...,))

getindex(a::TaylorN, n::Int) = a.coeffs[n+1]
getindex(a::TaylorN, u::UnitRange) = view(a.coeffs, u .+ 1)
getindex(a::TaylorN, c::Colon) = view(a.coeffs, c)

function setindex!(a::TaylorN{T}, x::HomogeneousPolynomial{T}, n::Int) where
        {T<:Number}
    @assert x.order == n
    a.coeffs[n+1] = x
end
setindex!(a::TaylorN{T}, x::T, n::Int) where {T<:Number} =
    a.coeffs[n+1] = HomogeneousPolynomial(x, n)
function setindex!(a::TaylorN{T}, x::T, u::UnitRange) where {T<:Number}
    for ind in u
        a[ind] = x
    end
    a[u]
end
function setindex!(a::TaylorN{T}, x::Array{HomogeneousPolynomial{T},1}, u::UnitRange) where {T<:Number}
    # a[u[:]] .= x[:]
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a[u[ind]] = x[ind]
    end
end
function setindex!(a::TaylorN{T}, x::Array{T,1}, u::UnitRange) where {T<:Number}
    # a[u] .= x
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a[u[ind]] = x[ind]
    end
end
setindex!(a::TaylorN{T}, x::T, ::Colon) where {T<:Number} =
    (a[0:end] = x; a[:])
setindex!(a::TaylorN{T}, x::Array{HomogeneousPolynomial{T},1}, ::Colon) where
    {T<:Number} = (a[0:end] = x; a[:])
setindex!(a::TaylorN{T}, x::Array{T,1}, ::Colon) where {T<:Number} =
    (a[0:end] = x; a[:])


## eltype, length, get_order ##
for T in (:Taylor1, :TaylorN)
    @eval begin
        @inline iterate(a::$T, state=0) = state > a.order ? nothing : (a.coeffs[state+1], state+1)
        @inline eachindex(a::$T) = firstindex(a):lastindex(a)
        # Base.iterate(rS::Iterators.Reverse{$T}, state=rS.itr.order) = state < 0 ? nothing : (a.coeffs[state], state-1)
        @inline eltype(::$T{S}) where {S<:Number} = S
        @inline length(a::$T) = length(a.coeffs)
        @inline size(a::$T) = (length(a),)
        @inline firstindex(a::$T) = 0
        @inline lastindex(a::$T) = a.order
        @inline get_order(a::$T) = a.order
        @inline axes(a::$T) = axes(a.coeffs)
    end
end

# Base.iterate(a::HomogeneousPolynomial, state=1) = state > a.order, nothing : (a.coeffs[state+1], state+1)
# Base.iterate(rS::Iterators.Reverse{$T}, state=rS.itr.order) = state < 0 ? nothing : (a.coeffs[state], state-1)
# Base.eachindex(a::HomogeneousPolynomial) = firstindex(a):lastindex(a)
@inline eltype(::HomogeneousPolynomial{S}) where {S<:Number} = S
@inline length(a::HomogeneousPolynomial) = size_table[a.order+1]#length(a.coeffs)
@inline size(a::HomogeneousPolynomial) = (length(a),)
@inline firstindex(a::HomogeneousPolynomial) = 1
@inline lastindex(a::HomogeneousPolynomial) = length(a)
@inline get_order(a::HomogeneousPolynomial) = a.order
@inline axes(a::HomogeneousPolynomial) = axes(a.coeffs)


## fixorder ##
for T in (:Taylor1, :TaylorN)
    @eval begin
        @inline function fixorder(a::$T, b::$T)
            a.order == b.order && return a, b
            a.order < b.order &&
                return $T(copy(a.coeffs), b.order), $T(copy(b.coeffs), b.order)
            return $T(copy(a.coeffs), a.order), $T(copy(b.coeffs), a.order)
        end
    end
end

function fixorder(a::HomogeneousPolynomial, b::HomogeneousPolynomial)
    @assert a.order == b.order
    return a, b
end


# Finds the first non zero entry; extended to Taylor1
function Base.findfirst(a::Taylor1{T}) where {T<:Number}
    first = findfirst(x->!iszero(x), a.coeffs)
    isa(first, Nothing) && return -1
    return first-1
end
# Finds the last non-zero entry; extended to Taylor1
function Base.findlast(a::Taylor1{T}) where {T<:Number}
    last = findlast(x->!iszero(x), a.coeffs)
    isa(last, Nothing) && return -1
    return last-1
end


## similar ##
similar(a::Taylor1) = Taylor1(similar(a.coeffs), a.order)
function similar(a::Array{Taylor1{T},1}) where {T}
    ret = Vector{Taylor1{T}}(undef, size(a,1))
    a1 = a[1].coeffs
    fill!(ret, similar(a1))
    return ret
end



"""
    constant_term(a)

Return the constant value (zero order coefficient) for `Taylor1`
and `TaylorN`. The fallback behavior is to return `a` itself if
`a::Number`, or `a[1]` when `a::Vector`.
"""
constant_term(a::Taylor1) = a[0]

constant_term(a::TaylorN) = a[0][1]

constant_term(a::Vector{T}) where {T<:Number}= a[1]

constant_term(a::Number) = a

"""
    linear_polynomial(a)

Return the linear part of `a` as a polynomial (`Taylor1`
or `TaylorN`). The fallback behavior is to return `a` itself if
`a::Number`, or `a[1]` when `a::Vector`.
"""
linear_polynomial(a::Taylor1) = Taylor1([zero(a[1]), a[1]])

linear_polynomial(a::TaylorN) = TaylorN([a[1]])

linear_polynomial(a::Vector{T}) where {T<:Number} = a[1]

linear_polynomial(a::Number) = a
