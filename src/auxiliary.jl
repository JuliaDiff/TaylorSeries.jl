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
        @simd for ord in lencoef+1:order+1
            @inbounds coeffs[ord] = z
        end
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
    @simd for ord in lencoef+1:num_coeffs
        @inbounds coeffs[ord] = z
    end
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
getindex(a::Taylor1, u::UnitRange{Int}) = view(a.coeffs, u .+ 1 )
getindex(a::Taylor1, c::Colon) = view(a.coeffs, c)
getindex(a::Taylor1{T}, u::StepRange{Int,Int}) where {T<:Number} =
    view(a.coeffs, u[:] .+ 1)

setindex!(a::Taylor1{T}, x::T, n::Int) where {T<:Number} = a.coeffs[n+1] = x
setindex!(a::Taylor1{T}, x::T, u::UnitRange{Int}) where {T<:Number} =
    a.coeffs[u .+ 1] .= x
function setindex!(a::Taylor1{T}, x::Array{T,1}, u::UnitRange{Int}) where {T<:Number}
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a.coeffs[u[ind]+1] = x[ind]
    end
end
setindex!(a::Taylor1{T}, x::T, c::Colon) where {T<:Number} = a.coeffs[c] .= x
setindex!(a::Taylor1{T}, x::Array{T,1}, c::Colon) where {T<:Number} = a.coeffs[c] .= x
setindex!(a::Taylor1{T}, x::T, u::StepRange{Int,Int}) where {T<:Number} =
    a.coeffs[u[:] .+ 1] .= x
function setindex!(a::Taylor1{T}, x::Array{T,1}, u::StepRange{Int,Int}) where {T<:Number}
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a.coeffs[u[ind]+1] = x[ind]
    end
end


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
getindex(a::HomogeneousPolynomial, n::UnitRange{Int}) = view(a.coeffs, n)
getindex(a::HomogeneousPolynomial, c::Colon) = view(a.coeffs, c)
getindex(a::HomogeneousPolynomial, u::StepRange{Int,Int}) = view(a.coeffs, u[:])

setindex!(a::HomogeneousPolynomial{T}, x::T, n::Int) where {T<:Number} =
    a.coeffs[n] = x
setindex!(a::HomogeneousPolynomial{T}, x::T, n::UnitRange{Int}) where {T<:Number} =
    a.coeffs[n] .= x
setindex!(a::HomogeneousPolynomial{T}, x::Array{T,1}, n::UnitRange{Int}) where {T<:Number} =
    a.coeffs[n] .= x
setindex!(a::HomogeneousPolynomial{T}, x::T, c::Colon) where {T<:Number} =
    a.coeffs[c] .= x
setindex!(a::HomogeneousPolynomial{T}, x::Array{T,1}, c::Colon) where {T<:Number} =
    a.coeffs[c] = x
setindex!(a::HomogeneousPolynomial{T}, x::T, u::StepRange{Int,Int}) where {T<:Number} =
    a.coeffs[u[:]] .= x
setindex!(a::HomogeneousPolynomial{T}, x::Array{T,1}, u::StepRange{Int,Int}) where {T<:Number} =
    a.coeffs[u[:]] .= x[:]


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
getindex(a::TaylorN, u::UnitRange{Int}) = view(a.coeffs, u .+ 1)
getindex(a::TaylorN, c::Colon) = view(a.coeffs, c)
getindex(a::TaylorN, u::StepRange{Int,Int}) = view(a.coeffs, u[:] .+ 1)

function setindex!(a::TaylorN{T}, x::HomogeneousPolynomial{T}, n::Int) where {T<:Number}
    @assert x.order == n
    a.coeffs[n+1] = x
end
setindex!(a::TaylorN{T}, x::T, n::Int) where {T<:Number} =
    a.coeffs[n+1] = HomogeneousPolynomial(x, n)
function setindex!(a::TaylorN{T}, x::T, u::UnitRange{Int}) where {T<:Number}
    for ind in u
        a[ind] = x
    end
    a[u]
end
function setindex!(a::TaylorN{T}, x::Array{HomogeneousPolynomial{T},1}, u::UnitRange{Int}) where {T<:Number}
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a[u[ind]] = x[ind]
    end
end
function setindex!(a::TaylorN{T}, x::Array{T,1}, u::UnitRange{Int}) where {T<:Number}
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
function setindex!(a::TaylorN{T}, x::T, u::StepRange{Int,Int}) where {T<:Number}
    for ind in u
        a[ind] = x
    end
    a[u]
end
function setindex!(a::TaylorN{T}, x::Array{HomogeneousPolynomial{T},1}, u::StepRange{Int,Int}) where {T<:Number}
    # a[u[:]] .= x[:]
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a[u[ind]] = x[ind]
    end
end
function setindex!(a::TaylorN{T}, x::Array{T,1}, u::StepRange{Int,Int}) where {T<:Number}
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a[u[ind]] = x[ind]
    end
end


## eltype, length, get_order, etc ##
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval begin
        if $T == HomogeneousPolynomial
            @inline iterate(a::$T, state=1) = state > length(a) ? nothing : (a.coeffs[state], state+1)
            # Base.iterate(rS::Iterators.Reverse{$T}, state=rS.itr.order) = state < 0 ? nothing : (a.coeffs[state], state-1)
            @inline length(a::$T) = size_table[a.order+1]
            @inline firstindex(a::$T) = 1
            @inline lastindex(a::$T) = length(a)
        else
            @inline iterate(a::$T, state=0) = state > a.order ? nothing : (a.coeffs[state+1], state+1)
            # Base.iterate(rS::Iterators.Reverse{$T}, state=rS.itr.order) = state < 0 ? nothing : (a.coeffs[state], state-1)
            @inline length(a::$T) = length(a.coeffs)
            @inline firstindex(a::$T) = 0
            @inline lastindex(a::$T) = a.order
        end
        @inline eachindex(a::$T) = firstindex(a):lastindex(a)
        @inline numtype(::$T{S}) where {S<:Number} = S
        @inline size(a::$T) = size(a.coeffs)
        @inline get_order(a::$T) = a.order
        @inline axes(a::$T) = ()
    end
end
numtype(a) = eltype(a)

@doc doc"""
    numtype(a::AbstractSeries)

Returns the type of the elements of the coefficients of `a`.
""" numtype

# Dumb methods included to properly export normalize_taylor (if IntervalArithmetic is loaded)
@inline normalize_taylor(a::AbstractSeries) = a


## fixorder ##
for T in (:Taylor1, :TaylorN)
    @eval begin
        @inline function fixorder(a::$T, b::$T)
            a.order == b.order && return a, b
            minorder, maxorder = minmax(a.order, b.order)
            if minorder ≤ 0
                minorder = maxorder
            end
            return $T(copy(a.coeffs), minorder), $T(copy(b.coeffs), minorder)
        end
    end
end

function fixorder(a::HomogeneousPolynomial, b::HomogeneousPolynomial)
    @assert a.order == b.order
    return a, b
end


# Finds the first non-zero entry; extended to Taylor1
function Base.findfirst(a::Taylor1{T}) where {T<:Number}
    first = findfirst(x->!iszero(x), a.coeffs)
    isnothing(first) && return -1
    return first-1
end

# Finds the first non-zero entry; extended to HomogeneousPolynomial and TaylorN
function Base.findfirst(a::HomogeneousPolynomial{T}) where {T<:Number}
    first = findfirst(x->!iszero(x), a.coeffs)
    isnothing(first) && return -1
    return first
end

function Base.findfirst(a::TaylorN{T}) where {T<:Number}
    ord = 1
    first = 1
    for ord in eachindex(a.coeffs)
        first = findfirst(a.coeffs[ord])
        first ≥ 1 && return (ord-1, first)
    end
    return (-1, first)
end

# Finds the last non-zero entry; extended to Taylor1
function Base.findlast(a::Taylor1{T}) where {T<:Number}
    last = findlast(x->!iszero(x), a.coeffs)
    isnothing(last) && return -1
    return last-1
end

# Finds the last non-zero entry; extended to HomogeneousPolynomial and TaylorN
function Base.findlast(a::HomogeneousPolynomial{T}) where {T<:Number}
    last = findlast(x->!iszero(x), a.coeffs)
    isnothing(last) && return -1
    return last
end
function Base.findlast(a::TaylorN{T}) where {T<:Number}
    ord = 1
    last = 1
    for ord in reverse(eachindex(a.coeffs))
        last = findlast(a.coeffs[ord])
        last ≥ 1 && return (ord-1, last)
    end
    return (-1, first)
end



## copyto! ##
# Inspired from base/abstractarray.jl, line 665
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval function copyto!(dst::$T{T}, src::$T{T}) where {T<:Number}
        length(dst) < length(src) && throw(ArgumentError(string("Destination has fewer elements than required; no copy performed")))
        destiter = eachindex(dst)
        y = iterate(destiter)
        for x in src
            dst[y[1]] = x
            y = iterate(destiter, y[2])
        end
        return dst
    end
end


"""
    constant_term(a)

Return the constant value (zero order coefficient) for `Taylor1`
and `TaylorN`. The fallback behavior is to return `a` itself if
`a::Number`, or `a[1]` when `a::Vector`.
"""
constant_term(a::Taylor1) = a[0]

constant_term(a::TaylorN) = a[0][1]

constant_term(a::Vector{T}) where {T<:Number} = constant_term.(a)

constant_term(a::Number) = a

"""
    linear_polynomial(a)

Returns the linear part of `a` as a polynomial (`Taylor1` or `TaylorN`),
*without* the constant term. The fallback behavior is to return `a` itself.
"""
linear_polynomial(a::Taylor1) = Taylor1([zero(a[1]), a[1]], a.order)

linear_polynomial(a::HomogeneousPolynomial) = HomogeneousPolynomial(a[1], a.order)

linear_polynomial(a::TaylorN) = TaylorN(a[1], a.order)

linear_polynomial(a::Vector{T}) where {T<:Number} = linear_polynomial.(a)

linear_polynomial(a::Number) = a

"""
    nonlinear_polynomial(a)

Returns the nonlinear part of `a`. The fallback behavior is to return `zero(a)`.
"""
nonlinear_polynomial(a::AbstractSeries) = a - constant_term(a) - linear_polynomial(a)

nonlinear_polynomial(a::Vector{T}) where {T<:Number} = nonlinear_polynomial.(a)

nonlinear_polynomial(a::Number) = zero(a)
