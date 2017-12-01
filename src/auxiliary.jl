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
        @__dot__ coeffs[lencoef+1:order+1] = z
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
    @__dot__ coeffs[lencoef+1:num_coeffs] = z
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

function getindex(a::Taylor1, n::Int)
    @assert 0 ≤ n ≤ length(a.coeffs)
    nn = n == length(a.coeffs) ? n : n+1
    return a.coeffs[nn]
end
function getindex(a::Taylor1, u::UnitRange)
    u_stop = u.stop == length(a.coeffs) ? u.stop : u.stop+1
    view(a.coeffs, (1+u.start):u_stop )
end
getindex(a::Taylor1, c::Colon) = view(a.coeffs, c)

function setindex!(a::Taylor1{T}, x::T, n::Int) where {T<:Number}
    @assert 0 ≤ n ≤ a.order
    a.coeffs[n+1] = x
end
function setindex!(a::Taylor1{T}, x::T, u::UnitRange) where {T<:Number}
    u_stop = u.stop == length(a.coeffs) ? u.stop : u.stop+1
    a.coeffs[(u.start+1):u_stop] = x
end
function setindex!(a::Taylor1{T}, x::Array{T,1}, u::UnitRange) where {T<:Number}
    u_stop = u.stop == length(a.coeffs) ? u.stop : u.stop+1
    a.coeffs[(u.start+1):u_stop] .= x
end
setindex!(a::Taylor1{T}, x::T, c::Colon) where {T<:Number} = a.coeffs[c] = x
setindex!(a::Taylor1{T}, x::Array{T,1}, c::Colon) where {T<:Number} = a.coeffs[c] = x


"""
    getcoeff(a, v)

Return the coefficient of `a::HomogeneousPolynomial`, specified by
`v::Array{Int,1}` which has the indices of the specific monomial.
"""
function getcoeff(a::HomogeneousPolynomial, v::Array{Int,1})
    @assert length(v) == get_numvars()
    kdic = in_base(get_order(),v)
    @inbounds n = pos_table[a.order+1][kdic]
    a[n]
end

getindex(a::HomogeneousPolynomial, n::Int) = a.coeffs[n]
getindex(a::HomogeneousPolynomial, n::UnitRange) = view(a.coeffs, n)
getindex(a::HomogeneousPolynomial, c::Colon) = view(a.coeffs, c)

setindex!(a::HomogeneousPolynomial{T}, x::T, n::Int) where {T<:Number} =
    a.coeffs[n] = x
setindex!(a::HomogeneousPolynomial{T}, x::T, n::UnitRange) where {T<:Number} =
    a.coeffs[n] = x
setindex!(a::HomogeneousPolynomial{T}, x::Array{T,1}, n::UnitRange) where {T<:Number} =
    a.coeffs[n] = x
setindex!(a::HomogeneousPolynomial{T}, x::T, c::Colon) where {T<:Number} =
    a.coeffs[c] = x
setindex!(a::HomogeneousPolynomial{T}, x::Array{T,1}, c::Colon) where {T<:Number} =
    a.coeffs[c] = x


"""
    getcoeff(a, v)

Return the coefficient of `a::TaylorN`, specified by
`v::Array{Int,1}` which has the indices of the specific monomial.
"""
function getcoeff(a::TaylorN, v::Array{Int,1})
    order = sum(v)
    @assert order ≤ a.order
    getcoeff(a[order], v)
end

function getindex(a::TaylorN, n::Int)
    # @assert 0 ≤ n ≤ length(a.coeffs)
    nn = n == length(a.coeffs) ? n : n+1
    return a.coeffs[nn]
end
function getindex(a::TaylorN, u::UnitRange)
    u_stop = u.stop == length(a.coeffs) ? u.stop : u.stop+1
    view(a.coeffs, (1+u.start):u_stop )
end
getindex(a::TaylorN, c::Colon) = view(a.coeffs, c)

function setindex!(a::TaylorN{T}, x::HomogeneousPolynomial{T}, n::Int) where
        {T<:Number}
    @assert x.order == n
    a.coeffs[n+1] = x
end
setindex!(a::TaylorN{T}, x::T, n::Int) where {T<:Number} =
    a.coeffs[n+1] = HomogeneousPolynomial(x, n)
# function setindex!(a::TaylorN{T}, x::HomogeneousPolynomial{T}, u::UnitRange) where {T<:Number}
#     u_stop = u.stop == length(a.coeffs) ? u.stop : u.stop+1
#     a.coeffs[(u.start+1):u_stop] = x
# end
function setindex!(a::TaylorN{T}, x::T, u::UnitRange) where {T<:Number}
    for ind in u
        a[ind] = x
    end
    a[u]
end
setindex!(a::TaylorN{T}, x::Array{HomogeneousPolynomial{T},1},
    u::UnitRange) where {T<:Number} = a[u[:]] .= x[:]
setindex!(a::TaylorN{T}, x::Array{T,1}, u::UnitRange) where {T<:Number} =
    a[u] .= x
# setindex!(a::TaylorN{T}, x::HomogeneousPolynomial{T}, c::Colon) where {T<:Number} =
#     a.coeffs[c] = x
setindex!(a::TaylorN{T}, x::T, ::Colon) where {T<:Number} =
    (a[0:end] = x; a[:])
setindex!(a::TaylorN{T}, x::Array{HomogeneousPolynomial{T},1}, ::Colon) where
    {T<:Number} = (a[0:end] = x; a[:])
setindex!(a::TaylorN{T}, x::Array{T,1}, ::Colon) where {T<:Number} =
    (a[0:end] = x; a[:])


## eltype, length, endof, get_order ##
for T in (:Taylor1, :TaylorN)
    @eval begin
        eltype(::$T{S}) where {S<:Number} = S
        length(a::$T) = length(a.coeffs)
        endof(a::$T) = a.order
        get_order(a::$T) = a.order

        start(a::$T) = start(a.coeffs)-1
        next(a::$T, ord) = ($T(a[ord], ord), ord+1)
        done(a::$T, ord) = ord > a.order
    end
end

eltype(::HomogeneousPolynomial{S}) where {S<:Number} = S
length(a::HomogeneousPolynomial) = length(a.coeffs)
endof(a::HomogeneousPolynomial) = length(a.coeffs)
get_order(a::HomogeneousPolynomial) = a.order

start(a::HomogeneousPolynomial) = start(a.coeffs)
next(a::HomogeneousPolynomial, ord) = (HomogeneousPolynomial(a[ord]), ord+1)
done(a::HomogeneousPolynomial, ord) = ord > length(a.coeffs)


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


"""
    zero_korder(a)

For `a::Taylor1` returns `zero(a[1])` while for `a::TaylorN` returns
a zero of a k-th order `HomogeneousPolynomial` of proper type.
"""
zero_korder(a::Taylor1, ::Int) = zero(a[0])

zero_korder(a::TaylorN, k::Int) = HomogeneousPolynomial(zero(constant_term(a)), k)


# Finds the first non zero entry; extended to Taylor1
Base.findfirst(a::Taylor1{T}) where {T<:Number} = findfirst(a.coeffs)-1


"""
    order_posTb(order, nv, ord)

Return a vector with the positions, in a `HomogeneousPolynomial` of
order `order`, where the variable `nv` has order `ord`.
"""
function order_posTb(order::Int, nv::Int, ord::Int)
    @assert order ≤ get_order()
    @inbounds indTb = coeff_table[order+1]
    @inbounds num_coeffs = size_table[order+1]
    posV = Int[]
    for pos = 1:num_coeffs
        @inbounds indTb[pos][nv] != ord && continue
        push!(posV, pos)
    end
    posV
end


"""
    constant_term(a)

Return the constant value (zero order coefficient) for `Taylor1`
and `TaylorN`.
"""
constant_term(a::Taylor1) = a[0]

constant_term(a::TaylorN) = a[0][1]

constant_term(a::Vector{T}) where {T<:Number}= a[1]

constant_term(a::Number) = a
