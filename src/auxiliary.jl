# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

## Auxiliary function ##

function check_taylor1_order!{T<:Number}(coeffs::Array{T,1}, order::Int)
    lencoef = length(coeffs)
    order = max(order, lencoef-1)
    order == lencoef-1 && return nothing
    resize!(coeffs, order+1)
    @simd for i = lencoef+1:order+1
        @inbounds coeffs[i] = zero(coeffs[1])
    end
    nothing
end

function check_hpoly_order!{T<:Number}(coeffs::Array{T,1}, order::Int)
    lencoef = length( coeffs )
    @inbounds num_coeffs = size_table[order+1]
    @assert order ≤ get_order() && lencoef ≤ num_coeffs
    num_coeffs == lencoef && return nothing
    resize!(coeffs, num_coeffs)
    @simd for i = lencoef+1:num_coeffs
        @inbounds coeffs[i] = zero(coeffs[1])
    end
    nothing
end

## Minimum order of an HomogeneousPolynomial compatible with the vector's length
function orderH{T}(coeffs::Array{T,1})
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
function maxorderH{T<:Number}(v::Array{HomogeneousPolynomial{T},1})
    m = 0
    @inbounds for i in eachindex(v)
        m = max(m, v[i].order)
    end
    return m
end



## get_coeff ##
"""
    get_coeff(a, n)

Return the coefficient of order `n::Int` of a `a::Taylor1` polynomial.
"""
get_coeff(a::Taylor1, n::Int) = (@assert 0 ≤ n ≤ a.order+1;
    return a.coeffs[n+1])

"""
    get_coeff(a, v)

Return the coefficient of `a::HomogeneousPolynomial`, specified by
`v::Array{Int,1}` which has the indices of the specific monomial.
"""
function get_coeff(a::HomogeneousPolynomial, v::Array{Int,1})
    @assert length(v) == get_numvars()
    kdic = in_base(get_order(),v)
    @inbounds n = pos_table[a.order+1][kdic]
    a.coeffs[n]
end

"""
    get_coeff(a, v)

Return the coefficient of `a::TaylorN`, specified by
`v::Array{Int,1}` which has the indices of the specific monomial.
"""
function get_coeff(a::TaylorN, v::Array{Int,1})
    order = sum(v)
    @assert order ≤ a.order
    get_coeff(a.coeffs[order+1], v)
end


## Type, length ##
eltype{T<:Number}(::Taylor1{T}) = T
length{T<:Number}(a::Taylor1{T}) = a.order

eltype{T<:Number}(::HomogeneousPolynomial{T}) = T
length(a::HomogeneousPolynomial) = length( a.coeffs )
get_order(a::HomogeneousPolynomial) = a.order

eltype{T<:Number}(::TaylorN{T}) = T
length(a::TaylorN) = length( a.coeffs )
get_order(x::TaylorN) = x.order


# Almost equivalent to `findfirst`
## Auxiliary function ##
function firstnonzero{T<:Number}(ac::Vector{T})
    # nonzero::Int = length(ac)
    for i in eachindex(ac)
        if ac[i] != zero(T)
            return i-1
            # nonzero = i-1
            # break
        end
    end
    # nonzero
    return length(ac)
end
firstnonzero{T<:Number}(a::Taylor1{T}) = firstnonzero(a.coeffs)

for T in (:Taylor1, :TaylorN)
    @eval begin
        fixorder{T<:Number}(a::$T{T}, order::Int64) = $T{T}(a.coeffs, order)
        function fixorder{R<:Number}(a::$T{R}, b::$T{R})
            a.order == b.order && return a, b
            a.order < b.order && return $T{R}(a.coeffs, b.order), b
            return a, $T{R}(b.coeffs, a.order)
        end
    end
end

function fixorder(a::HomogeneousPolynomial, b::HomogeneousPolynomial)
    @assert a.order == b.order
    return a, b
end



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
