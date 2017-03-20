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
function resize_coeffs1!{T<:Number}(coeffs::Array{T,1}, order::Int)
    lencoef = length(coeffs)
    order ≤ lencoef-1 && return nothing
    resize!(coeffs, order+1)
    coeffs[lencoef+1:end] .= zero(coeffs[1])
    return nothing
end

"""
    resize_coeffsHP!{T<Number}(coeffs::Array{T,1}, order::Int)

If the length of `coeffs` is smaller than the number of coefficients
correspondinf to `order` (given by `size_table[order+1]`), it resizes
`coeffs` appropriately filling it with zeros.
"""
function resize_coeffsHP!{T<:Number}(coeffs::Array{T,1}, order::Int)
    lencoef = length( coeffs )
    @inbounds num_coeffs = size_table[order+1]
    @assert order ≤ get_order() && lencoef ≤ num_coeffs
    num_coeffs == lencoef && return nothing
    resize!(coeffs, num_coeffs)
    coeffs[lencoef+1:end] .= zero(coeffs[1])
    return nothing
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
    return a[n+1])

function getindex(a::Taylor1, n::Int)
    (1 ≤ n ≤ length(a.coeffs)) && return a.coeffs[n]
    return zero(a.coeffs[1])
end
getindex(a::Taylor1, n::UnitRange) = a.coeffs[n]
setindex!{T<:Number}(a::Taylor1{T}, x::T, n::Int) = a.coeffs[n] = x
setindex!{T<:Number}(a::Taylor1{T}, x::T, n::UnitRange) = a.coeffs[n] = x



"""
    get_coeff(a, v)

Return the coefficient of `a::HomogeneousPolynomial`, specified by
`v::Array{Int,1}` which has the indices of the specific monomial.
"""
function get_coeff(a::HomogeneousPolynomial, v::Array{Int,1})
    @assert length(v) == get_numvars()
    kdic = in_base(get_order(),v)
    @inbounds n = pos_table[a.order+1][kdic]
    a[n]
end

getindex(a::HomogeneousPolynomial, n::Int) = a.coeffs[n]
getindex(a::HomogeneousPolynomial, n::UnitRange) = a.coeffs[n]
setindex!{T<:Number}(a::HomogeneousPolynomial{T}, x::T, n::Int) =
    a.coeffs[n] = x
setindex!{T<:Number}(a::HomogeneousPolynomial{T}, x::T, n::UnitRange) =
    a.coeffs[n] = x


"""
    get_coeff(a, v)

Return the coefficient of `a::TaylorN`, specified by
`v::Array{Int,1}` which has the indices of the specific monomial.
"""
function get_coeff(a::TaylorN, v::Array{Int,1})
    order = sum(v)
    @assert order ≤ a.order
    get_coeff(a[order+1], v)
end

function getindex(a::TaylorN, n::Int)
    1 ≤ n ≤ length(a.coeffs) && return a.coeffs[n]
    n ≤ get_order()+1 && return zero_korder(a, n-1)
    throw(BoundsError(a.coeffs,n))
end
getindex(a::TaylorN, n::UnitRange) = a.coeffs[n]
setindex!{T<:Number}(a::TaylorN{T}, x::HomogeneousPolynomial{T}, n::Int) =
    a.coeffs[n] = x
setindex!{T<:Number}(a::TaylorN{T}, x::T, n::Int) =
    a.coeffs[n] = HomogeneousPolynomial(x, n-1)
setindex!{T<:Number}(a::TaylorN{T}, x::HomogeneousPolynomial{T}, n::UnitRange) =
    a.coeffs[n] = x
setindex!{T<:Number}(a::TaylorN{T}, x::T, n::UnitRange) =
    a.coeffs[n] = x


## eltype, length, endof, get_order ##
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval begin
        eltype{T<:Number}(::$T{T}) = T

        length(a::$T) = length(a.coeffs)

        endof(a::$T) = length(a.coeffs)

        get_order(a::$T) = a.order
    end
end


"""
    max_order(a)

Returns the maximum order of a `Taylor1` or `TaylorN`.

For `a::Taylor1` returns `a.order` while for `a::TaylorN` returns
`get_order()`.
"""
max_order(a::Taylor1) = a.order

max_order(a::TaylorN) = get_order()


"""
    zero_korder(a)

For `a::Taylor1` returns `zero(a[1])` while for `a::TaylorN` returns
a zero of a k-th order `HomogeneousPolynomial` of proper type.
"""
zero_korder(a::Taylor1, k::Int) = zero(a[1])

zero_korder(a::TaylorN, k::Int) = HomogeneousPolynomial(zero(a[1][1]), k)


# Finds the first non zero entry; extended to Taylor1
Base.findfirst{T<:Number}(a::Taylor1{T}) = findfirst(a.coeffs)-1


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
constant_term(a::Taylor1) = a[1]

constant_term(a::TaylorN) = a[1][1]
