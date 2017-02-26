# This file is part of the Taylor1Series.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#


## Constructors ##

## Taylor1
doc"""
    Taylor1{T<:Number} <: Number

DataType for polynomial expansions in one independent variable.

**Fields:**

- `coeffs :: Array{T,1}` Expansion coefficients; the $i$-th
    component is the coefficient of degree $i-1$ of the expansion.
- `order  :: Int64` Maximum order (degree) of the polynomial.
"""
immutable Taylor1{T<:Number} <: Number
    coeffs :: Array{T,1}
    order :: Int

    ## Inner constructor ##
    function Taylor1(coeffs::Array{T,1}, order::Int)
        lencoef = length(coeffs)
        order = max(order, lencoef-1)
        if order == lencoef-1
            return new( coeffs, order)
        else
            resize!(coeffs, order+1)
            @inbounds for i = lencoef+1:order+1
                coeffs[i] = zero(T)
            end
            return new( coeffs, order)
        end
    end
end

## Outer constructors ##
Taylor1{T<:Number}(x::Taylor1{T}, order::Int) = Taylor1{T}(x.coeffs, order)
Taylor1{T<:Number}(x::Taylor1{T}) = x
Taylor1{T<:Number}(coeffs::Array{T,1}, order::Int) = Taylor1{T}(coeffs, order)
Taylor1{T<:Number}(coeffs::Array{T,1}) = Taylor1{T}(coeffs, length(coeffs)-1)
Taylor1{T<:Number}(x::T, order::Int) = Taylor1{T}([x], order)

# Shortcut to define Taylor1 independent variables
doc"""
    Taylor1([T::Type=Float64], [order::Int=1])

Shortcut to define the independent variable of a `Taylor1{T}` polynomial of
given `order`. The default type for `T` is `Float64`.

```jldoctest
julia> Taylor1(16)
 1.0 t + 𝒪(t¹⁷)

julia> Taylor1(Rational{Int},4)
 1//1 t + 𝒪(t⁵)

```
"""
Taylor1{T<:Number}(::Type{T}, order::Int=1) = Taylor1{T}( [zero(T), one(T)], order)
Taylor1(order::Int=1) = Taylor1(Float64, order)



## HomogeneousPolynomial (homogeneous polynomial) constructors ##
doc"""
    HomogeneousPolynomial{T<:Number} <: Number

DataType for homogenous polynomials in many (>1) independent variables.

**Fields:**

- `coeffs  :: Array{T,1}` Expansion coefficients of the homogeneous
polynomial; the $i$-th component is related to a monomial, where the degrees
of the independent variables are specified by `coeff_table[order+1][i]`.
- `order   :: Int` order (degree) of the homogenous polynomial.
"""
immutable HomogeneousPolynomial{T<:Number} <: Number
    coeffs  :: Array{T,1}
    order   :: Int

    function HomogeneousPolynomial( coeffs::Array{T,1}, order::Int )
        @assert order <= get_order()
        lencoef = length( coeffs )
        @inbounds num_coeffs = size_table[order+1]
        @assert lencoef <= num_coeffs
        num_coeffs == lencoef && return new(coeffs, order)
        resize!(coeffs, num_coeffs)
        @simd for i = lencoef+1:num_coeffs
            @inbounds coeffs[i] = zero(coeffs[1])
        end
        new(coeffs, order)
    end
end

HomogeneousPolynomial{T<:Number}(x::HomogeneousPolynomial{T}, order::Int) =
    HomogeneousPolynomial{T}(x.coeffs, order)
HomogeneousPolynomial{T<:Number}(x::HomogeneousPolynomial{T}) = x
HomogeneousPolynomial{T<:Number}(coeffs::Array{T,1}, order::Int) =
    HomogeneousPolynomial{T}(coeffs, order)
HomogeneousPolynomial{T<:Number}(coeffs::Array{T,1}) =
    HomogeneousPolynomial{T}(coeffs, orderH(coeffs))
HomogeneousPolynomial{T<:Number}(x::T, order::Int) =
    HomogeneousPolynomial{T}([x], order)
HomogeneousPolynomial{T<:Number}(x::T) = HomogeneousPolynomial{T}([x], 0)



doc"""
    TaylorN{T<:Number} <: Number

DataType for polynomial expansions in many (>1) independent variables.

**Fields:**

- `coeffs  :: Array{HomogeneousPolynomial{T},1}` Vector containing the
`HomogeneousPolynomial` entries. The $i$-th component corresponds to the
homogeneous polynomial of degree $i-1$.
- `order   :: Int`  maximum order of the polynomial expansion.
"""
immutable TaylorN{T<:Number} <: Number
    coeffs  :: Array{HomogeneousPolynomial{T},1}
    order   :: Int

    function TaylorN( v::Array{HomogeneousPolynomial{T},1}, order::Int )
        m = maxorderH(v)
        order = max( m, order )
        coeffs = zeros(HomogeneousPolynomial{T}, order)
        @simd for i in eachindex(v)
            @inbounds ord = v[i].order
            @inbounds coeffs[ord+1] += v[i]
        end
        new(coeffs, order)
    end
end

TaylorN{T<:Number}(x::TaylorN{T}, order::Int) = TaylorN{T}(x.coeffs, order)
TaylorN{T<:Number}(x::TaylorN{T}) = x
TaylorN{T<:Number}(x::Array{HomogeneousPolynomial{T},1}, order::Int) =
    TaylorN{T}(x, order)
TaylorN{T<:Number}(x::Array{HomogeneousPolynomial{T},1}) = TaylorN{T}(x,0)
TaylorN{T<:Number}(x::HomogeneousPolynomial{T},order::Int) =
    TaylorN{T}([x], order)
TaylorN{T<:Number}(x::HomogeneousPolynomial{T}) = TaylorN{T}([x], x.order)
TaylorN{T<:Number}(x::T, order::Int) =
    TaylorN{T}([HomogeneousPolynomial(x)],order)
# TaylorN{T<:Number}(x::T) = TaylorN{T}([HomogeneousPolynomial(x)], 0)

## Shortcut to define TaylorN independent variables
"""
    TaylorN([T::Type=Float64], nv::Int; [order::Int=get_order()])

Shortcut to define the `nv`-th independent `TaylorN{T}` variable as a
polynomial. The order is defined through the keyword parameter `order`,
whose default corresponds to `get_order()`. The default of `T::Type`
`Float64`.
"""
function TaylorN{T<:Number}(::Type{T}, nv::Int; order::Int=get_order())
    @assert 0 < nv <= get_numvars()
    v = zeros(T, get_numvars())
    @inbounds v[nv] = one(T)
    return TaylorN( HomogeneousPolynomial(v,1), order )
end
TaylorN(nv::Int; order::Int=get_order()) = TaylorN(Float64, nv; order=order)





## NCoeffType
# """
#     NCoeffType{T}
#
# Non exported trait
# """
# immutable NCoeffType{T} end
