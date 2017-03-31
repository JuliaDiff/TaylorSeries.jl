# This file is part of the Taylor1Series.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

"Abstract type for Taylor1, HomogeneousPolynomial and TaylorN"
@compat abstract type AbstractSeries{T<:Number} <: Number end


## Constructors ##

######################### Taylor1
doc"""
    Taylor1{T<:Number} <: AbstractSeries{T}

DataType for polynomial expansions in one independent variable.

**Fields:**

- `coeffs :: Array{T,1}` Expansion coefficients; the $i$-th
    component is the coefficient of degree $i-1$ of the expansion.
- `order  :: Int64` Maximum order (degree) of the polynomial.
"""
immutable Taylor1{T<:Number} <: AbstractSeries{T}
    coeffs :: Array{T,1}
    order :: Int

    ## Inner constructor ##
    function (::Type{Taylor1{T}}){T}(coeffs::Array{T,1}, order::Int)
        resize_coeffs1!(coeffs, order)
        return new{T}(coeffs, order)
    end
end

## Outer constructors ##
Taylor1{T<:Number}(x::Taylor1{T}) = x
Taylor1{T<:Number}(coeffs::Array{T,1}, order::Int) = Taylor1{T}(coeffs, order)
Taylor1{T<:Number}(coeffs::Array{T,1}) = Taylor1(coeffs, length(coeffs)-1)
function Taylor1{T<:Number}(x::T, order::Int)
    v = zeros(T, order+1)
    v[1] = x
    return Taylor1(v, order)
end

# Shortcut to define Taylor1 independent variables
doc"""
    Taylor1([T::Type=Float64], [order::Int=1])

Shortcut to define the independent variable of a `Taylor1{T}` polynomial of
given `order`. The default type for `T` is `Float64`.

```julia
julia> Taylor1(16)
 1.0 t + 𝒪(t¹⁷)

julia> Taylor1(Rational{Int}, 4)
 1//1 t + 𝒪(t⁵)
```
"""
Taylor1{T<:Number}(::Type{T}, order::Int=1) = Taylor1{T}( [zero(T), one(T)], order)
Taylor1(order::Int=1) = Taylor1(Float64, order)


######################### HomogeneousPolynomial
doc"""
    HomogeneousPolynomial{T<:Number} <: AbstractSeries{T}

DataType for homogenous polynomials in many (>1) independent variables.

**Fields:**

- `coeffs  :: Array{T,1}` Expansion coefficients of the homogeneous
polynomial; the $i$-th component is related to a monomial, where the degrees
of the independent variables are specified by `coeff_table[order+1][i]`.
- `order   :: Int` order (degree) of the homogenous polynomial.
"""
immutable HomogeneousPolynomial{T<:Number} <: AbstractSeries{T}
    coeffs  :: Array{T,1}
    order   :: Int

    function (::Type{HomogeneousPolynomial{T}}){T}(coeffs::Array{T,1}, order::Int)
        resize_coeffsHP!(coeffs, order)
        return new{T}(coeffs, order)
    end
end

HomogeneousPolynomial{T<:Number}(x::HomogeneousPolynomial{T}) = x
HomogeneousPolynomial{T<:Number}(coeffs::Array{T,1}, order::Int) =
    HomogeneousPolynomial{T}(coeffs, order)
HomogeneousPolynomial{T<:Number}(coeffs::Array{T,1}) =
    HomogeneousPolynomial{T}(coeffs, orderH(coeffs))
HomogeneousPolynomial{T<:Number}(x::T, order::Int) =
    HomogeneousPolynomial{T}([x], order)

# Shortcut to define HomogeneousPolynomial independent variable
"""
    HomogeneousPolynomial([T::Type=Float64], nv::Int])

Shortcut to define the `nv`-th independent `HomogeneousPolynomial{T}`.
The default type for `T` is `Float64`.

```julia
julia> HomogeneousPolynomial(1)
 1.0 x₁

julia> HomogeneousPolynomial(Rational{Int}, 2)
 1//1 x₂
```
"""
function HomogeneousPolynomial{T<:Number}(::Type{T}, nv::Int)
    @assert 0 < nv ≤ get_numvars()
    v = zeros(T, get_numvars())
    @inbounds v[nv] = one(T)
    return HomogeneousPolynomial(v, 1)
end
HomogeneousPolynomial(nv::Int) = HomogeneousPolynomial(Float64, nv)



######################### TaylorN
doc"""
    TaylorN{T<:Number} <: AbstractSeries{T}

DataType for polynomial expansions in many (>1) independent variables.

**Fields:**

- `coeffs  :: Array{HomogeneousPolynomial{T},1}` Vector containing the
`HomogeneousPolynomial` entries. The $i$-th component corresponds to the
homogeneous polynomial of degree $i-1$.
- `order   :: Int`  maximum order of the polynomial expansion.
"""
immutable TaylorN{T<:Number} <: AbstractSeries{T}
    coeffs  :: Array{HomogeneousPolynomial{T},1}
    order   :: Int

    function (::Type{TaylorN{T}}){T}(v::Array{HomogeneousPolynomial{T},1}, order::Int)
        m = maxorderH(v)
        order = max( m, order )
        coeffs = zeros(HomogeneousPolynomial{T}, order)
        @simd for i in eachindex(v)
            @inbounds ord = v[i].order
            @inbounds coeffs[ord+1] += v[i]
        end
        new{T}(coeffs, order)
    end
end

TaylorN{T<:Number}(x::TaylorN{T}) = x
TaylorN{T<:Number}(x::Array{HomogeneousPolynomial{T},1}, order::Int) =
    TaylorN{T}(x, order)
TaylorN{T<:Number}(x::Array{HomogeneousPolynomial{T},1}) = TaylorN{T}(x,0)
TaylorN{T<:Number}(x::HomogeneousPolynomial{T},order::Int) =
    TaylorN{T}([x], order)
TaylorN{T<:Number}(x::HomogeneousPolynomial{T}) = TaylorN{T}([x], x.order)
TaylorN{T<:Number}(x::T, order::Int) =
    TaylorN{T}([HomogeneousPolynomial([x], 0)], order)

# Shortcut to define TaylorN independent variables
"""
    TaylorN([T::Type=Float64], nv::Int; [order::Int=get_order()])

Shortcut to define the `nv`-th independent `TaylorN{T}` variable as a
polynomial. The order is defined through the keyword parameter `order`,
whose default corresponds to `get_order()`. The default of type for
`T` is `Float64`.

```julia
julia> TaylorN(1)
 1.0 x₁ + 𝒪(‖x‖⁷)
 
julia> TaylorN(Rational{Int},2)
 1//1 x₂ + 𝒪(‖x‖⁷)
```
"""
TaylorN{T<:Number}(::Type{T}, nv::Int; order::Int=get_order()) =
    return TaylorN( [HomogeneousPolynomial(T, nv)], order )
TaylorN(nv::Int; order::Int=get_order()) = TaylorN(Float64, nv, order=order)


# A `Number` which is not an `AbstractSeries`
const NumberNotSeries = Union{setdiff(subtypes(Number), [AbstractSeries])...}

# A `Number` which is not a TaylorN nor a HomogeneousPolynomial
const NumberNotSeriesN =
    Union{setdiff(subtypes(Number), [AbstractSeries])..., Taylor1}
