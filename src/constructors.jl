# This file is part of the Taylor1Series.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

"""
    AbstractSeries{T<:Number} <: Number

Parameterized abstract type for [`Taylor1`](@ref),
[`HomogeneousPolynomial`](@ref) and [`TaylorN`](@ref).
"""
abstract type AbstractSeries{T<:Number} <: Number end


## Constructors ##

######################### Taylor1
"""
    Taylor1{T<:Number} <: AbstractSeries{T}

DataType for polynomial expansions in one independent variable.

**Fields:**

- `coeffs :: Array{T,1}` Expansion coefficients; the ``i``-th
    component is the coefficient of degree ``i-1`` of the expansion.
- `order  :: Int` Maximum order (degree) of the polynomial.

Note that `Taylor1` variables are callable. For more information, see
[`evaluate`](@ref).
"""
struct Taylor1{T<:Number} <: AbstractSeries{T}
    coeffs :: Array{T,1}
    order :: Int

    ## Inner constructor ##
    function Taylor1{T}(coeffs::Array{T,1}, order::Int) where T<:Number
        resize_coeffs1!(coeffs, order)
        return new{T}(coeffs, order)
    end
end

## Outer constructors ##
Taylor1(x::Taylor1{T}) where {T<:Number} = x
Taylor1(coeffs::Array{T,1}, order::Int) where {T<:Number} = Taylor1{T}(coeffs, order)
Taylor1(coeffs::Array{T,1}) where {T<:Number} = Taylor1(coeffs, length(coeffs)-1)
function Taylor1(x::T, order::Int) where {T<:Number}
    v = fill(zero(x), order+1)
    v[1] = x
    return Taylor1(v, order)
end

# Methods using 1-d views to create Taylor1's
Taylor1(a::SubArray{T,1}, order::Int) where {T<:Number} = Taylor1(a.parent[a.indices...], order)
Taylor1(a::SubArray{T,1}) where {T<:Number} = Taylor1(a.parent[a.indices...])


# Shortcut to define Taylor1 independent variables
"""
    Taylor1([T::Type=Float64], order::Int)

Shortcut to define the independent variable of a `Taylor1{T}` polynomial of
given `order`. The default type for `T` is `Float64`.

```julia
julia> Taylor1(16)
 1.0 t + 𝒪(t¹⁷)

julia> Taylor1(Rational{Int}, 4)
 1//1 t + 𝒪(t⁵)
```
"""
Taylor1(::Type{T}, order::Int) where {T<:Number} = Taylor1( [zero(T), one(T)], order)
Taylor1(order::Int) = Taylor1(Float64, order)

######################### HomogeneousPolynomial
"""
    HomogeneousPolynomial{T<:Number} <: AbstractSeries{T}

DataType for homogenous polynomials in many (>1) independent variables.

**Fields:**

- `coeffs  :: Array{T,1}` Expansion coefficients of the homogeneous
polynomial; the ``i``-th component is related to a monomial, where the degrees
of the independent variables are specified by `coeff_table[order+1][i]`.
- `order   :: Int` order (degree) of the homogenous polynomial.

Note that `HomogeneousPolynomial` variables are callable. For more information,
see [`evaluate`](@ref).
"""
struct HomogeneousPolynomial{T<:Number} <: AbstractSeries{T}
    coeffs  :: Array{T,1}
    order   :: Int

    function HomogeneousPolynomial{T}(coeffs::Array{T,1}, order::Int) where T<:Number
        resize_coeffsHP!(coeffs, order)
        return new{T}(coeffs, order)
    end
end

HomogeneousPolynomial(x::HomogeneousPolynomial{T}) where {T<:Number} = x
HomogeneousPolynomial(coeffs::Array{T,1}, order::Int) where {T<:Number} =
    HomogeneousPolynomial{T}(coeffs, order)
HomogeneousPolynomial(coeffs::Array{T,1}) where {T<:Number} =
    HomogeneousPolynomial(coeffs, orderH(coeffs))
HomogeneousPolynomial(x::T, order::Int) where {T<:Number} =
    HomogeneousPolynomial([x], order)

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
function HomogeneousPolynomial(::Type{T}, nv::Int) where {T<:Number}
    @assert 0 < nv ≤ get_numvars()
    v = zeros(T, get_numvars())
    @inbounds v[nv] = one(T)
    return HomogeneousPolynomial(v, 1)
end
HomogeneousPolynomial(nv::Int) = HomogeneousPolynomial(Float64, nv)



######################### TaylorN
"""
    TaylorN{T<:Number} <: AbstractSeries{T}

DataType for polynomial expansions in many (>1) independent variables.

**Fields:**

- `coeffs  :: Array{HomogeneousPolynomial{T},1}` Vector containing the
`HomogeneousPolynomial` entries. The ``i``-th component corresponds to the
homogeneous polynomial of degree ``i-1``.
- `order   :: Int`  maximum order of the polynomial expansion.

Note that `TaylorN` variables are callable. For more information, see
[`evaluate`](@ref).
"""
struct TaylorN{T<:Number} <: AbstractSeries{T}
    coeffs  :: Array{HomogeneousPolynomial{T},1}
    order   :: Int

    function TaylorN{T}(v::Array{HomogeneousPolynomial{T},1}, order::Int) where T<:Number
        coeffs = zeros(HomogeneousPolynomial{T}, order)
        @inbounds for i in eachindex(v)
            ord = v[i].order
            if ord ≤ order
                coeffs[ord+1] += v[i]
            end
        end
        new{T}(coeffs, order)
    end
end

TaylorN(x::TaylorN{T}) where T<:Number = x
function TaylorN(x::Array{HomogeneousPolynomial{T},1}, order::Int) where {T<:Number}
    if order == 0
        order = maxorderH(x)
    end
    return TaylorN{T}(x, order)
end
TaylorN(x::Array{HomogeneousPolynomial{T},1}) where {T<:Number} =
    TaylorN(x, maxorderH(x))
TaylorN(x::HomogeneousPolynomial{T}, order::Int) where {T<:Number} =
    TaylorN( [x], order )
TaylorN(x::HomogeneousPolynomial{T}) where {T<:Number} = TaylorN(x, x.order)
TaylorN(x::T, order::Int) where {T<:Number} =
    TaylorN(HomogeneousPolynomial([x], 0), order)

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
TaylorN(::Type{T}, nv::Int; order::Int=get_order()) where {T<:Number} =
    TaylorN( HomogeneousPolynomial(T, nv), order )
TaylorN(nv::Int; order::Int=get_order()) = TaylorN(Float64, nv, order=order)



# A `Number` which is not an `AbstractSeries`
const NumberNotSeries = Union{setdiff(subtypes(Number), [AbstractSeries])...}

# A `Number` which is not `TaylorN` nor a `HomogeneousPolynomial`
const NumberNotSeriesN =
    Union{setdiff(subtypes(Number), [AbstractSeries])..., Taylor1}

## Additional Taylor1 outer constructor ##
Taylor1{T}(x::S) where {T<:Number, S<:NumberNotSeries} = Taylor1([convert(T,x)], 0)
