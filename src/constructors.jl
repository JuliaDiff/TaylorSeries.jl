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

- `coeffs :: FixedSizeVectorDefault{T}` Expansion coefficients; the ``i``-th
    component is the coefficient of degree ``i-1`` of the expansion.

Note that `Taylor1` variables are callable. For more information, see
[`evaluate`](@ref).
"""
struct Taylor1{T<:Number} <: AbstractSeries{T}
        coeffs :: FixedSizeVectorDefault{T}
    ## Inner constructors ##
    function Taylor1{T}(coeffs::FixedSizeVectorDefault{T}) where {T<:Number}
        return new{T}(coeffs)
    end
    function Taylor1{T}(coeffs::AbstractVector{T}, order::Int) where {T<:Number}
        v = FixedSizeVectorDefault{T}(undef, order+1)
        minrange = min(eachindex(v), eachindex(coeffs))
        for ord in minrange
            v[ord] = coeffs[ord]
        end
        for ord in last(minrange)+1:order+1
            v[ord] = zero(coeffs[1])
        end
        return new{T}(v)
    end
end

## Outer constructors ##
Taylor1(x::Taylor1{T}) where {T<:Number} = x
Taylor1(coeffs::AbstractArray{T,1}, order::Int) where {T<:Number} =
    Taylor1{T}(coeffs, order)
Taylor1(coeffs::AbstractArray{T,1}) where {T<:Number} =
    Taylor1{T}(FixedSizeVectorDefault(coeffs))
Taylor1(coeffs::FixedSizeVectorDefault{T}) where {T<:Number} = Taylor1{T}(coeffs)
Taylor1(coeffs::FixedSizeVectorDefault{T}, order::Int) where {T<:Number} =
    Taylor1{T}(coeffs, order)
function Taylor1(x::T, order::Int) where {T<:Number}
    v = FixedSizeVectorDefault{T}(undef, order+1)
    v .= zero.(x)
    v[1] = deepcopy(x)
    return Taylor1{T}(v)
end


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
# Taylor1(::Type{T}, order::Int) where {T<:Number} = Taylor1( [zero(T), one(T)], order)
function Taylor1(::Type{T}, order::Int) where {T<:Number}
    coeffs = FixedSizeVectorDefault{T}(undef, order+1)
    coeffs .= zero(T)
    coeffs[2] = one(T)
    return Taylor1{T}(coeffs)
end
# Taylor1(order::Int) = Taylor1([0.0, 1.0], order)
function Taylor1(order::Int)
    coeffs = FixedSizeVectorDefault{Float64}(undef, order+1)
    coeffs .= 0.0
    coeffs[2] = 1.0
    return Taylor1{Float64}(coeffs)
end


######################### HomogeneousPolynomial
"""
    HomogeneousPolynomial{T<:Number} <: AbstractSeries{T}

DataType for homogeneous polynomials in many (>1) independent variables.

**Fields:**

- `coeffs  :: Array{T,1}` Expansion coefficients of the homogeneous
polynomial; the ``i``-th component is related to a monomial, where the degrees
of the independent variables are specified by `coeff_table[order+1][i]`.
- `order   :: Int` order (degree) of the homogeneous polynomial.

Note that `HomogeneousPolynomial` variables are callable. For more information,
see [`evaluate`](@ref).
"""
struct HomogeneousPolynomial{T<:Number} <: AbstractSeries{T}
        coeffs  :: FixedSizeVectorDefault{T}
        order   :: Int
    ## Inner constructors ##
    HomogeneousPolynomial{T}(x::T, order::Int) where {T<:Number} =
        new{T}(_coeffsHP(x, order), order)
    HomogeneousPolynomial{T}(coeffs::AbstractArray{T,1}, order::Int) where {T<:Number} =
        new{T}(_coeffsHP(coeffs, order), order)
end

HomogeneousPolynomial(x::HomogeneousPolynomial{T}) where {T<:Number} = x
HomogeneousPolynomial(coeffs::AbstractArray{T,1}, order::Int) where {T<:Number} =
    HomogeneousPolynomial{T}(coeffs, order)
HomogeneousPolynomial(coeffs::AbstractArray{T,1}) where {T<:Number} =
    HomogeneousPolynomial{T}(coeffs, orderH(coeffs))
HomogeneousPolynomial(x::T, order::Int) where {T<:Number} =
    HomogeneousPolynomial{T}(x, order)

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
    v = FixedSizeVectorDefault{T}(undef, get_numvars())
    v .= zero(T)
    v[nv] = one(T)
    return HomogeneousPolynomial(v, 1)
end
function HomogeneousPolynomial(nv::Int)
    res = HomogeneousPolynomial(0.0, 1)
    res[nv] = 1.0
    return res
end



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
        coeffs  :: FixedSizeVectorDefault{HomogeneousPolynomial{T}}
    ## Inner constructor ##
    function TaylorN{T}(v::AbstractArray{HomogeneousPolynomial{T},1},
            order::Int) where {T<:Number}
        isempty(v) &&
            return new{T}(zeros(HomogeneousPolynomial(zero(T),0), order))
        coeffs = _coeffsTN(v, order)
        return new{T}(coeffs)
    end
end

TaylorN(x::TaylorN{T}) where {T<:Number} = x
function TaylorN(x::AbstractArray{HomogeneousPolynomial{T},1},
        order::Int) where {T<:Number}
    if order == 0
        order = maxorderH(x)
    end
    return TaylorN{T}(x, order)
end
TaylorN(x::AbstractArray{HomogeneousPolynomial{T},1}) where {T<:Number} =
    TaylorN{T}(x, maxorderH(x))
TaylorN(x::HomogeneousPolynomial{T}, order::Int) where {T<:Number} =
    TaylorN{T}([x], order )
TaylorN(x::HomogeneousPolynomial{T}) where {T<:Number} =
    TaylorN{T}([x], get_order(x))
TaylorN(x::T, order::Int) where {T<:Number} =
    TaylorN(HomogeneousPolynomial(x, 0), order)

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
const NumberNotSeries = Union{Real,Complex}

# A `Number` which is not `TaylorN` nor a `HomogeneousPolynomial`
const NumberNotSeriesN = Union{Real,Complex,Taylor1}

## Additional Taylor1 and TaylorN outer constructor ##
Taylor1{T}(x::S) where {T<:Number,S<:NumberNotSeries} = Taylor1([convert(T,x)], 0)
TaylorN{T}(x::S) where {T<:Number,S<:NumberNotSeries} = TaylorN(convert(T, x), TaylorSeries.get_order())


# """
#     get_numvars
#
# Return the number of variables of a `Taylor1`, `HomogeneousPolynomial`
# or `TaylorN` object.
# """
get_numvars(t::Number) = 0
get_numvars(t::Taylor1) = 1
get_numvars(t::Taylor1{Taylor1{T}}) where {T<:Number} = get_numvars(t[0])+1
get_numvars(::T) where {T<:Union{HomogeneousPolynomial, TaylorN}} = _params_TaylorN_.num_vars
