# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

## Equality ##
for T in (:Taylor1, :TaylorN)
    @eval begin
        ==(a::$T{T}, b::$T{S}) where {T<:Number,S<:Number} = ==(promote(a,b)...)

        function ==(a::$T{T}, b::$T{T}) where {T<:Number}
            if get_order(a) != get_order(b)
                a, b = fixorder(a, b)
            end
            return a.coeffs == b.coeffs
        end
    end
end

function ==(a::Taylor1{TaylorN{T}}, b::TaylorN{Taylor1{S}}) where {T, S}
    R = promote_type(T, S)
    return a == convert(Taylor1{TaylorN{R}}, b)
end
==(b::TaylorN{Taylor1{S}}, a::Taylor1{TaylorN{T}}) where {T, S} = a == b

function ==(a::HomogeneousPolynomial, b::HomogeneousPolynomial)
    get_order(a) == get_order(b) && return a.coeffs == b.coeffs
    return iszero(a.coeffs) && iszero(b.coeffs)
end


## Total ordering ##
for T in (:Taylor1, :TaylorN)
    @eval begin
        @inline function isless(a::$T{<:Number}, b::Real)
            a0 = constant_term(a)
            a0 != b && return isless(a0, b)
            nz = findfirst(a-b)
            if nz == -1
                return isless(a0, b)
            else
                return isless(a[nz], zero(b))
            end
        end
        @inline function isless(b::Real, a::$T{<:Number})
            a0 = constant_term(a)
            a0 != b && return isless(b, a0)
            nz = findfirst(b-a)
            if nz == -1
                return isless(b, a0)
            else
                return isless(zero(b), a[nz])
            end
        end
        #
        @inline isless(a::$T{<:NumberNotSeries}, b::$T{<:NumberNotSeries}) =
            isless(promote(a,b)...)
        @inline isless(a::$T{T}, b::$T{T}) where {T<:NumberNotSeries} =
            isless(a - b, zero(constant_term(a)))
    end
end


#=
# The following works for nested Taylor1s; iss #326.
ti = Taylor1(9)
to = Taylor1([zero(ti), one(ti)], 3)
too = Taylor1([zero(to), one(to)], 2)
tito = ti * to
ti > to > 0        # ok
ti^2 > tito > to^2 # ok
ti > ti^2 > to     # ok, in the sense that ti is "constant term" of series in to
=#
@inline isless(a::Taylor1{T}, b::Taylor1{S}) where {T<:Number, S<:Number} =
    isless(promote(a, b)...)
@inline isless(a::Taylor1{Taylor1{T}}, b::Taylor1{Taylor1{T}}) where
    {T<:Number} = isless(a - b, _zero_abstractfloat(a))

# Helper function to get the correct zero to use in `isless`
@inline _zero_abstractfloat(a::Taylor1{T}) where {T<:NumberNotSeries} = zero(T)
_zero_abstractfloat(a::Taylor1{Taylor1{T}}) where {T<:Number} =
    _zero_abstractfloat(constant_term(a))


@inline function isless(a::HomogeneousPolynomial{<:Number}, b::Real)
    get_order(a) == 0 && return isless(a[1], b)
    !iszero(b) && return isless(zero(a[1]), b)
    nz = max(findfirst(a), 1)
    return isless(a[nz], b)
end
@inline function isless(b::Real, a::HomogeneousPolynomial{<:Number})
    get_order(a) == 0 && return isless(b, a[1])
    !iszero(b) && return isless(b, zero(a[1]))
    nz = max(findfirst(a),1)
    return isless(b, a[nz])
end

@inline isless(a::HomogeneousPolynomial{T}, b::HomogeneousPolynomial{S}) where
    {T<:Number, S<:Number} = isless(promote(a,b)...)
@inline function isless(a::HomogeneousPolynomial{T},
        b::HomogeneousPolynomial{T}) where {T<:Number}
    orda = get_order(a)
    ordb = get_order(b)
    if orda == ordb
        return isless(a-b, zero(a[1]))
    elseif orda < ordb
        return isless(a, zero(a[1]))
    else
        return isless(-b, zero(a[1]))
    end
end

# Mixtures
@inline isless(a::Taylor1{TaylorN{<:NumberNotSeries}},
    b::Taylor1{TaylorN{<:NumberNotSeries}}) = isless(promote(a, b)...)
@inline isless(a::Taylor1{TaylorN{T}}, b::Taylor1{TaylorN{T}}) where
    {T<:NumberNotSeries} = isless(a - b, zero(T))

@inline isless(a::HomogeneousPolynomial{Taylor1{<:NumberNotSeries}},
        b::HomogeneousPolynomial{Taylor1{<:NumberNotSeries}}) =
    isless(promote(a, b)...)
@inline function isless(a::HomogeneousPolynomial{Taylor1{T}},
        b::HomogeneousPolynomial{Taylor1{T}}) where {T<:NumberNotSeries}
    orda = get_order(a)
    ordb = get_order(b)
    if orda == ordb
        return isless(a-b, zero(T))
    elseif orda < ordb
        return isless(a, zero(T))
    else
        return isless(-b, zero(T))
    end
end

@inline isless(a::TaylorN{Taylor1{<:NumberNotSeries}},
    b::TaylorN{Taylor1{<:NumberNotSeries}}) = isless(promote(a, b)...)
@inline isless(a::TaylorN{Taylor1{T}}, b::TaylorN{Taylor1{T}}) where
    {T<:NumberNotSeries} = isless(a - b, zero(T))


@doc doc"""
    isless(a::Taylor1{<:Real}, b::Real)
    isless(a::TaylorN{<:Real}, b::Real)

Compute `isless` by comparing the `constant_term(a)` and `b`. If they are equal,
returns `a[nz] < 0`, with `nz` the first
non-zero coefficient after the constant term. This defines a total order.

For many variables, the ordering includes a lexicographical convention in order to be
total. We have opted for the simplest one: the *larger* variable appears *before*
for the `TaylorN` variables are defined (e.g., through [`set_variables`](@ref)).

For nested `Taylor1{Taylor1{...}}`s, the ordering is established by which one
is a `constant_term` of the other.

Refs:
- M. Berz, AIP Conference Proceedings 177, 275 (1988); https://doi.org/10.1063/1.37800
- M. Berz, "Automatic Differentiation as Nonarchimedean Analysis", Computer Arithmetic and
    Enclosure Methods, (1992), Elsevier, 439-450.

---

    isless(a::Taylor1{T}, b::Taylor1{T})
    isless(a::TaylorN{T}, b::Taylor1{T})

Returns `isless(a - b, zero(T))`.
""" isless

