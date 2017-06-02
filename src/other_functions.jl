# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)

    ## real, imag, conj and ctranspose ##
    for f in (:real, :imag, :conj)
        @eval ($f)(a::$T) = $T(($f).(a.coeffs), a.order)
    end

    # @eval real(x::Type{$T{T}}) where {T<:Number} = typeof(real(zero(x)))

    @eval ctranspose(a::$T) = conj.(a)

    ## isinf and isnan ##
    @eval isinf(a::$T) = any( isinf.(a.coeffs) )

    @eval isnan(a::$T) = any( isnan.(a.coeffs) )
end


## Division functions: rem and mod ##
for op in (:mod, :rem)
    for T in (:Taylor1, :TaylorN)
        @eval begin
            function ($op)(a::$T{T}, x::T) where T<:Real
                coeffs = copy(a.coeffs)
                @inbounds coeffs[1] = ($op)(constant_term(a), x)
                return $T(coeffs, a.order)
            end

            function ($op)(a::$T{T}, x::S) where {T<:Real, S<:Real}
                R = promote_type(T, S)
                a = convert($T{R}, a)
                return ($op)(a, convert(R,x))
            end
        end
    end

    @eval begin
        function ($op)(a::TaylorN{Taylor1{T}}, x::T) where T<:Real
            coeffs = copy(a.coeffs)
            @inbounds coeffs[1] = ($op)(constant_term(a), x)
            return TaylorN( coeffs, a.order )
        end

        function ($op)(a::TaylorN{Taylor1{T}}, x::S) where {T<:Real, S<:Real}
            R = promote_type(T,S)
            a = convert(TaylorN{Taylor1{R}}, a)
            return ($op)(a, convert(R,x))
        end

        function ($op)(a::Taylor1{TaylorN{T}}, x::T) where T<:Real
            coeffs = copy(a.coeffs)
            @inbounds coeffs[1] = ($op)(constant_term(a), x)
            return Taylor1( coeffs, a.order )
        end

        @inbounds function ($op)(a::Taylor1{TaylorN{T}}, x::S) where {T<:Real, S<:Real}
            R = promote_type(T,S)
            a = convert(Taylor1{TaylorN{R}}, a)
            return ($op)(a, convert(R,x))
        end
    end
end


## mod2pi and abs ##
for T in (:Taylor1, :TaylorN)
    @eval begin
        function mod2pi(a::$T{T}) where T<:Real
            coeffs = copy(a.coeffs)
            @inbounds coeffs[1] = mod2pi( constant_term(a) )
            return $T( coeffs, a.order)
        end

        function abs(a::$T{T}) where T<:Real
            if constant_term(a) > zero(T)
                return a
            elseif constant_term(a) < zero(T)
                return -a
            else
                throw(ArgumentError(
                """The 0th order Taylor1 coefficient must be non-zero
                (abs(x) is not differentiable at x=0)."""))
            end
        end
    end
end

function mod2pi(a::TaylorN{Taylor1{T}}) where T<:Real
    coeffs = copy(a.coeffs)
    @inbounds coeffs[1] = mod2pi( constant_term(a) )
    return TaylorN( coeffs, a.order )
end

function mod2pi(a::Taylor1{TaylorN{T}}) where T<:Real
    coeffs = copy(a.coeffs)
    @inbounds coeffs[1] = mod2pi( constant_term(a) )
    return Taylor1( coeffs, a.order )
end

function abs(a::TaylorN{Taylor1{T}}) where T<:Real
    if constant_term(a)[1] > zero(T)
        return a
    elseif constant_term(a)[1] < zero(T)
        return -a
    else
        throw(ArgumentError(
        """The 0th order TaylorN{Taylor1{T}} coefficient must be non-zero
        (abs(x) is not differentiable at x=0)."""))
    end
end

function abs(a::Taylor1{TaylorN{T}}) where {T<:Real}
    if constant_term(a[1]) > zero(T)
        return a
    elseif constant_term(a[1]) < zero(T)
        return -a
    else
        throw(ArgumentError(
        """The 0th order Taylor1{TaylorN{T}} coefficient must be non-zero
        (abs(x) is not differentiable at x=0)."""))
    end
end

doc"""
    abs(a)

Returns `a` if `constant_term(a) > 0` and `-a` if `constant_term(a) < 0` for
`a <:Union{Taylor1,TaylorN}`.
Notice that `typeof(abs(a)) <: AbstractSeries`.

""" abs


#norm
doc"""
    norm(x::AbstractSeries,p::Real)

Computes the p-norm of an `AbstractSeries` defined by
``P(\vec{x}) = \sum_k a_k \vec{x}^k`` as

```math
||P||_p =  \left( \sum_k ||a_k||_p^p \right)^{\frac{1}{p}}
```
which returns a non-negative number.

"""
norm(x::AbstractSeries, p::Real=2) = norm( norm.(x.coeffs, p), p)
norm(x::Union{Taylor1{T},HomogeneousPolynomial{T}}, p::Real=2) where
    {T<:NumberNotSeries} = norm(x.coeffs, p)

# rtoldefault
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval rtoldefault(::Type{$T{T}}) where {T<:Number} = rtoldefault(T)
end

# isfinite
doc"""
    isfinite(x::AbstractSeries) -> Bool

Test whether the coefficients of the polynomial `x` are finite.
"""
isfinite(x::AbstractSeries) = !isnan(x) && !isinf(x)

# isapprox; modified from Julia's Base.isapprox
doc"""
`isapprox(x::AbstractSeries, y::AbstractSeries;
    [rtol::Real=sqrt(eps), atol::Real=0, nans::Bool=false])`

Inexact equality comparison between polynomials: returns `true` if
`norm(x-y,1) <= atol + rtol*max(norm(x,1), norm(y,1))`, where `x` and `y` are
polynomials. For more details, see [`Base.isapprox`](@ref).
"""
function isapprox(x::T, y::S; rtol::Real=rtoldefault(x,y), atol::Real=0,
        nans::Bool=false) where {T<:AbstractSeries, S<:AbstractSeries}

    x == y || (isfinite(x) && isfinite(y) &&
        norm(x-y,1) <= atol + rtol*max(norm(x,1), norm(y,1))) ||
        (nans && isnan(x) && isnan(y))
end


#taylor_expand function for Taylor1
doc"""
    taylor_expand(f ,x0 ;order)

Makes a Taylor expansion of the function `f` around the point `x0`.

If x0 is a scalar, a `Taylor1` expansion will be done. If `x0` is a vector, a
`TaylorN` expansion will be computed. If the dimension of x0 (`length(x0)`)
is different from the variables set for `TaylorN` (`get_numvars()`), an
`AssertionError` will be thrown.
"""
function taylor_expand(f::Function; order::Int64=15)
   a = Taylor1(order)
   return f(a)
end

function taylor_expand(f::Function, x0::T; order::Int64=15) where {T<:Number}
   a = Taylor1([x0, one(T)], order)
   return f(a)
end

#taylor_expand function for TaylorN
function taylor_expand(f::Function, x0::Vector{T};
        order::Int64=get_order()) where {T<:Number}

    ll = length(x0)
    @assert ll == get_numvars() && order <= get_order()
    X = Array{TaylorN{T}}(ll)

    for i in eachindex(X)
        X[i] = x0[i] + TaylorN(T, i, order=order)
    end

    return f( X )
end

function taylor_expand(f::Function, x0...; order::Int64=get_order())
    x0 = promote(x0...)
    T = eltype(x0[1])
    ll = length(x0)
    @assert ll == get_numvars() && order <= get_order()
    X = Array{TaylorN{T}}(ll)

    for i in eachindex(X)
        X[i] = x0[i] + TaylorN(T, i, order=order)
    end

    return f( X... )
end

#update! function for Taylor1
doc"""
    update!(a,x0;order)

Takes `a <: Union{Taylo1,TaylorN}` and expands it around the coordinate `x0`.
"""
function update!(a::Taylor1, x0::T) where {T<:Number}
    a.coeffs .= evaluate(a, Taylor1([x0,one(x0)], a.order) ).coeffs
    nothing
end

#update! function for TaylorN
function update!(a::TaylorN,vals::Vector{T}) where {T<:Number}
    a.coeffs .= evaluate(a,get_variables().+vals).coeffs
    nothing
end

function update!(a::Union{Taylor1,TaylorN})
    #shifting around zero shouldn't change anything...
    nothing
end
