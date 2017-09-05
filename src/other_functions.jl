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

    @eval real{T<:Number}(x::Type{$T{T}}) = typeof(real(zero(x)))

    @eval ctranspose(a::$T) = conj.(a)

    ## isinf and isnan ##
    @eval isinf(a::$T) = any( isinf.(a.coeffs) )

    @eval isnan(a::$T) = any( isnan.(a.coeffs) )
end


## Division functions: rem and mod ##
for op in (:mod, :rem)
    for T in (:Taylor1, :TaylorN)
        @eval begin
            function ($op){T<:Real}(a::$T{T}, x::T)
                coeffs = copy(a.coeffs)
                @inbounds coeffs[1] = ($op)(constant_term(a), x)
                return $T(coeffs, a.order)
            end

            function ($op){T<:Real, S<:Real}(a::$T{T}, x::S)
                R = promote_type(T, S)
                a = convert($T{R}, a)
                return ($op)(a, convert(R,x))
            end
        end
    end

    @eval begin
        function ($op){T<:Real}(a::TaylorN{Taylor1{T}}, x::T)
            coeffs = copy(a.coeffs)
            @inbounds coeffs[1] = ($op)(constant_term(a), x)
            return TaylorN( coeffs, a.order )
        end

        function ($op){T<:Real,S<:Real}(a::TaylorN{Taylor1{T}}, x::S)
            R = promote_type(T,S)
            a = convert(TaylorN{Taylor1{R}}, a)
            return ($op)(a, convert(R,x))
        end

        function ($op){T<:Real}(a::Taylor1{TaylorN{T}}, x::T)
            coeffs = copy(a.coeffs)
            @inbounds coeffs[1] = ($op)(constant_term(a), x)
            return Taylor1( coeffs, a.order )
        end

        @inbounds function ($op){T<:Real,S<:Real}(a::Taylor1{TaylorN{T}}, x::S)
            R = promote_type(T,S)
            a = convert(Taylor1{TaylorN{R}}, a)
            return ($op)(a, convert(R,x))
        end
    end
end


## mod2pi and abs ##
for T in (:Taylor1, :TaylorN)
    @eval begin
        function mod2pi{T<:Real}(a::$T{T})
            coeffs = copy(a.coeffs)
            @inbounds coeffs[1] = mod2pi( constant_term(a) )
            return $T( coeffs, a.order)
        end

        function abs{T<:Real}(a::$T{T})
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


function mod2pi{T<:Real}(a::TaylorN{Taylor1{T}})
    coeffs = copy(a.coeffs)
    @inbounds coeffs[1] = mod2pi( constant_term(a) )
    return TaylorN( coeffs, a.order )
end

function mod2pi{T<:Real}(a::Taylor1{TaylorN{T}})
    coeffs = copy(a.coeffs)
    @inbounds coeffs[1] = mod2pi( constant_term(a) )
    return Taylor1( coeffs, a.order )
end

function abs{T<:Real}(a::TaylorN{Taylor1{T}})
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

function abs{T<:Real}(a::Taylor1{TaylorN{T}})
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

Computes the p-norm of an `AbstractSeries` defined by ``P(\vec{x}) = \sum_k a_k \vec{x}^k`` as

```math
||P||_p =  \left( \sum_k ||a_k||_p^p \right)^{\frac{1}{p}}
```
which returns a non-negative number.

"""
norm(x::AbstractSeries, p::Real=2) = norm( norm.(x.coeffs, p), p)
norm{T<:NumberNotSeries}(x::Union{Taylor1{T}, HomogeneousPolynomial{T}}, p::Real=2) = norm(x.coeffs, p)

# rtoldefault
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval rtoldefault{T<:Number}(::Type{$T{T}}) = rtoldefault(T)
end

# isfinite
doc"""
    isfinite(x::AbstractSeries) -> Bool

Test whether the coefficients of the polynomial `x` are finite.
"""
isfinite(x::AbstractSeries) = !isnan(x) && !isinf(x)

# isapprox; modified from Julia's Base.isapprox
doc"""
    isapprox(x::AbstractSeries, y::AbstractSeries; [rtol::Real=sqrt(eps), atol::Real=0, nans::Bool=false])

Inexact equality comparison between polynomials: returns `true` if
`norm(x-y,1) <= atol + rtol*max(norm(x,1), norm(y,1))`, where `x` and `y` are
polynomials. For more details, see [`Base.isapprox`](@ref).
"""
function isapprox{T<:AbstractSeries,S<:AbstractSeries}(x::T, y::S; rtol::Real=rtoldefault(x,y), atol::Real=0, nans::Bool=false)
    x == y || (isfinite(x) && isfinite(y) && norm(x-y,1) <= atol + rtol*max(norm(x,1), norm(y,1))) || (nans && isnan(x) && isnan(y))


#taylor_expand function for Taylor1
function taylor_expand(f::Function; order::Int64=15)
   a = Taylor1(order)
   return f(a)
end

function taylor_expand{T<:Number}(f::Function, x0::T; order::Int64=15)
   a = Taylor1([x0, one(T)], order)
   return f(a)
end

#taylor_expand function for TaylorN
function taylor_expand{T<:Number}(f::Function, x0::Vector{T}; order::Int64=get_order()) #a Taylor expansion around x0
    ll = length(x0)
    ll == get_numvars() ? X = get_variables() : begin
        X = set_variables("x",order=order,numvars=ll)
        warn("Changed number of TaylorN variables to $ll.")
        end
    return f(X.+x0)
end

function taylor_expand(f::Function, x0...; order::Int64=get_order()) #a Taylor expansion around x0
    ll = length(x0)
    ll == get_numvars() ? X = get_variables() : begin
        X = set_variables("x",order=order,numvars=ll)
        warn("Changed number of TaylorN variables to $ll.")
        end

    return f(x0 .+ X...)
end
