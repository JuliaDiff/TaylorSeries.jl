# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

## real, imag, conj and ctranspose ##
for f in (:real, :imag, :conj)
    @eval ($f){T<:Number}(a::Taylor1{T}) = Taylor1(($f)(a.coeffs), a.order)
end
ctranspose{T<:Number}(a::Taylor1{T}) = conj(a)


## real, imag, conj and ctranspose ##
for TT in (:HomogeneousPolynomial, :TaylorN), f in (:real, :imag, :conj)
    @eval ($f){T<:Number}(a::($TT){T}) = ($TT)(($f)(a.coeffs), a.order)
end

ctranspose{T<:Number}(a::HomogeneousPolynomial{T}) = conj(a)
ctranspose{T<:Number}(a::TaylorN{T}) = conj(a)


# Tests `isinf` and `isnan` for *any* of the polynomial coefficients
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval begin
        isinf(a::$T) = any( isinf.(a.coeffs) )
        isnan(a::$T) = any( isnan.(a.coeffs) )
    end
end

## Division functions: rem and mod
for op in (:mod, :rem)
    @eval begin
        function ($op){T<:Real}(a::Taylor1{T}, x::T)
            coeffs = copy(a.coeffs)
            @inbounds coeffs[1] = ($op)(a[1], x)
            return Taylor1(coeffs, a.order)
        end
        function ($op){T<:Real, S<:Real}(a::Taylor1{T}, x::S)
            R = promote_type(T,S)
            a = convert(Taylor1{R}, a)
            return ($op)(a, convert(R,x))
        end
    end
end

function mod2pi{T<:Real}(a::Taylor1{T})
    coeffs = copy(a.coeffs)
    @inbounds coeffs[1] = mod2pi( a[1] )
    return Taylor1( coeffs, a.order)
end

## abs function ##
"""
    abs(a)

Return `a` or `-a` depending on the 0-th order coefficient of the
`Taylor1` polynomial `a`.
If `a[1]` is zero, an `ArgumentError` is thrown.
"""
function abs{T<:Real}(a::Taylor1{T})
    if a[1] > zero(T)
        return a
    elseif a[1] < zero(T)
        return -a
    else
        throw(ArgumentError(
        """The 0th order Taylor1 coefficient must be non-zero
        (abs(x) is not differentiable at x=0)."""))
    end
end



## Division functions: rem and mod
for op in (:mod, :rem)
    @eval begin
        @inbounds function ($op){T<:Real,S<:Real}(a::TaylorN{T}, x::S)
            coeffs = copy(a.coeffs)
            y = ($op)(a[1][1], x)
            coeffs[1] = HomogeneousPolynomial([y], 0)
            return TaylorN( coeffs, a.order )
        end
        @inbounds function ($op){T<:Real}(a::TaylorN{Taylor1{T}}, x::T)
            coeffs = copy(a.coeffs)
            y = ($op)(a[1][1], x)
            coeffs[1] = HomogeneousPolynomial([y], 0)
            return TaylorN( coeffs, a.order )
        end
        @inbounds function ($op){T<:Real,S<:Real}(a::TaylorN{Taylor1{T}}, x::S)
            R = promote_type(T,S)
            a = convert(TaylorN{Taylor1{R}}, a)
            return ($op)(a, convert(R,x))
        end
        @inbounds function ($op){T<:Real}(a::Taylor1{TaylorN{T}}, x::T)
            coeffs = copy(a.coeffs)
            y = ($op)(a[1], x)
            coeffs[1] = y
            return Taylor1( coeffs, a.order )
        end
        @inbounds function ($op){T<:Real,S<:Real}(a::Taylor1{TaylorN{T}}, x::S)
            R = promote_type(T,S)
            a = convert(Taylor1{TaylorN{R}}, a)
            return ($op)(a, convert(R,x))
        end
    end
end

function mod2pi{T<:Real}(a::TaylorN{T})
    coeffs = copy(a.coeffs)
    @inbounds y = mod2pi(a[1][1])
    @inbounds coeffs[1] = HomogeneousPolynomial([y], 0)
    return TaylorN( coeffs, a.order )
end
function mod2pi{T<:Real}(a::TaylorN{Taylor1{T}})
    coeffs = copy(a.coeffs)
    @inbounds y = mod2pi(a[1][1])
    @inbounds coeffs[1] = HomogeneousPolynomial([y], 0)
    return TaylorN( coeffs, a.order )
end
function mod2pi{T<:Real}(a::Taylor1{TaylorN{T}})
    coeffs = copy(a.coeffs)
    @inbounds y = mod2pi(a[1])
    @inbounds coeffs[1] = y
    return Taylor1( coeffs, a.order )
end

## abs function ##
"""
    abs(a::TaylorN)

Absolute value of a `TaylorN` polynomial, using the 0-th order coefficient.

Return `a` or `-a`, depending on the 0-th order coefficient of `a`.
If it is zero, it throws an `ArgumentError`.
"""
function abs{T<:Real}(a::TaylorN{T})
    if a[1][1] > zero(T)
        return a
    elseif a[1][1] < zero(T)
        return -a
    else
        throw(ArgumentError(
        """The 0th order TaylorN coefficient must be non-zero
        (`abs(x)` is not differentiable at zero)."""))
    end
end
