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

    @eval ctranspose(a::$T) = conj.(a)

    ## isinf and isnan ##
    @eval isinf(a::$T) = any( isinf.(a.coeffs) )

    @eval isnan(a::$T) = any( isnan.(a.coeffs) )
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

        function ($op){T<:Real,S<:Real}(a::TaylorN{T}, x::S)
            coeffs = copy(a.coeffs)
            # y = ($op)(a[1][1], x)
            @inbounds coeffs[1] = HomogeneousPolynomial([($op)(a[1][1], x)], 0)
            return TaylorN( coeffs, a.order )
        end

        function ($op){T<:Real}(a::TaylorN{Taylor1{T}}, x::T)
            coeffs = copy(a.coeffs)
            # y = ($op)(a[1][1], x)
            @inbounds coeffs[1] = HomogeneousPolynomial([($op)(a[1][1], x)], 0)
            return TaylorN( coeffs, a.order )
        end

        function ($op){T<:Real,S<:Real}(a::TaylorN{Taylor1{T}}, x::S)
            R = promote_type(T,S)
            a = convert(TaylorN{Taylor1{R}}, a)
            return ($op)(a, convert(R,x))
        end

        function ($op){T<:Real}(a::Taylor1{TaylorN{T}}, x::T)
            coeffs = copy(a.coeffs)
            @inbounds y = ($op)(a[1], x)
            @inbounds coeffs[1] = y
            return Taylor1( coeffs, a.order )
        end

        @inbounds function ($op){T<:Real,S<:Real}(a::Taylor1{TaylorN{T}}, x::S)
            R = promote_type(T,S)
            a = convert(Taylor1{TaylorN{R}}, a)
            return ($op)(a, convert(R,x))
        end
    end
end


function mod2pi{T<:Real}(a::Taylor1{T})
    coeffs = copy(a.coeffs)
    @inbounds coeffs[1] = mod2pi( a[1] )
    return Taylor1( coeffs, a.order)
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
