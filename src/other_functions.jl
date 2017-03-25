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
