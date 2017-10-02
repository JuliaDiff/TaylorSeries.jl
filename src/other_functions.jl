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

#taylor_expand! function for Taylor1
function taylor_expand!(a::Taylor1)
    #shifting around zero shouldn't change anything...
    nothing
end

function taylor_expand!{T<:Number}(a::Taylor1, x0::T)
    a.coeffs .= evaluate(a, Taylor1([x0,one(x0)], a.order) ).coeffs
    nothing
end

#taylor_expand! function for TaylorN
function taylor_expand!{T<:Number}(a::TaylorN,vals::Vector{T})
    a.coeffs .= evaluate(a,get_variables().+vals).coeffs
    nothing
end

function taylor_expand!(a::TaylorN)
    #shifting around zero shouldn't change anything...
    nothing
end
