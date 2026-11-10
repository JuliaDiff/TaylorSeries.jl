# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#


## real, imag, conj and ctranspose ##
for f in (:real, :imag, :conj)
    @eval ($f)(a::Taylor1) = Taylor1(($f)(a.coeffs), order(a))
    for T in (:HomogeneousPolynomial, :TaylorN)
        @eval ($f)(a::$T) = $T(a.space, ($f)(a.coeffs), order(a))
    end
end

for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval adjoint(a::$T) = conj(a)

    ## isinf and isnan ##
    @eval isinf(a::$T) = any(isinf, a.coeffs)

    @eval isnan(a::$T) = any(isnan, a.coeffs)
end


# Rounding
round(x::AbstractSeries, r::RoundingMode; kwargs...) = round(constant_term(x), r; kwargs...)
round(::Type{T}, x::AbstractSeries; kwargs...) where {T <: NumberNotSeries} =
    round(T, x, RoundNearest; kwargs...)
round(::Type{T}, x::AbstractSeries, r::RoundingMode; kwargs...) where {T <: NumberNotSeries} =
    round(T, constant_term(x), r; kwargs...)

ceil(x::AbstractSeries; kwargs...) = ceil(constant_term(x); kwargs...)
ceil(::Type{T}, x::AbstractSeries; kwargs...) where {T <: NumberNotSeries} =
    ceil(T, constant_term(x); kwargs...)

floor(x::AbstractSeries; kwargs...) = floor(constant_term(x); kwargs...)
floor(::Type{T}, x::AbstractSeries; kwargs...) where {T <: NumberNotSeries} =
    floor(T, constant_term(x); kwargs...)


## Division functions: rem and mod ##
for op in (:mod, :rem)
    for T in (:Taylor1, :TaylorN)
        @eval begin
            function ($op)(a::$T{T}, x::T) where {T<:Real}
                c = $T(copy(a.coeffs), order(a))
                constant_term!(c, ($op)(constant_term(a), x))
                return c
            end

            function ($op)(a::$T{T}, x::S) where {T<:Real,S<:Real}
                R = promote_type(T, S)
                b = convert($T{R}, a)
                return ($op)(b, convert(R,x))
            end
        end
    end

    @eval begin
        function ($op)(a::TaylorN{Taylor1{T}}, x::T) where {T<:Real}
            c = TaylorN(copy(a.coeffs), order(a))
            constant_term!(c, ($op)(constant_term(a), x))
            return c
        end

        function ($op)(a::TaylorN{Taylor1{T}}, x::S) where {T<:Real,S<:Real}
            R = promote_type(T,S)
            a = convert(TaylorN{Taylor1{R}}, a)
            return ($op)(a, convert(R,x))
        end

        function ($op)(a::Taylor1{TaylorN{T}}, x::T) where {T<:Real}
            c = Taylor1(copy(a.coeffs), order(a))
            constant_term!(c, ($op)(constant_term(a), x))
            return c
        end

        @inbounds function ($op)(a::Taylor1{TaylorN{T}}, x::S) where {T<:Real,S<:Real}
            R = promote_type(T,S)
            a = convert(Taylor1{TaylorN{R}}, a)
            return ($op)(a, convert(R,x))
        end
    end
end


## mod2pi and abs ##
for T in (:Taylor1, :TaylorN)
    @eval begin
        function mod2pi(a::$T{T}) where {T<:Real}
            c = $T(copy(a.coeffs), order(a))
            constant_term!(c, mod2pi(constant_term(a)))
            return c
        end

        function abs(a::$T{T}) where {T<:Real}
            if constant_term(a) > 0
                return a
            elseif constant_term(a) < 0
                return -a
            else
                throw(DomainError(a,
                """The 0th order Taylor1 coefficient must be non-zero
                (abs(x) is not differentiable at x=0)."""))
            end
        end

        abs2(a::$T) = real(a)^2 + imag(a)^2

        abs(x::$T{T}) where {T<:Complex} = sqrt(abs2(x))
    end
end

function mod2pi(a::TaylorN{Taylor1{T}}) where {T<:Real}
    c = TaylorN(copy(a.coeffs), order(a))
    constant_term!(c, mod2pi(constant_term(a)))
    return c
end

function mod2pi(a::Taylor1{TaylorN{T}}) where {T<:Real}
    c = Taylor1(copy(a.coeffs), order(a))
    constant_term!(c, mod2pi(constant_term(a)))
    return c
end

function abs(a::TaylorN{Taylor1{T}}) where {T<:Real}
    if constant_term(a)[0] > 0
        return a
    elseif constant_term(a)[0] < 0
        return -a
    else
        throw(DomainError(a,
        """The 0th order TaylorN{Taylor1{T}} coefficient must be non-zero
        (abs(x) is not differentiable at x=0)."""))
    end
end

function abs(a::Taylor1{TaylorN{T}}) where {T<:Real}
    if constant_term(a[0]) > 0
        return a
    elseif constant_term(a[0]) < 0
        return -a
    else
        throw(DomainError(a,
        """The 0th order Taylor1{TaylorN{T}} coefficient must be non-zero
        (abs(x) is not differentiable at x=0)."""))
    end
end

abs(x::Taylor1{TaylorN{T}}) where {T<:NumberNotSeries} = sqrt(abs2(x))
abs(x::TaylorN{Taylor1{T}}) where {T<:NumberNotSeries} = sqrt(abs2(x))
abs(x::Taylor1{Taylor1{T}}) where {T<:Number} = sqrt(abs2(x))

@doc doc"""
    abs(a::Taylor1{T})
    abs(a::TaylorN{T})

For `T<:Real` type returns `a` if `constant_term(a) > 0` and `-a` if
`constant_term(a) < 0`. For other `Number` types, it returns `sqrt(abs2(a))`.

""" abs


#norm
@doc doc"""
    norm(x::AbstractSeries, p::Real)

Returns the p-norm of an `x::AbstractSeries`, defined by

```math
\begin{equation*}
\left\Vert x \right\Vert_p =  \left( \sum_k | x_k |^p \right)^{\frac{1}{p}},
\end{equation*}
```
which returns a non-negative number.

""" norm

norm(x::AbstractSeries, p::Real=2) = norm( norm.(x.coeffs, p), p)
#norm for Taylor vectors
norm(v::Vector{T}, p::Real=2) where {T<:AbstractSeries} = norm( norm.(v, p), p)

# rtoldefault
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval rtoldefault(::Type{$T{T}}) where {T<:Number} = rtoldefault(T)
    @eval rtoldefault(::$T{T}) where {T<:Number} = rtoldefault(T)
end


# isfinite
"""
    isfinite(x::AbstractSeries) -> Bool

Test whether the coefficients of the polynomial `x` are finite.
"""
isfinite(x::AbstractSeries) = !isnan(x) && !isinf(x)


# isapprox; modified from Julia's Base.isapprox
"""
    isapprox(x::AbstractSeries, y::AbstractSeries; rtol::Real=sqrt(eps), atol::Real=0, nans::Bool=false)

Inexact equality comparison between polynomials: returns `true` if
`norm(x-y,1) <= atol + rtol*max(norm(x,1), norm(y,1))`, where `x` and `y` are
polynomials. For more details, see [`Base.isapprox`](@ref).
"""
function isapprox(x::T, y::S; rtol::Real=rtoldefault(x,y,0), atol::Real=0.0,
        nans::Bool=false) where {T<:AbstractSeries,S<:AbstractSeries}

    x == y || (isfinite(x) && isfinite(y) &&
        norm(x-y,1) <= atol + rtol*max(norm(x,1), norm(y,1))) ||
        (nans && isnan(x) && isnan(y))
end
#isapprox for vectors of Taylors
function isapprox(x::Vector{T}, y::Vector{S}; rtol::Real=rtoldefault(T,S,0), atol::Real=0.0,
        nans::Bool=false) where {T<:AbstractSeries,S<:AbstractSeries}

    x == y || norm(x-y,1) <= atol + rtol*max(norm(x,1), norm(y,1)) ||
        (nans && isnan(x) && isnan(y))
end


#taylor_expand function for Taylor1
"""
    taylor_expand(f, x0; order)

Computes the Taylor expansion of the function `f` around the point `x0`.

If `x0` is a scalar, a `Taylor1` expansion will be returned. If `x0` is a vector,
a `TaylorN` expansion will be computed. If the dimension of x0 (`length(x0)`)
is different from the variables set for `TaylorN` (`get_numvars()`), an
`AssertionError` will be thrown.
"""
function taylor_expand(f::F; order::Int=15) where {F}
   a = Taylor1(order)
   return f(a)
end

function taylor_expand(f::F, x0::T; order::Int=15) where {T<:Number, F}
   a = Taylor1(T[x0, one(x0)], order)
   return f(a)
end

#taylor_expand function for TaylorN
function taylor_expand(f::F, x0::Vector{T}; order::Int=TS.order()) where {T<:Number, F}
    ll = length(x0)
    @assert ll == get_numvars() && order <= TS.order()
    X = Array{TaylorN{T}}(undef, ll)
    for i in eachindex(X)
        X[i] = x0[i] + TaylorN(T, i, order=order)
    end
    return f( X )
end

function taylor_expand(f::F, x1::Vararg{Number,N}; order::Int=TS.order()) where {F, N}
    x0 = promote(x1...)
    ll = length(x0)
    T = eltype(x0[1])
    @assert ll == get_numvars() && order <= TS.order()
    X = Array{TaylorN{T}}(undef, ll)
    for i in eachindex(X)
        X[i] = x0[i] + TaylorN(T, i, order=order)
    end
    return f( X... )
end


#update! function for Taylor1
"""
    update!(a, x0)

Takes `a <: Union{Taylo1,TaylorN}` and expands it around the coordinate `x0`.
"""
function update!(a::Taylor1{T}, x0::T) where {T<:Number}
    a.coeffs .= evaluate(a, Taylor1([x0, one(x0)], order(a)) ).coeffs
    return nothing
end
function update!(a::Taylor1{T}, x0::S) where {T<:Number, S<:Number}
    xx0 = convert(T, x0)
    return update!(a, xx0)
end

#update! function for TaylorN
function update!(a::TaylorN{T}, vals::Vector{T}) where {T<:Number}
    a.coeffs .= evaluate(a, variables(a.space; order=order(a)) .+ vals).coeffs
    return nothing
end
function update!(a::TaylorN{T}, vals::Vector{S}) where {T<:Number, S<:Number}
    vv = convert(Vector{T}, vals)
    return update!(a, vv)
end

function update!(a::Union{Taylor1,TaylorN})
    #shifting around zero shouldn't change anything...
    nothing
end


for T in (:Taylor1, :TaylorN)
    @eval deg2rad(z::$T{T}) where {T<:AbstractFloat} = z * (convert(T, pi) / 180)
    @eval deg2rad(z::$T{T}) where {T<:Real} = z * (convert(float(T), pi) / 180)

    @eval rad2deg(z::$T{T}) where {T<:AbstractFloat} = z * (180 / convert(T, pi))
    @eval rad2deg(z::$T{T}) where {T<:Real} = z * (180 / convert(float(T), pi))
end

# Internal mutating deg2rad!, rad2deg! functions
for T in (:Taylor1, :TaylorN)
    @eval @inline function deg2rad!(v::$T{T}, a::$T{T}, k::Int) where {T<:AbstractFloat}
        @inbounds v[k] = a[k] * (convert(T, pi) / 180)
        return nothing
    end
    @eval @inline function deg2rad!(v::$T{S}, a::$T{T}, k::Int) where {S<:AbstractFloat,T<:Real}
        @inbounds v[k] = a[k] * (convert(float(T), pi) / 180)
        return nothing
    end
    @eval @inline function rad2deg!(v::$T{T}, a::$T{T}, k::Int) where {T<:AbstractFloat}
        @inbounds v[k] = a[k] * (180 / convert(T, pi))
        return nothing
    end
    @eval @inline function rad2deg!(v::$T{S}, a::$T{T}, k::Int) where {S<:AbstractFloat,T<:Real}
        @inbounds v[k] = a[k] * (180 / convert(float(T), pi))
        return nothing
    end
end


# In-place abs and abs2
for T in (:Taylor1, :TaylorN)
    @eval begin
        function abs!(c::$T{T}, a::$T{T}, k::Int) where {T<:Number}
            z = zero(constant_term(a))
            if constant_term(a) > constant_term(z)
                # return add!(c, a, k)
                return identity!(c, a, k)
            elseif constant_term(a) < constant_term(z)
                return subst!(c, a, k)
            else
                throw(DomainError(a,
                    """The 0th order coefficient must be non-zero
                    (abs(x) is not differentiable at x=0)."""))
            end
            return nothing
        end

    end
end

abs2!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where {T<:Number} =
    sqr!(c, a, zero(a[0]), k)
abs2!(c::TaylorN{T}, a::TaylorN{T}, k::Int) where {T<:Number} =
    sqr!(c, a, zero(a[0][1]), k)


# Inverse of a series and inverse map
@doc doc"""
    inverse(f)

Return the Taylor expansion of ``f^{-1}(t)``, of order `N = order(f)`,
for `f::Taylor1` polynomial, assuming the first coefficient of `f` is zero.
Otherwise, a `DomainError` is thrown.

The algorithm implements Lagrange inversion at ``t=0`` if ``f(0)=0``:
```math
\begin{equation*}
f^{-1}(t) = \sum_{n=1}^{N} \frac{t^n}{n!} \left.
    \frac{{\rm d}^{n-1}}{{\rm d} z^{n-1}}\left(\frac{z}{f(z)}\right)^n
    \right\vert_{z=0}.
\end{equation*}
```

""" inverse

function inverse(f::Taylor1{T}) where {T<:Number}
    if !iszero(f[0])
        throw(DomainError(f,
        """
        Evaluation of Taylor1 series at 0 is non-zero; revert
        a Taylor1 series with constant coefficient 0 and re-expand about f(0).
        """))
    end
    z = Taylor1(T, order(f))
    zdivf = z/f
    zdivfpown = zdivf
    res = Taylor1(zero(TS.numtype(zdivf)), order(f))

    @inbounds for ord in 1:order(f)
        res[ord] = zdivfpown[ord-1]/ord
        zdivfpown *= zdivf
    end
    return res
end


@doc doc"""
    inverse_map(f)

Return the Taylor expansion of ``f^{-1}(t)``, of order `N = order(f)`,
for `Taylor1` or `TaylorN` polynomials, assuming the first coefficient of `f` is zero.
Otherwise, a `DomainError` is thrown.

This method is based in the algorithm by M. Berz, Modern map methods in
Particle Beam Physics, Academic Press (1999), Sect 2.3.1.
See [`inverse`](@ref) (for `f::Taylor1`).

"""
function inverse_map(p::Taylor1)
    if !iszero(constant_term(p))
        throw(DomainError(p,
        """
        Evaluation of Taylor1 series at 0 is non-zero; revert
        a Taylor1 series with constant coefficient 0 and re-expand about f(0).
        """))
    end
    inv_m_pol = inv(linear_polynomial(p)[1])
    n_pol = inv_m_pol * nonlinear_polynomial(p)
    scaled_ident = inv_m_pol * Taylor1(order(p))
    res = scaled_ident
    aux1 = zero(res)
    aux2 = zero(res)
    for ord in 1:order(p)
        _horner!(aux2, n_pol, res, aux1)
        subst!(res, scaled_ident, aux2, ord)
    end
    return res
end

function inverse_map(p::Vector{TaylorN{T}}) where {T<:NumberNotSeries}
    if !iszero(constant_term(p))
        throw(DomainError(p,
        """
        Evaluation of Taylor1 series at 0 is non-zero; revert
        a Taylor1 series with constant coefficient 0 and re-expand about f(0).
        """))
    end
    @assert length(p) == get_numvars(p[1])
    _check_same_space_all(p[1], p)
    inv_m_pol = inv(jacobian(p))
    n_pol = inv_m_pol * nonlinear_polynomial(p)
    scaled_ident = inv_m_pol * TaylorN.(Ref(p[1].space), 1:get_numvars(p[1]), order=order(p[1]))
    res = deepcopy(scaled_ident)
    aux = zero.(res)
    auxvec = [zero(res[1]) for val in eachindex(res[1])]
    valscache = [zero(val) for val in res]
    aaux = zero(res[1])
    for ord in 1:order(p[1])
        t_res = (res...,)
        for i = 1:get_numvars(p[1])
            zero!.(auxvec)
            _evaluate!(auxvec, n_pol[i], t_res, valscache, aaux)
            aux[i] = sum( auxvec )
            subst!(res[i], scaled_ident[i], aux[i], ord)
        end
    end
    return res
end