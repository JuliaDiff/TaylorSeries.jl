# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

## Differentiating ##
"""
    differentiate(a)

Return the `Taylor1` polynomial of the differential of `a::Taylor1`.
The order of the result is `a.order-1`.

The function `derivative` is an exact synonym of `differentiate`.
"""
function differentiate(a::Taylor1)
    res = Taylor1(zero(a[0]), get_order(a)-1)
    for ord in eachindex(res)
        differentiate!(res, a, ord)
    end
    return res
end

"""
    derivative

An exact synonym of [`differentiate`](@ref).
"""
const derivative = differentiate

"""
    differentiate!(res, a) --> nothing

In-place version of `differentiate`. Compute the `Taylor1` polynomial of the
differential of `a::Taylor1` and return it as `res` (order of `res` remains
unchanged).
"""
function differentiate!(res::Taylor1, a::Taylor1)
    for ord in eachindex(res)
        differentiate!(res, a, ord)
    end
    return nothing
end

"""
    differentiate!(p, a, k) --> nothing

Update in-place the `k-th` expansion coefficient `p[k]` of `p = differentiate(a)`
for both `p` and `a` `Taylor1`.

The coefficients are given by

```math
p_k = (k+1) a_{k+1}.
```

"""
function differentiate!(p::Taylor1{T}, a::Taylor1{T}, k::Int) where
        {T<:NumberNotSeries}
    k >= a.order && return nothing
    @inbounds p[k] = (k+1)*a[k+1]
    return nothing
end
function differentiate!(p::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}}, k::Int) where
        {T<:NumberNotSeries}
    k >= a.order && return nothing
    @inbounds for ord in eachindex(p[k])
        zero!(p[k], ord)
        mul!(p[k], a[k+1], k+1, ord)
    end
    return nothing
end
function differentiate!(p::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}}, k::Int) where
        {T<:NumberNotSeriesN}
    k >= a.order && return nothing
    # @inbounds p[k] = (k+1)*a[k+1]
    @inbounds for ord in eachindex(p[k])
        mul!(p[k], a[k+1], k+1, ord)
    end
    return nothing
end


"""
    differentiate(a, n)

Compute recursively the `Taylor1` polynomial of the n-th derivative of
`a::Taylor1`. The order of the result is `a.order-n`.
"""
function differentiate(a::Taylor1{T}, n::Int) where {T <: Number}
    if n > a.order
        return Taylor1(zero(T), 0)
    elseif n == a.order
        return Taylor1(differentiate(n, a), 0)
    elseif n==0
        return a
    end
    res = differentiate(a)
    for i = 2:n
        differentiate!(res, res)
    end
    return Taylor1(res.coeffs[1:a.order-n+1])
end

"""
    differentiate(n, a)

Return the value of the `n`-th differentiate of the polynomial `a`.
"""
function differentiate(n::Int, a::Taylor1{T}) where {T<:Number}
    @assert a.order ≥ n ≥ 0
    return factorial( widen(n) ) * a[n] :: T
end


## Integrating ##
"""
    integrate(a, [x])

Return the integral of `a::Taylor1`. The constant of integration
(0-th order coefficient) is set to `x`, which is zero if omitted.
Note that the order of the result is `a.order+1`.
"""
function integrate(a::Taylor1{T}, x::S) where {T<:Number, S<:Number}
    order = get_order(a)
    aa = a[0]/1 + zero(x)
    R = typeof(aa)
    coeffs = Array{typeof(aa)}(undef, order+1)
    # fill!(coeffs, zero(aa))
    @inbounds for i = 1:order
        coeffs[i+1] = a[i-1] / i
    end
    @inbounds coeffs[1] = convert(R, x)
    return Taylor1(coeffs, a.order)
end
integrate(a::Taylor1{T}) where {T<:Number} = integrate(a, zero(a[0]))



## Differentiation ##
"""
    differentiate(a, r)

Partial differentiation of `a::HomogeneousPolynomial` series with respect
to the `r`-th variable.
"""
function differentiate(a::HomogeneousPolynomial, r::Int)
    @assert 1 ≤ r ≤ get_numvars()
    T = TS.numtype(a)
    a.order == 0 && return HomogeneousPolynomial(zero(a[1]), 0)
    @inbounds num_coeffs = size_table[a.order]
    coeffs = zeros(T, num_coeffs)
    @inbounds posTb = pos_table[a.order]
    @inbounds num_coeffs = size_table[a.order+1]
    ct = deepcopy(coeff_table[a.order+1])
    @inbounds for i = 1:num_coeffs
        # iind = @isonethread coeff_table[a.order+1][i]
        iind = ct[i]
        n = iind[r]
        n == 0 && continue
        iind[r] -= 1
        kdic = in_base(get_order(), iind)
        pos = posTb[kdic]
        coeffs[pos] = n * a[i]
        iind[r] += 1
    end

    return HomogeneousPolynomial{T}(coeffs, a.order-1)
end
differentiate(a::HomogeneousPolynomial, s::Symbol) = differentiate(a, lookupvar(s))

"""
    differentiate(a, r)

Partial differentiation of `a::TaylorN` series with respect
to the `r`-th variable. The `r`-th variable may be also
specified through its symbol.
"""
function differentiate(a::TaylorN, r=1::Int)
    T = TS.numtype(a)
    coeffs = Array{HomogeneousPolynomial{T}}(undef, a.order)

    @inbounds for ord = 1:a.order
        coeffs[ord] = differentiate( a[ord], r)
    end
    return TaylorN{T}( coeffs, a.order )
end
differentiate(a::TaylorN, s::Symbol) = differentiate(a, lookupvar(s))

"""
    differentiate(a::TaylorN{T}, ntup::NTuple{N,Int})

Return a `TaylorN` with the partial derivative of `a` defined
by `ntup::NTuple{N,Int}`, where the first entry is the number
of derivatives with respect to the first variable, the second is
the number of derivatives with respect to the second, and so on.
"""
function differentiate(a::TaylorN, ntup::NTuple{N,Int}) where {N}

    @assert N == get_numvars() && all(ntup .>= 0)

    sum(ntup) > a.order && return zero(a)
    sum(ntup) == 0 && return copy(a)

    aa = copy(a)
    for nvar in 1:get_numvars()
        for _ in 1:ntup[nvar]
            aa = differentiate(aa, nvar)
        end
    end

    return aa
end

"""
    differentiate(ntup::NTuple{N,Int}, a::TaylorN{T})

Returns the value of the coefficient of `a` specified by
`ntup::NTuple{N,Int}`, multiplied by the corresponding
factorials.
"""
function differentiate(ntup::NTuple{N,Int}, a::TaylorN) where {N}

    @assert N == get_numvars() && all(ntup .>= 0)

    c = getcoeff(a, [ntup...])
    for ind = 1:get_numvars()
        c *= factorial(ntup[ind])
    end

    return c
end


## Gradient, jacobian and hessian
"""
```
    gradient(f)
    ∇(f)
```

Compute the gradient of the polynomial `f::TaylorN`.
"""
function gradient(f::TaylorN)
    T = TS.numtype(f)
    numVars = get_numvars()
    grad = Array{TaylorN{T}}(undef, numVars)
    @inbounds for nv = 1:numVars
        grad[nv] = differentiate(f, nv)
    end
    return grad
end
const ∇ = TS.gradient

"""
```
    jacobian(vf)
    jacobian(vf, [vals])
```

Compute the jacobian matrix of `vf`, a vector of `TaylorN` polynomials,
evaluated at the vector `vals`. If `vals` is omitted, it is evaluated at zero.
"""
function jacobian(vf::Array{TaylorN{T},1}) where {T<:Number}
    numVars = get_numvars()
    jac = Array{T}(undef, numVars, length(vf))

    @inbounds for comp in eachindex(vf)
        jac[:,comp] = vf[comp][1][1:end]
    end

    return transpose(jac)
end
function jacobian(vf::Array{TaylorN{T},1}, vals::Array{S,1}) where {T<:Number,S<:Number}
    R = promote_type(T,S)
    numVars = get_numvars()
    @assert numVars == length(vals)
    jac = Array{R}(undef, numVars, length(vf))

    for comp in eachindex(vf)
        @inbounds grad = gradient( vf[comp] )
        @inbounds for nv = 1:numVars
            jac[nv,comp] = evaluate(grad[nv], vals)
        end
    end

    return transpose(jac)
end

function jacobian(vf::Array{Taylor1{TaylorN{T}},1}) where {T<:Number}
    vv = convert(Array{TaylorN{Taylor1{T}},1}, vf)
    jacobian(vv)
end

"""
```
    jacobian!(jac, vf)
    jacobian!(jac, vf, [vals])
```

Compute the jacobian matrix of `vf`, a vector of `TaylorN` polynomials
evaluated at the vector `vals`, and write results to `jac`. If `vals` is omitted,
it is evaluated at zero.
"""
function jacobian!(jac::Array{T,2}, vf::Array{TaylorN{T},1}) where {T<:Number}
    numVars = get_numvars()
    @assert (length(vf), numVars) == size(jac)
    for comp2 = 1:numVars
        for comp1 in eachindex(vf)
            @inbounds jac[comp1,comp2] = vf[comp1][1][comp2]
        end
    end
    nothing
end
function jacobian!(jac::Array{T,2}, vf::Array{TaylorN{T},1},
        vals::Array{T,1}) where {T<:Number}
    numVars = get_numvars()
    @assert numVars == length(vals)
    @assert (length(vf), numVars) == size(jac)
    for comp = 1:numVars
        @inbounds for nv in eachindex(vf)
            jac[nv,comp] = evaluate(differentiate(vf[nv], comp), vals)
        end
    end
    nothing
end


"""
```
    hessian(f)
    hessian(f, [vals])
```

Return the hessian matrix (jacobian of the gradient) of `f::TaylorN`,
evaluated at the vector `vals`. If `vals` is omitted, it is evaluated at
zero.
"""
hessian(f::TaylorN{T}, vals::Array{S,1}) where {T<:Number,S<:Number} =
    (R = promote_type(T,S); jacobian( gradient(f), vals::Array{R,1}) )

hessian(f::TaylorN{T}) where {T<:Number} = hessian( f, zeros(T, get_numvars()) )

"""
```
    hessian!(hes, f)
    hessian!(hes, f, [vals])
```

Return the hessian matrix (jacobian of the gradient) of `f::TaylorN`,
evaluated at the vector `vals`, and write results to `hes`. If `vals` is
omitted, it is evaluated at zero.
"""
hessian!(hes::Array{T,2}, f::TaylorN{T}, vals::Array{T,1}) where {T<:Number} =
    jacobian!(hes, gradient(f), vals)

hessian!(hes::Array{T,2}, f::TaylorN{T}) where {T<:Number} =
    jacobian!(hes, gradient(f))


##Integration
"""
    integrate(a, r)

Integrate the `a::HomogeneousPolynomial` with respect to the `r`-th
variable. The returned `HomogeneousPolynomial` has no added constant of
integration. If the order of a corresponds to `get_order()`, a zero
`HomogeneousPolynomial` of 0-th order is returned.

"""
function integrate(a::HomogeneousPolynomial, r::Int)
    @assert 1 ≤ r ≤ get_numvars()

    order_max = get_order()
    # NOTE: the following returns order 0, but could be get_order(), or get_order(a)
    a.order == order_max && return HomogeneousPolynomial(zero(a[1]/1), 0)

    @inbounds posTb = pos_table[a.order+2]
    @inbounds num_coeffs = size_table[a.order+1]

    T = promote_type(TS.numtype(a), TS.numtype(a[1]/1))
    coeffs = zeros(T, size_table[a.order+2])
    ct = deepcopy(coeff_table[a.order+1])
    @inbounds for i = 1:num_coeffs
        # iind = @isonethread coeff_table[a.order+1][i]
        iind = ct[i]
        n = iind[r]
        n == order_max && continue
        iind[r] += 1
        kdic = in_base(get_order(), iind)
        pos = posTb[kdic]
        coeffs[pos] = a[i] / (n+1)
        iind[r] -= 1
    end

    return HomogeneousPolynomial(coeffs, a.order+1)
end
integrate(a::HomogeneousPolynomial, s::Symbol) = integrate(a, lookupvar(s))


"""
    integrate(a, r, [x0])

Integrate the `a::TaylorN` series with respect to the `r`-th variable,
where `x0` the integration constant and must be independent
of the `r`-th variable; if `x0` is omitted, it is taken as zero.
"""
function integrate(a::TaylorN, r::Int)
    T = promote_type(TS.numtype(a), TS.numtype(a[0]/1))
    order_max = min(get_order(), a.order+1)
    coeffs = zeros(HomogeneousPolynomial{T}, order_max)

    @inbounds for ord = 0:order_max-1
        coeffs[ord+1] = integrate( a[ord], r)
    end

    return TaylorN(coeffs)
end
function integrate(a::TaylorN, r::Int, x0::TaylorN)
    # Check constant of integration is independent of re
    @assert differentiate(x0, r) == 0.0 """
    The integration constant ($x0) must be independent of the
    $(_params_TaylorN_.variable_names[r]) variable"""

    res = integrate(a, r)
    return x0+res
end
integrate(a::TaylorN, r::Int, x0) =
    integrate(a,r,TaylorN(HomogeneousPolynomial([convert(TS.numtype(a),x0)], 0)))

integrate(a::TaylorN, s::Symbol) = integrate(a, lookupvar(s))
integrate(a::TaylorN, s::Symbol, x0::TaylorN) = integrate(a, lookupvar(s), x0)
integrate(a::TaylorN, s::Symbol, x0) = integrate(a, lookupvar(s), x0)
