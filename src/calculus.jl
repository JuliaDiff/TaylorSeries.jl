# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

## Differentiating ##
"""
    derivative(a)

Return the `Taylor1` polynomial of the differential of `a::Taylor1`.
The last coefficient is set to zero.

The function `differentiate` is an exact synonym of `derivative`.
"""
function derivative(a::Taylor1)
    res = zero(a)
    @inbounds for ord in eachindex(a)
        derivative!(res, a, ord)
        # ord == a.order-1 && break
    end
    return res
end

"""
    differentiate

An exact synonym of [`derivative`](@ref).
"""
const differentiate = derivative

"""
    derivative!(res, a) --> nothing

In-place version of `derivative`. Compute the `Taylor1` polynomial of the
differential of `a::Taylor1` and save it into `res`. The last coefficient is
set to zero.
"""
function derivative!(res::Taylor1, a::Taylor1)
    @inbounds for ord in eachindex(a)
        derivative!(res, a, ord)
        # ord == a.order-1 && break
    end
    res[a.order] = zero(a[0])
    nothing
end

"""
    derivative!(p, a, k) --> nothing

Update in-place the `k-th` expansion coefficient `p[k]` of `p = derivative(a)`
for both `p` and `a` `Taylor1`.

The coefficients are given by

```math
p_k = (k+1)a_{k+1}.
```

"""
derivative!(p::Taylor1, a::Taylor1, k::Int) = k < a.order ? p[k] = (k+1)*a[k+1] : nothing


"""
    derivative(a, n)

Compute recursively the `Taylor1` polynomial of the n-th derivative of
`a::Taylor1`.
"""
function derivative(a::Taylor1{T}, n::Int) where {T <: Number}
    @assert a.order ≥ n ≥ 0
    if n==0
        return a
    else
        res = deepcopy(a)
        for i in 1:n
            derivative!(res, res)
        end
        return res
    end
end

"""
    derivative(n, a)

Return the value of the `n`-th derivative of the polynomial `a`.
"""
function derivative(n::Int, a::Taylor1{T}) where {T<:Number}
    @assert a.order ≥ n ≥ 0
    factorial( widen(n) ) * a[n] :: T
end

## Integrating ##
"""
    integrate(a, [x])

Return the integral of `a::Taylor1`. The constant of integration
(0-th order coefficient) is set to `x`, which is zero if ommitted.
"""
function integrate(a::Taylor1{T}, x::S) where {T<:Number, S<:Number}
    order = get_order(a)
    aa = a[0]/1 + zero(x)
    R = typeof(aa)
    coeffs = Array{typeof(aa)}(undef, order+1)
    fill!(coeffs, zero(aa))
    @inbounds for i = 1:order
        coeffs[i+1] = a[i-1] / i
    end
    @inbounds coeffs[1] = convert(R, x)
    return Taylor1(coeffs, a.order)
end
integrate(a::Taylor1{T}) where {T<:Number} = integrate(a, zero(a[0]))



## Differentiation ##
"""
    derivative(a, r)

Partial differentiation of `a::HomogeneousPolynomial` series with respect
to the `r`-th variable.
"""
function derivative(a::HomogeneousPolynomial, r::Int)
    @assert 1 ≤ r ≤ get_numvars()
    T = eltype(a)
    a.order == 0 && return HomogeneousPolynomial([zero(T)], 0)
    @inbounds num_coeffs = size_table[a.order]
    coeffs = zeros(T,num_coeffs)
    @inbounds posTb = pos_table[a.order]
    @inbounds num_coeffs = size_table[a.order+1]

    @inbounds for i = 1:num_coeffs
        iind = coeff_table[a.order+1][i]
        n = iind[r]
        n == 0 && continue
        iind[r] -= 1
        kdic = in_base(get_order(),iind)
        pos = posTb[kdic]
        coeffs[pos] = n * a[i]
        iind[r] += 1
    end

    return HomogeneousPolynomial{T}(coeffs, a.order-1)
end
derivative(a::HomogeneousPolynomial, s::Symbol) = derivative(a, lookupvar(s))

"""
    derivative(a, r)

Partial differentiation of `a::TaylorN` series with respect
to the `r`-th variable. The `r`-th variable may be also
specified through its symbol.
"""
function derivative(a::TaylorN, r=1::Int)
    T = eltype(a)
    coeffs = Array{HomogeneousPolynomial{T}}(undef, a.order)

    @inbounds for ord = 1:a.order
        coeffs[ord] = derivative( a[ord], r)
    end
    return TaylorN{T}( coeffs, a.order )
end
derivative(a::TaylorN, s::Symbol) = derivative(a, lookupvar(s))

"""
    derivative(a::TaylorN{T}, ntup::NTuple{N,Int})

Return a `TaylorN` with the partial derivative of `a` defined
by `ntup::NTuple{N,Int}`, where the first entry is the number
of derivatives with respect to the first variable, the second is
the number of derivatives with respect to the second, and so on.
"""
function derivative(a::TaylorN, ntup::NTuple{N,Int}) where {N}

    @assert N == get_numvars() && all(ntup .>= 0)

    sum(ntup) > a.order && return zero(a)
    sum(ntup) == 0 && return copy(a)

    aa = copy(a)
    for nvar in 1:get_numvars()
        for numder in 1:ntup[nvar]
            aa = derivative(aa, nvar)
        end
    end

    return aa
end

"""
    derivative(ntup::NTuple{N,Int}, a::TaylorN{T})

Returns the value of the coefficient of `a` specified by
`ntup::NTuple{N,Int}`, multiplied by the corresponding
factorials.
"""
function derivative(ntup::NTuple{N,Int}, a::TaylorN) where {N}

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
    T = eltype(f)
    numVars = get_numvars()
    grad = Array{TaylorN{T}}(undef, numVars)
    @inbounds for nv = 1:numVars
        grad[nv] = derivative(f, nv)
    end
    return grad
end
const ∇ = TaylorSeries.gradient

"""
```
    jacobian(vf)
    jacobian(vf, [vals])
```

Compute the jacobian matrix of `vf`, a vector of `TaylorN` polynomials,
evaluated at the vector `vals`. If `vals` is ommited, it is evaluated at zero.
"""
function jacobian(vf::Array{TaylorN{T},1}) where {T<:Number}
    numVars = get_numvars()
    @assert length(vf) == numVars
    jac = Array{T}(undef, numVars, numVars)

    @inbounds for comp = 1:numVars
        jac[:,comp] = vf[comp][1][1:end]
    end

    return transpose(jac)
end
function jacobian(vf::Array{TaylorN{T},1}, vals::Array{S,1}) where
        {T<:Number, S<:Number}

    R = promote_type(T,S)
    numVars = get_numvars()
    @assert length(vf) == numVars == length(vals)
    jac = Array{R}(undef, numVars, numVars)

    for comp = 1:numVars
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
evaluated at the vector `vals`, and write results to `jac`. If `vals` is ommited,
it is evaluated at zero.
"""
function jacobian!(jac::Array{T,2}, vf::Array{TaylorN{T},1}) where {T<:Number}
    numVars = get_numvars()
    @assert length(vf) == numVars
    @assert (numVars, numVars) == size(jac)
    for comp2 = 1:numVars
        for comp1 = 1:numVars
            @inbounds jac[comp1,comp2] = vf[comp1][1][comp2]
        end
    end
    nothing
end
function jacobian!(jac::Array{T,2}, vf::Array{TaylorN{T},1},
        vals::Array{T,1}) where {T<:Number}

    numVars = get_numvars()
    @assert length(vf) == numVars == length(vals)
    @assert (numVars, numVars) == size(jac)
    for comp = 1:numVars
        @inbounds for nv = 1:numVars
            jac[nv,comp] = evaluate(derivative(vf[nv], comp), vals)
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
evaluated at the vector `vals`. If `vals` is ommited, it is evaluated at
zero.
"""
hessian(f::TaylorN{T}, vals::Array{S,1}) where {T<:Number, S<:Number} =
    (R = promote_type(T,S); jacobian( gradient(f), vals::Array{R,1}) )

hessian(f::TaylorN{T}) where {T<:Number} = hessian( f, zeros(T, get_numvars()) )

"""
```
    hessian!(hes, f)
    hessian!(hes, f, [vals])
```

Return the hessian matrix (jacobian of the gradient) of `f::TaylorN`,
evaluated at the vector `vals`, and write results to `hes`. If `vals` is
ommited, it is evaluated at zero.
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
    T = promote_type(eltype(a), eltype(a[1]/1))
    a.order == order_max && return HomogeneousPolynomial(zero(T), 0)

    @inbounds posTb = pos_table[a.order+2]
    @inbounds num_coeffs = size_table[a.order+1]

    coeffs = zeros(T, size_table[a.order+2])

    @inbounds for i = 1:num_coeffs
        iind = coeff_table[a.order+1][i]
        n = iind[r]
        n == order_max && continue
        iind[r] += 1
        kdic = in_base(get_order(), iind)
        pos = posTb[kdic]
        coeffs[pos] = a[i] / (n+1)
        iind[r] -= 1
    end

    return HomogeneousPolynomial{T}(coeffs, a.order+1)
end
integrate(a::HomogeneousPolynomial, s::Symbol) = integrate(a, lookupvar(s))


"""
    integrate(a, r, [x0])

Integrate the `a::TaylorN` series with respect to the `r`-th variable,
where `x0` the integration constant and must be independent
of the `r`-th variable; if `x0` is ommitted, it is taken as zero.
"""
function integrate(a::TaylorN, r::Int)
    T = promote_type(eltype(a), eltype(a[0]/1))
    order_max = min(get_order(), a.order+1)
    coeffs = zeros(HomogeneousPolynomial{T}, order_max)

    @inbounds for ord = 0:order_max-1
        coeffs[ord+1] = integrate( a[ord], r)
    end

    return TaylorN(coeffs)
end
function integrate(a::TaylorN, r::Int, x0::TaylorN)
    # Check constant of integration is independent of re
    @assert derivative(x0, r) == 0.0 """
    The integration constant ($x0) must be independent of the
    $(_params_TaylorN_.variable_names[r]) variable"""

    res = integrate(a, r)
    return x0+res
end
integrate(a::TaylorN, r::Int, x0::NumberNotSeries) =
    integrate(a,r,TaylorN(HomogeneousPolynomial([convert(eltype(a),x0)], 0)))

integrate(a::TaylorN, s::Symbol) = integrate(a, lookupvar(s))
integrate(a::TaylorN, s::Symbol, x0::TaylorN) = integrate(a, lookupvar(s), x0)
integrate(a::TaylorN, s::Symbol, x0::NumberNotSeries) = integrate(a, lookupvar(s), x0)
