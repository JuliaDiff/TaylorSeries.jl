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
"""
function derivative(a::Taylor1)
    coeffs = zero(a.coeffs)
    @inbounds coeffs[1] = a[1]
    @inbounds for i = 1:a.order
        coeffs[i] = i*a[i]
    end
    return Taylor1(coeffs, a.order)
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
    R = promote_type(T, typeof(a[0] / 1), S)
    coeffs = zeros(R, a.order+1)
    @inbounds for i = 1:a.order
        coeffs[i+1] = a[i-1] / i
    end
    @inbounds coeffs[1] = convert(R, x)
    return Taylor1(coeffs, a.order)
end
integrate(a::Taylor1{T}) where {T<:Number} = integrate(a, zero(T))



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

"""
    derivative(a, [r=1])

Partial differentiation of `a::TaylorN` series with respect
to the `r`-th variable.
"""
function derivative(a::TaylorN, r=1::Int)
    T = eltype(a)
    coeffs = Array{HomogeneousPolynomial{T}}(a.order)

    @inbounds for ord = 0:a.order
        coeffs[ord+1] = derivative( a[ord], r)
    end
    return TaylorN{T}( coeffs, a.order )
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
    grad = Array{TaylorN{T}}(numVars)
    @inbounds for nv = 1:numVars
        grad[nv] = derivative(f, nv)
    end
    return grad
end
const ∇ = gradient

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
    jac = Array{T}(numVars,numVars)

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
    jac = Array{R}(numVars,numVars)

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


## TODO: Integration...
