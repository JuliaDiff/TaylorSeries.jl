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
    @inbounds coeffs[1] = a.coeffs[2]
    @inbounds for i = 1:a.order
        coeffs[i] = i*a.coeffs[i+1]
    end
    return Taylor1(coeffs, a.order)
end

"""
    derivative(n, a)

Return the value of the `n`-th derivative of the polynomial `a`.
"""
function derivative{T<:Number}(n::Int, a::Taylor1{T})
    @assert a.order >= n >= 0
    factorial( widen(n) ) * a.coeffs[n+1] :: T
end

## Integrating ##
"""
    integrate(a, [x])

Return the integral of `a::Taylor1`. The constant of integration
(0-th order coefficient) is set to `x`, which is zero if ommitted.
"""
function integrate{T<:Number, S<:Number}(a::Taylor1{T}, x::S)
    R = promote_type(T, typeof(a.coeffs[1] / 1), S)
    coeffs = zeros(R, a.order+1)
    @inbounds for i = 1:a.order
        coeffs[i+1] = a.coeffs[i] / i
    end
    @inbounds coeffs[1] = convert(R, x)
    return Taylor1(coeffs, a.order)
end
integrate{T<:Number}(a::Taylor1{T}) = integrate(a, zero(T))



## Differentiation ##
"""
    derivative(a, r)

Partial differentiation of `a::HomogeneousPolynomial` series with respect
to the `r`-th variable.
"""
function derivative(a::HomogeneousPolynomial, r::Int)
    @assert 1 <= r <= get_numvars()
    T = eltype(a)
    a.order == 0 && return HomogeneousPolynomial(zero(T))
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
        coeffs[pos] = n * a.coeffs[i]
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

    @inbounds for ord in eachindex(coeffs)
        ord == a.order+1 && continue
        coeffs[ord] = derivative( a.coeffs[ord+1], r)
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
function jacobian{T<:Number}(vf::Array{TaylorN{T},1})
    numVars = get_numvars()
    @assert length(vf) == numVars
    jac = Array{T}(numVars,numVars)

    @inbounds for comp = 1:numVars
        jac[:,comp] = vf[comp].coeffs[2].coeffs[1:end]
    end

    return transpose(jac)
end
function jacobian{T<:Number,S<:Number}(vf::Array{TaylorN{T},1},vals::Array{S,1})
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

function jacobian{T<:Number}(vf::Array{Taylor1{TaylorN{T}},1})
    vv = convert(Array{TaylorN{Taylor1{T}},1}, vf)
    jacobian(vv)
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
hessian{T<:Number,S<:Number}(f::TaylorN{T}, vals::Array{S,1}) =
    (R = promote_type(T,S); jacobian( gradient(f), vals::Array{R,1}) )
hessian{T<:Number}(f::TaylorN{T}) = hessian( f, zeros(T, get_numvars()) )

## TODO: Integration...
