# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#


## Square root ##
function sqrt(a::Taylor1{T}) where {T<:Number}
    # First non-zero coefficient
    l0nz = findfirst(a)
    aaux = zero(sqrt( constant_term(a) ))
    if l0nz < 0
        return Taylor1(aaux, order(a))
    elseif isodd(l0nz) # l0nz must be pair
        throw(DomainError(a,
            """First non-vanishing Taylor1 coefficient must correspond
            to an **even power** in order to expand `sqrt` around 0."""))
    end
    # The last l0nz coefficients are dropped.
    lnull = l0nz >> 1 # integer division by 2
    c_order = l0nz == 0 ? order(a) : order(a) >> 1
    c = Taylor1( aaux, c_order )
    aa = convert(Taylor1{eltype(aaux)}, a)
    aux = zero(aa)
    for k in eachindex(c)
        sqrt!(c, aa, aux, k, lnull)
    end
    return c
end

function sqrt(a::TaylorN{T}) where {T<:Number}
    p0 = sqrt( constant_term(a) )
    if TS._isthinzero(p0)
        throw(DomainError(a,
            """The 0-th order TaylorN coefficient must be non-zero
            in order to expand `sqrt` around 0."""))
    end
    c = TaylorN(a.space, p0, order(a))
    aa = convert(TaylorN{eltype(p0)}, a)
    aux = zero(aa)
    for k in eachindex(c)
        sqrt!(c, aa, aux, k)
    end
    return c
end

function sqrt(a::Taylor1{TaylorN{T}}) where {T<:NumberNotSeries}
    # First non-zero coefficient
    l0nz = findfirst(a)
    aux = zero(a)
    if l0nz < 0
        return Taylor1( aux.coeffs[1], order(a) )
    elseif isodd(l0nz) # l0nz must be pair
        throw(DomainError(a,
            """First non-vanishing Taylor1 coefficient must correspond
            to an **even power** in order to expand `sqrt` around 0."""))
    end
    # The last l0nz coefficients are dropped.
    lnull = l0nz >> 1 # integer division by 2
    c_order = l0nz == 0 ? order(a) : order(a) >> 1
    c = Taylor1( aux.coeffs[1], c_order )
    aa = convert(Taylor1{eltype(aux.coeffs[1])}, a)
    for k in eachindex(c)
        sqrt!(c, aa, aux, k, lnull)
    end
    return c
end


# Homogeneous coefficients for the square-root
@doc doc"""
    sqrt!(c, a, aux, k::Int, k0::Int=0)

Compute the `k-th` expansion coefficient `c[k]` of `c = sqrt(a)`
for both`c` and `a` either `Taylor1` or `TaylorN`.

The coefficients are given by

```math
\begin{aligned}
c_k &= \frac{1}{2 c_0} \big( a_k - 2 \sum_{j=1}^{(k-1)/2} c_{k-j}c_j\big),
    \text{ if $k$ is odd,} \\
c_k &= \frac{1}{2 c_0} \big( a_k - 2 \sum_{j=1}^{(k-2)/2} c_{k-j}c_j
    - (c_{k/2})^2\big), \text{ if $k$ is even.}
\end{aligned}
```

For `Taylor1` polynomials, `k0` is the order of the first non-zero
coefficient, which must be even.

""" sqrt!

function sqrt!(c::Taylor1{T}, a::Taylor1{T}, ::Taylor1{T}, k::Int,
        k0::Int=0) where {T<:NumberNotSeries}
    k < k0 && return nothing
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    if k == k0
        @inbounds c_coeffs[k+1] = sqrt(a_coeffs[2*k0+1])
        return nothing
    end
    # Recursion formula
    kodd = (k - k0)%2
    # kend = div(k - k0 - 2 + kodd, 2)
    kend = (k - k0 - 2 + kodd) >> 1
    imax = min(k0+kend, order(a))
    imin = max(k0+1, k+k0-order(a))
    kk = k+1
    @inbounds acc = zero(c_coeffs[kk])
    if k+k0 ≤ order(a)
        @inbounds acc = a_coeffs[kk+k0]
    end
    if kodd == 0
        @inbounds acc -= (c_coeffs[kend+k0+2])^2
    end
    imin ≤ imax && ( @inbounds acc -= 2 * c_coeffs[imin+1] * c_coeffs[kk+k0-imin] )
    @inbounds for i = imin+1:imax
        acc -= 2 * c_coeffs[i+1] * c_coeffs[kk+k0-i]
    end
    @inbounds c_coeffs[kk] = acc / (2*c_coeffs[k0+1])
    return nothing
end

function sqrt!(c::TaylorN{T}, a::TaylorN{T}, ::TaylorN{T}, k::Int) where
        {T<:NumberNotSeriesN}
    if k == 0
        @inbounds c.coeffs[1].coeffs[1] = sqrt( constant_term(a) )
        return nothing
    end
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    kk = k + 1
    # Recursion formula
    kodd = k%2
    kend = (k - 2 + kodd) >> 1
    # c[k] <- a[k]
    @inbounds for i in eachindex(c_coeffs[kk])
        c_coeffs[kk].coeffs[i] = a_coeffs[kk].coeffs[i]
    end
    if kodd == 0
        # @inbounds c[k] <- c[k] - (c[kend+1])^2
        @inbounds mul_scalar!(c_coeffs[kk], -1, c_coeffs[kend+2], c_coeffs[kend+2])
    end
    @inbounds for i = 1:kend
        # c[k] <- c[k] - 2*c[i]*c[k-i]
        mul_scalar!(c_coeffs[kk], -2, c_coeffs[i+1], c_coeffs[kk-i])
    end
    # @inbounds c[k] <- c[k] / (2*c[0])
    div!(c_coeffs[kk], c_coeffs[kk], 2*constant_term(c))
    return nothing
end

function sqrt!(c::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}}, aux0::Taylor1{TaylorN{T}},
        k::Int, k0::Int=0) where {T<:NumberNotSeries}
    k < k0 && return nothing
    if k == k0
        @inbounds for l in eachindex(c.coeffs[k+1])
            sqrt!(c.coeffs[k+1], a.coeffs[2*k0+1], aux0.coeffs[k+1], l)
        end
        return nothing
    end
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    kk = k + 1
    # Recursion formula
    kodd = (k - k0)%2
    # kend = div(k - k0 - 2 + kodd, 2)
    kend = (k - k0 - 2 + kodd) >> 1
    imax = min(k0+kend, order(a))
    imin = max(k0+1, k+k0-order(a))
    if k+k0 ≤ order(a)
        # @inbounds c[k] += a[k+k0]
        ### TODO: add in-place add! method for Taylor1, TaylorN and mixtures: c[k] += a[k] -> add!(c, a, k)
        ###       and/or add identity! method such that each coeff is copied individually,
        ###       otherwise memory-mixing issues happen
        identity!(c_coeffs[kk], a_coeffs[kk+k0])
    end
    if kodd == 0
        # c[k] <- c[k] - c[kend+1]^2
        # TODO: use accsqr! here?
        @inbounds mul_scalar!(c_coeffs[kk], -1, c_coeffs[kend+k0+2], c_coeffs[kend+k0+2])
    end
    @inbounds for i = imin:imax
        # c[k] <- c[k] - 2 * c[i] * c[k+k0-i]
        mul_scalar!(c_coeffs[kk], -2, c_coeffs[i+1], c_coeffs[kk+k0-i])
    end
    # @inbounds c[k] <- c[k] / (2*c[k0])
    @inbounds div_scalar!(c_coeffs[kk], 0.5, c_coeffs[k0+1])

    return nothing
end

function sqrt!(c::Taylor1{Taylor1{T}}, a::Taylor1{Taylor1{T}}, aux::Taylor1{Taylor1{T}},
        k::Int, k0::Int=0) where {T<:Number}
    k < k0 && return nothing
    c_coeffs = c.coeffs
    a_coeffs = a.coeffs
    if k == k0
        # Compute sqrt(a[2*k0])
        # First non-zero coefficient
        dosk01 = 2*k0 + 1
        l0nz = findfirst(a_coeffs[dosk01])
        _aaux = zero(sqrt( constant_term(a_coeffs[dosk01]) ))
        if l0nz < 0
            return Taylor1(_aaux, order(a_coeffs[dosk01]))
        elseif isodd(l0nz) # l0nz must be pair
            throw(DomainError(a_coeffs[dosk01],
                """First non-vanishing Taylor1 coefficient must correspond
                to an **even power** in order to expand `sqrt` around 0."""))
        end
        # The last l0nz coefficients are dropped.
        lnull = l0nz >> 1 # integer division by 2
        # aux.coeffs[1] = zero(a_coeffs[dosk01])
        zero!(aux.coeffs[1], 0)
        for j in eachindex(c_coeffs[k+1])
            sqrt!(c_coeffs[k+1], a_coeffs[dosk01], aux.coeffs[1], j, lnull)
        end
        return nothing
    end
    # Recursion formula
    kodd = (k - k0)%2
    kend = (k - k0 - 2 + kodd) >> 1
    imax = min(k0+kend, order(a))
    imin = max(k0+1, k+k0-order(a))
    kk = k + 1
    if k+k0 ≤ order(a)
        # @inbounds c[k] = a[k+k0]
        for j in eachindex(c_coeffs[kk])
            @inbounds identity!(c_coeffs[kk], a_coeffs[kk+k0], j)
        end
    end
    aux_coeffs = aux.coeffs
    zero!(aux)
    if kodd == 0
        aaux = zero(aux_coeffs[1].coeffs[1])
        # @inbounds c[k] -= (c[kend+k0+1])^2
        @inbounds for j in eachindex(c_coeffs[kk])
            sqr!(aux_coeffs[kk], c_coeffs[kend+k0+2], aaux, j)
            subst!(c_coeffs[kk], c_coeffs[kk], aux_coeffs[kk], j)
            # zero!(aaux)
        end
    end
    @inbounds for i = imin:imax
        # c[k] -= 2 * c[i] * c[k+k0-i]
        for j in eachindex(c_coeffs[kk])
            zero!(aux_coeffs[kk], j)
            mul_scalar!(aux_coeffs[kk], 2, c_coeffs[i+1], c_coeffs[kk+k0-i], j)
            subst!(c_coeffs[kk], c_coeffs[kk], aux_coeffs[kk], j)
        end
    end
    # @inbounds c[k] = c[k] / (2*c[k0])
    @inbounds for j in eachindex(c_coeffs[kk])
        identity!(aux_coeffs[kk], c_coeffs[kk], j)
    end
    @inbounds for j in eachindex(c[k0])
        div!(c_coeffs[kk], aux_coeffs[kk], c_coeffs[k0+1], j)
    end
    @inbounds for j in eachindex(c[k0])
        div!(c_coeffs[kk], c_coeffs[kk], 2, j)
    end
    return nothing
end
