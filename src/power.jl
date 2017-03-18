# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#


function ^(a::HomogeneousPolynomial, n::Integer)
    n == 0 && return one(a)
    n == 1 && return a
    n == 2 && return square(a)
    n < 0 && throw(DomainError())
    return Base.power_by_squaring(a, n)
end

for T in (:Taylor1, :TaylorN)
    @eval function ^{T<:Number}(a::$T{T}, n::Integer)
        n == 0 && return one(a)
        n == 1 && return a
        n == 2 && return square(a)
        n < 0 && return inv( a^(-n) )
        return Base.power_by_squaring(a, n)
    end

    @eval function ^{T<:Integer}(a::$T{T}, n::Integer)
        n == 0 && return one(a)
        n == 1 && return a
        n == 2 && return square(a)
        n < 0 && throw(DomainError())
        return Base.power_by_squaring(a, n)
    end

    @eval ^(a::$T, x::Rational) = a^(x.num/x.den)

    @eval ^(a::$T, b::$T) = exp( b*log(a) )

    @eval ^{T<:Complex}(a::$T, x::T) = exp( x*log(a) )
end


## Real power ##
function ^{S<:Real}(a::Taylor1, r::S)
    r == zero(r) && return one(a)
    r == one(r)/2 && return sqrt(a)
    isinteger(r) && return a^round(Int,r)

    # First non-zero coefficient
    l0nz = findfirst(a)
    l0nz > a.order && return zero(a)

    # The first non-zero coefficient of the result; must be integer
    !isinteger(r*l0nz) && throw(ArgumentError(
        """The 0th order Taylor1 coefficient must be non-zero
        to raise the Taylor1 polynomial to a non-integer exponent."""))
    lnull = trunc(Int, r*l0nz )

    @inbounds aux = ( a[l0nz+1] )^r
    k0 = lnull+l0nz
    c = Taylor1( zero(aux), a.order)
    @inbounds c[lnull+1] = aux
    for k = k0+1:a.order
        pow!(c, a, r, k, l0nz)
    end

    return c
end

## Real power ##
function ^{S<:Real}(a::TaylorN, r::S)
    r == zero(r) && return TaylorN( one(eltype(a)), 0 )
    r == one(r)/2 && return sqrt(a)
    isinteger(r) && return a^round(Int,r)

    @inbounds a0 = a[1][1]
    @assert a0 != zero(a0)
    aux = ( a0 )^r

    c = TaylorN( aux, a.order)
    for ord in 1:a.order
        pow!(c, a, r, ord)
    end

    return c
end


# Homogeneous coefficients for real power
doc"""
    pow!(c, a, r::Real, k::Int, k0::Int=0)

Update the `k-th` expansion coefficient `c[k+1]` of `c = a^r`, for
both `c` and `a` either `Taylor1` or `TaylorN`.

The coefficients are given by

\begin{equation*}
c_k = \frac{1}{k a_0} \sum_{j=0}^{k-1} \big(r(k-j) -j\big)a_{k-j} c_j.
\end{equation*}

For `Taylor1` polynomials, `k0` is the order of the first non-zero
coefficient of `a`.

"""
function pow!{S<:Real}(c::Taylor1, a::Taylor1, r::S, k::Int, l0::Int=0)
    if k == l0
        @inbounds c[1] = ( a[l0+1] )^r
        return nothing
    end

    coefhomog = zero(eltype(c))
    for i = 0:k-l0-1
        aux = r*(k-i) - i
        @inbounds coefhomog += aux * a[k-i+1] * c[i+1]
    end
    aux = k - l0*(r+1)
    @inbounds coefhomog = coefhomog / (aux * a[l0+1])

    c[k-l0+1] = coefhomog
    return nothing
end

function pow!{S<:Real}(c::TaylorN, a::TaylorN, r::S, k::Int)
    if k == 0
        @inbounds c[1] = ( a[1][1] )^r
        return nothing
    end

    coefhomog = zero(eltype(c))
    coefhomog = HomogeneousPolynomial(zero(eltype(c)), k)
    for i = 0:k-1
        aux = r*(k-i) - i
        @inbounds coefhomog += aux * a[k-i+1] * c[i+1]
    end
    @inbounds c[k+1] = coefhomog / (k * a[1][1])

    return nothing
end



## Square ##
"""
    square(a::AbstractSeries) --> typeof(a)

Return `a^2`; see [TaylorSeries.sqr!](@ref).
""" square

for T in (:Taylor1, :TaylorN)
    @eval function square(a::$T)
        if $T == Taylor1
            c = Taylor1(Array{eltype(a)}(a.order+1))
        else
            c = TaylorN( zeros(HomogeneousPolynomial{eltype(a)}, a.order) )
        end

        @inbounds for k = 0:a.order
            sqr!(c, a, k)
        end
        return c
    end
end

function square(a::HomogeneousPolynomial)
    T = eltype(a)

    order = 2*a.order
    order > get_order() && return HomogeneousPolynomial([zero(T)], get_order())

    res = HomogeneousPolynomial([zero(T)], order)
    sqr!(res, a)
    return res
end


# Homogeneous coefficients for square
doc"""
    sqr!(c, a, k::Int) --> nothing

Update the `k-th` expansion coefficient `c[k+1]` of `c = a^2`, for
both `c` and `a` either `Taylor1` or `TaylorN`.

The coefficients are given by

\\begin{eqnarray*}
c_k & = & 2 \\sum_{j=0}^{(k-1)/2} a_{k-j} a_j,
    \\text{ if k is odd, or }  \\\\
c_k & = & 2 \\sum_{j=0}^{(k-2)/2} a_{k-j} a_j + (a_{k/2})^2,
    \\text{ if k is even. }
\\end{eqnarray*}

""" sqr!

for T = (:Taylor1, :TaylorN)
    @eval function sqr!(c::$T, a::$T, k::Int)
        if k == 0
            @inbounds c[1] = a[1]^2
            return nothing
        end

        if $T == Taylor1
            @inbounds c[k+1] = zero(c[1])
        else
            c[k+1] = HomogeneousPolynomial(zero(c[1][1]), k)
        end

        kodd = k%2
        kend = div(k - 2 + kodd, 2)
        for i = 0:kend
            @inbounds c[k+1] += a[i+1]*a[k-i+1]
        end
        c[k+1] = 2 * c[k+1]
        kodd == 1 && return nothing
        @inbounds c[k+1] += a[div(k,2)+1]^2

        return nothing
    end
end


"""
    sqr!(c, a)

Return `c = a*a` with no allocation; all parameters are `HomogeneousPolynomial`.

"""
function sqr!(c::HomogeneousPolynomial, a::HomogeneousPolynomial)
    iszero(a) && return nothing

    T = eltype(c)
    @inbounds num_coeffs_a = size_table[a.order+1]
    @inbounds num_coeffs  = size_table[c.order+1]

    @inbounds posTb = pos_table[c.order+1]

    for na = 1:num_coeffs_a
        @inbounds ca = a[na]
        ca == zero(T) && continue
        @inbounds inda = index_table[a.order+1][na]
        @inbounds pos = posTb[2*inda]
        @inbounds c[pos] += ca * ca
        @inbounds for nb = na+1:num_coeffs_a
            cb = a[nb]
            cb == zero(T) && continue
            indb = index_table[a.order+1][nb]
            pos = posTb[inda+indb]
            c[pos] += 2 * ca * cb
        end
    end

    return nothing
end

function squareHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1})
    kcoef == 0 && return ac[1]^2
    coefhomog = zero(T)
    kodd = kcoef%2
    kend = div(kcoef - 2 + kodd, 2)
    @inbounds for i = 0:kend
        coefhomog += ac[i+1]*ac[kcoef-i+1]
    end
    coefhomog = convert(T,2) * coefhomog
    if kodd == 0
        @inbounds coefhomog += ac[div(kcoef,2)+1]^2
    end
    coefhomog
end



## Square root ##
function sqrt(a::Taylor1)

    # First non-zero coefficient
    l0nz = findfirst(a)
    if l0nz > a.order
        return zero(a)
    elseif l0nz%2 == 1 # l0nz must be pair
        throw(ArgumentError(
        """First non-vanishing Taylor1 coefficient must correspond
        to an **even power** in order to expand `sqrt` around 0."""))
    end

    # The last l0nz coefficients are set to zero.
    lnull = div(l0nz, 2)
    @inbounds aux = sqrt( a[l0nz+1] )
    T = typeof(aux)

    c = Taylor1( zeros(T, a.order+1) )
    @inbounds c[lnull+1] = aux
    @inbounds for k = lnull+1:a.order-l0nz
        sqrt!(c, a, k, lnull)
    end

    return c
end

function sqrt(a::TaylorN)
    @inbounds p0 = sqrt( a[1][1] )
    if p0 == zero(p0)
        throw(ArgumentError(
        """The 0-th order TaylorN coefficient must be non-zero
        in order to expand `sqrt` around 0."""))
    end

    c = TaylorN( p0, a.order)
    @inbounds for k = 1:a.order
        sqrt!(c, a, k)
    end

    return c
end

# Homogeneous coefficients for the square-root
doc"""
    sqrt!(c, a, k::Int, k0::Int=0)

Compute the `k-th` expansion coefficient `c[k+1]` of `c = sqrt(a)`
for both`c` and `a` either `Taylor1` or `TaylorN`.

The coefficients are given by

\begin{eqnarray*}
c_k &=& \frac{1}{2 c_0} ( a_k - 2 \sum_{j=0}^{(k-1)/2} c_{k-j}c_j),
\text{ if $k$ is odd, or } \\\\
c_k &=& \frac{1}{2 c_0} ( a_k - 2 \sum_{j=0}^{(k-2)/2} c_{k-j}c_j) - (c_{k/2})^2,
\text{ if $k$ is even.}
\end{eqnarray*}

For `Taylor1` polynomials, `k0` is the order of the first non-zero
coefficient, which must be even.

"""
function sqrt!(c::Taylor1, a::Taylor1, k::Int, k0::Int=0)

    if k == k0
        @inbounds c[k+1] = sqrt(a[2*k0+1])
        return nothing
    end

    kodd = (k - k0)%2
    kend = div(k - k0 - 2 + kodd, 2)
    c[k+1] = zero(c[1])
    @inbounds for i = k0+1:k0+kend
        c[k+1] += c[i+1]*c[k+k0-i+1]
    end
    @inbounds aux = a[k+k0+1] - 2*c[k+1]
    if kodd == 0
        @inbounds aux = aux - (c[kend+k0+2])^2
    end
    @inbounds c[k+1] = aux / (2*c[k0+1])

    return nothing
end

function sqrt!(c::TaylorN, a::TaylorN, k::Int)

    if k == 0
        c[1] = sqrt(a[1][1])
        return nothing
    end

    kodd = k%2
    kend = div(k - 2 + kodd, 2)
    c[k+1] = HomogeneousPolynomial(zero(c[1][1]), k)
    @inbounds for i = 1:kend
        c[k+1] += c[i+1]*c[k-i+1]
    end
    @inbounds aux = a[k+1] - 2*c[k+1]
    if kodd == 0
        @inbounds aux = aux - (c[kend+2])^2
    end
    @inbounds c[k+1] = aux / (2*c[1][1])

    return nothing
end
