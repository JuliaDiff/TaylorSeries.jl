# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#


## Constructors ##
@doc """
    Taylor1{T<:Number} <: Number

DataType for polynomial expansions in one independent variable.

**Fields:**

- `coeffs :: Array{T,1}` Expansion coefficients; the \$i\$-th
component is the coefficient of degree \$i-1\$ of the expansion.
- `order  :: Int64` Maximum order (degree) of the polynomial.
""" ->
immutable Taylor1{T<:Number} <: Number
    coeffs :: Array{T,1}
    order :: Int

    ## Inner constructor ##
    function Taylor1(coeffs::Array{T,1}, order::Int)
        lencoef = length(coeffs)
        order = max(order, lencoef-1)
        order == lencoef-1 && return new(coeffs, order)
        resize!(coeffs, order+1)
        for i = lencoef+1:order+1
            coeffs[i] = zero(T)
        end
        new(coeffs, order)
    end
end

## Outer constructors ##
Taylor1{T<:Number}(x::Taylor1{T}, order::Int) = Taylor1{T}(x.coeffs, order)
Taylor1{T<:Number}(x::Taylor1{T}) = x
Taylor1{T<:Number}(coeffs::Array{T,1}, order::Int) = Taylor1{T}(coeffs, order)
Taylor1{T<:Number}(coeffs::Array{T,1}) = Taylor1{T}(coeffs, length(coeffs)-1)
Taylor1{T<:Number}(x::T, order::Int) = Taylor1{T}([x], order)
Taylor1{T<:Number}(x::T) = Taylor1{T}([x], 0)

# Shortcut to define Taylor1 independent variables
taylor1_variable(T::Type, order::Int=1) = Taylor1{T}( [zero(T), one(T)], order)
taylor1_variable(order::Int=1) = taylor1_variable(Float64, order)
@doc """
    taylor1_variable(T, [order=1])
    taylor1_variable([order=1])

Short-cut to define the independent variable as a `Taylor1` polynomial of
given `order`. If `T::Type` is ommitted, `Float64` is assumend.
""" taylor1_variable

## get_coeff
@doc """
    get_coeff(a, n)

Return the coefficient of order `n::Int` of a `a::Taylor1` polynomial.
""" ->
get_coeff(a::Taylor1, n::Int) = (@assert 0 <= n <= a.order+1;
    return a.coeffs[n+1])

## Type, length ##
eltype{T<:Number}(::Taylor1{T}) = T
length{T<:Number}(a::Taylor1{T}) = a.order

## Conversion and promotion rules ##
convert{T<:Number}(::Type{Taylor1{T}}, a::Taylor1) =
    Taylor1(convert(Array{T,1}, a.coeffs), a.order)
function convert{T<:Integer, S<:AbstractFloat}(::Type{Taylor1{Rational{T}}},
    a::Taylor1{S})
    v = Array(Rational{T}, length(a.coeffs))
    for i in eachindex(v)
        # v[i] = convert(Rational{T}, a.coeffs[i])
        v[i] = rationalize(a.coeffs[i], tol=eps(one(S)))
    end
    Taylor1(v)
end
convert{T<:Number, S<:Number}(::Type{Taylor1{T}}, b::Array{S,1}) =
    Taylor1(convert(Array{T,1},b))
convert{T<:Number, S<:Number}(::Type{Taylor1{T}}, b::S) = Taylor1([convert(T,b)], 0)
convert{T<:Number}(::Type{Taylor1{T}}, b::T) = Taylor1([b], 0)

promote_rule{T<:Number, S<:Number}(::Type{Taylor1{T}}, ::Type{Taylor1{S}}) =
    Taylor1{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{Taylor1{T}}, ::Type{Array{S,1}}) =
    Taylor1{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{Taylor1{T}}, ::Type{S}) =
    Taylor1{promote_type(T, S)}

## Auxiliary function ##
function firstnonzero{T<:Number}(a::Taylor1{T})
    nonzero::Int = a.order+1
    for i in eachindex(a.coeffs)
        if a.coeffs[i] != zero(T)
            nonzero = i-1
            break
        end
    end
    nonzero
end

function fixshape{T<:Number, S<:Number}(a::Taylor1{T}, b::Taylor1{S})
    eltype(a) == eltype(b) && a.order == b.order && return a, b
    if eltype(a) != eltype(b)
        a, b = promote(a, b)
    end
    if a.order == b.order
        return a, b
    elseif a.order < b.order
        return Taylor1(a, b.order), b
    end
    return a, Taylor1(b, a.order)
end

## real, imag, conj and ctranspose ##
for f in (:real, :imag, :conj)
    @eval ($f){T<:Number}(a::Taylor1{T}) = Taylor1(($f)(a.coeffs), a.order)
end
ctranspose{T<:Number}(a::Taylor1{T}) = conj(a)

## zero and one ##
zero{T<:Number}(a::Taylor1{T}) = Taylor1(zero(T), a.order)
one{T<:Number}(a::Taylor1{T}) = Taylor1(one(T), a.order)

## Equality ##
function ==(a::Taylor1, b::Taylor1)
    a, b = fixshape(a, b)
    return a.coeffs == b.coeffs
end

# Tests `isinf` and `isnan` for all polynomial coefficients
for f in (:isinf, :isnan)
    @eval begin
        function ($f)(a::Taylor1)
            test = false
            for i in eachindex(a.coeffs)
                @inbounds test = ($f)(a.coeffs[i])
                test && break
            end
            test
        end
    end
end

## Addition and substraction ##
for f in (:+, :-)
    @eval begin
        function ($f)(a::Taylor1, b::Taylor1)
            a, b = fixshape(a, b)
            v = similar(a.coeffs)
            @simd for i in eachindex(a.coeffs)
                @inbounds v[i] = ($f)(a.coeffs[i], b.coeffs[i])
            end
            return Taylor1(v, a.order)
        end
        function ($f)(a::Taylor1)
            v = similar(a.coeffs)
            @simd for i in eachindex(v)
                @inbounds v[i] = ($f)(a.coeffs[i])
            end
            return Taylor1(v, a.order)
        end
        function ($f)(a::Taylor1, b::Union{Real,Complex})
            @inbounds aux = ($f)(a.coeffs[1], b)
            v = Array(typeof(aux), length(a.coeffs))
            @simd for i in eachindex(v)
                @inbounds v[i] = a.coeffs[i]
            end
            @inbounds v[1] = aux
            Taylor1(v, a.order)
        end
        function ($f)(a::Union{Real,Complex}, b::Taylor1)
            @inbounds aux = ($f)(a, b.coeffs[1])
            v = Array(typeof(aux), length(b.coeffs))
            @simd for i in eachindex(v)
                @inbounds v[i] = ($f)(b.coeffs[i])
            end
            @inbounds v[1] = aux
            Taylor1(v, b.order)
        end
    end
end

## Multiplication ##
*(a::Bool, b::Taylor1) = *(promote(a,b)...)
*(a::Taylor1, b::Bool) = b*a
function *(a::Union{Real,Complex}, b::Taylor1)
    @inbounds aux = a * b.coeffs[1]
    v = Array(typeof(aux), length(b.coeffs))
    @simd for i in eachindex(v)
        @inbounds v[i] = a * b.coeffs[i]
    end
    Taylor1(v, b.order)
end
*(a::Taylor1, b::Union{Real,Complex}) = b * a
function *(a::Taylor1, b::Taylor1)
    a, b = fixshape(a, b)
    coeffs = similar(a.coeffs)
    @inbounds coeffs[1] = a.coeffs[1] * b.coeffs[1]
    @inbounds for k = 1:a.order
        coeffs[k+1] = mulHomogCoef(k, a.coeffs, b.coeffs)
    end
    Taylor1(coeffs, a.order)
end

# Homogeneous coefficient for the multiplication
"""
Let be ``f(x)``  and ``g(x)`` analitical functions, then, the **`r-th` expansion coefficient** of ``p(x) = f(x) g(x)`` is

<center>
``pᵣ = ∑ᵏ fⱼ gᵣ₋ⱼ``.
"""
function mulHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, bc::Array{T,1})
    kcoef == 0 && return ac[1] * bc[1]
    coefhomog = zero(T)
    @inbounds for i = 0:kcoef
        coefhomog += ac[i+1] * bc[kcoef-i+1]
    end
    coefhomog
end

## Division ##
/{T<:Real}(a::Taylor1, b::T) = a * inv(b)
/{T<:Complex}(a::Taylor1, b::T) = a * inv(b)

function /(a::Taylor1, b::Taylor1)
    a, b = fixshape(a, b)
    # order and coefficient of first factorized term
    orddivfact, cdivfact = divfactorization(a, b)
    T = typeof(cdivfact)
    v1 = convert(Array{T,1}, a.coeffs)
    v2 = convert(Array{T,1}, b.coeffs)
    coeffs = zeros(T, a.order+1)
    @inbounds coeffs[1] = cdivfact
    @inbounds for k = orddivfact+1:a.order
        coeffs[k-orddivfact+1] = divHomogCoef(k, v1, v2, coeffs, orddivfact)
    end
    Taylor1(coeffs, a.order)
end

function divfactorization(a1::Taylor1, b1::Taylor1)
    # order of first factorized term; a1 and b1 assumed to be of the same order
    a1nz = firstnonzero(a1)
    b1nz = firstnonzero(b1)
    orddivfact = min(a1nz, b1nz)
    if orddivfact > a1.order
        orddivfact = a1.order
    end
    cdivfact = a1.coeffs[orddivfact+1] / b1.coeffs[orddivfact+1]

    # Is the polynomial factorizable?
    if isinf(cdivfact) || isnan(cdivfact)
        throw(ArgumentError(
        """Division does not define a Taylor1 polynomial
        or its first non-zero coefficient is Inf/NaN.
        Order k=$(orddivfact) => coeff[$(orddivfact+1)]=$(cdivfact)."""))
    end

    return orddivfact, cdivfact
end

# Homogeneous coefficient for the division
"""
Let be ``f(x)``  and ``g(x)`` analitical functions, then, the **`r-th` expansion coefficient** of ``d(x) = f(x) / g(x)`` is

<center>
``dᵣ = 1/g₀ (fᵣ  -  ∑ᵏ⁻ⁱ dⱼ gᵣ₋ⱼ)``.
"""
function divHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, bc::Array{T,1},
    coeffs::Array{T,1}, ordfact::Int)
    #
    @inbounds kcoef == ordfact && return ac[ordfact+1] / bc[ordfact+1]
    coefhomog = mulHomogCoef(kcoef, coeffs, bc)
    @inbounds coefhomog = (ac[kcoef+1]-coefhomog) / bc[ordfact+1]
    coefhomog
end

## Division functions: rem and mod
for op in (:mod, :rem)
    @eval begin
        function ($op){T<:Real, S<:Real}(a::Taylor1{T}, x::S)
            coeffs = copy(a.coeffs)
            @inbounds coeffs[1] = ($op)(a.coeffs[1], x)
            return Taylor1(coeffs, a.order)
        end
    end
end

function mod2pi{T<:Real}(a::Taylor1{T})
    coeffs = copy(a.coeffs)
    @inbounds coeffs[1] = mod2pi( a.coeffs[1] )
    return Taylor1( coeffs, a.order)
end

## abs function ##
"""
    abs(a)

For `a::Taylor1`, return `a` or `-a` depending on the 0-th order coefficient
of `a`. If it is zero, it throws an `ArgumentError`.
"""
function abs{T<:Real}(a::Taylor1{T})
    if a.coeffs[1] > zero(T)
        return a
    elseif a.coeffs[1] < zero(T)
        return -a
    else
        throw(ArgumentError(
        """The 0th order Taylor1 coefficient must be non-zero
        (abs(x) is not differentiable at x=0)."""))
    end
end

## Int power ##
"""
    ^(a, x)

For `a::Taylor1` and `x::Number`, return `a^x` as an `Taylor1` object.
"""
function ^{T<:Number}(a::Taylor1{T}, n::Integer)
    n == 0 && return one(a)
    n == 1 && return a
    n == 2 && return square(a)
    n < 0 && return inv( a^(-n) )
    return power_by_squaring(a, n)
end

function ^{T<:Integer}(a::Taylor1{T}, n::Integer)
    n == 0 && return one(a)
    n == 1 && return a
    n == 2 && return square(a)
    n < 0 && throw(DomainError())
    return power_by_squaring(a, n)
end

## power_by_squaring; slightly modified from base/intfuncs.jl
## Licensed under MIT "Expat"
function power_by_squaring(x::Taylor1, p::Integer)
    p == 1 && return copy(x)
    p == 0 && return one(x)
    p == 2 && return square(x)
    t = trailing_zeros(p) + 1
    p >>= t

    while (t -= 1) > 0
        x = square(x)
    end

    y = x
    while p > 0
        t = trailing_zeros(p) + 1
        p >>= t
        while (t -= 1) >= 0
            x = square(x)
        end
        y *= x
    end

    return y
end

## Rational power ##
^{T<:Integer}(a::Taylor1,x::Rational{T}) = a^(x.num/x.den)

^(a::Taylor1, b::Taylor1) = exp( b*log(a) )

## Real power ##
"""
If `x::Real` and the 0th order coefficient is non-zero, an `ArgumentError` is thrown.
"""
function ^{S<:Real}(a::Taylor1, x::S)
    x == zero(x) && return one(a)
    x == one(x)/2 && return sqrt(a)
    isinteger(x) && return a^round(Int,x)

    # First non-zero coefficient
    l0nz = firstnonzero(a)
    l0nz > a.order && return zero(a)

    # The first non-zero coefficient of the result; must be integer
    lnull = x*l0nz
    !isinteger(lnull) &&
        throw(ArgumentError(
        """The 0th order Taylor1 coefficient must be non-zero
        to raise the Taylor1 polynomial to a non-integer exponent."""))

    # Reaching this point, it is possible to implement the power of the Taylor1
    # polynomial. The last l0nz coefficients are set to zero.
    lnull = trunc(Int,lnull)
    #l0nz > 0 && warn("The last k=$(l0nz) Taylor1 coefficients ARE SET to 0.")
    @inbounds aux = (a.coeffs[l0nz+1])^x
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    coeffs = zeros(T, a.order+1)
    @inbounds coeffs[lnull+1] = aux
    k0 = lnull+l0nz
    @inbounds for k = k0+1:a.order
        coeffs[k-l0nz+1] = powHomogCoef(k, v, x, coeffs, l0nz)
    end

    Taylor1(coeffs,a.order)
end
^{T<:Complex}(a::Taylor1, x::T) = exp( x*log(a) )

# Homogeneous coefficients for real power
"""
Let be ``f(x)`` an analitical function, then, the **`r-th` expansion coefficient** of ``p(x) = f(x)ᵝ`` is

<center>
``pᵣ = 1/(kf₀) ∑ᵏ⁻¹(β(k-j) -j)fᵣ₋ⱼ pⱼ)``.
"""
function powHomogCoef{T<:Number, S<:Real}(kcoef::Int, ac::Array{T,1}, x::S,
    coeffs::Array{T,1}, knull::Int)

    kcoef == knull && return (ac[knull+1])^x
    coefhomog = zero(T)
    for i = 0:kcoef-knull-1
        aux = x*(kcoef-i)-i
        @inbounds coefhomog += aux*ac[kcoef-i+1]*coeffs[i+1]
    end
    aux = kcoef - knull*(x+1)
    @inbounds coefhomog = coefhomog / (aux*ac[knull+1])

    coefhomog
end

## Square ##
function square(a::Taylor1)
    coeffs = Array(eltype(a), a.order+1)
    coeffs[1] = a.coeffs[1]^2
    @inbounds for k = 1:a.order
        coeffs[k+1] = squareHomogCoef(k, a.coeffs)
    end
    Taylor1(coeffs,a.order)
end

# Homogeneous coefficients for square
"""
Let be ``f(x)`` an analitical function, then, the **`r-th` expansion coefficient** of ``s(x) = f(x)²`` is

<center>
``sᵣ = 2 ∑ᶹfᵣ₋ⱼ fⱼ``
</center>

when `r`is **odd**, with `υ = (r-1)/2`

<center>
``sᵣ = 2 ∑ᵠ(fᵣ₋ⱼ fⱼ + (f_₍ᵣ/₂₎)²``
</center>

when `r`is **even**, with `φ = (r-2)/2`.

"""
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
"""
    sqrt(a)

For `a::Taylor1`, returns the square root of `a` expansion of order `a.order` as an `Taylor1` object.

If the first non-vanishing coefficient of `a` is an **odd power**, and `ArgumentError` will be thrown.  
"""
function sqrt(a::Taylor1)
    # First non-zero coefficient
    l0nz = firstnonzero(a)
    if l0nz > a.order
        return zero(a)
    elseif l0nz%2 == 1 # l0nz must be pair
        throw(ArgumentError(
        """First non-vanishing Taylor1 coefficient must correspond
        to an **even power** in order to expand `sqrt` around 0."""))
    end

    # Reaching this point, it is possible to implement the sqrt of the Taylor1 polynomial.
    # The last l0nz coefficients are set to zero.
    ##l0nz > 0 && warn("The last k=$(l0nz) Taylor1 coefficients ARE SET to 0.")
    lnull = div(l0nz, 2)
    @inbounds aux = sqrt(a.coeffs[l0nz+1])
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    coeffs = zeros(T, a.order+1)
    @inbounds coeffs[lnull+1] = aux
    @inbounds for k = lnull+1:a.order-l0nz
        coeffs[k+1] = sqrtHomogCoef(k, v, coeffs, lnull)
    end
    Taylor1(coeffs, a.order)
end

# Homogeneous coefficients for the square-root
"""
Let be ``f(x)`` an analitical function, then, the **`r-th` expansion coefficient** of ``s(x) = √f(x)`` is

<center>
``sᵣ = 1/(2 s₀) ( fᵣ - 2 ∑ᶹsᵣ₋ⱼ sⱼ )``
</center>

when `r`is **odd**, with `υ = (r-1)/2`

<center>
``sᵣ = 1/(2 s₀) ( fᵣ - 2 ∑ᵠ(sᵣ₋ⱼ sⱼ - (s_₍ᵣ/₂₎)² )``
</center>

when `r`is **even**, with `φ = (r-2)/2`.

"""
function sqrtHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, coeffs::Array{T,1}, knull::Int)
    kcoef == knull && return sqrt(ac[2*knull+1])
    coefhomog = zero(T)
    kodd = (kcoef - knull)%2
    kend = div(kcoef - knull - 2 + kodd, 2)
    two = convert(T,2)
    @inbounds for i = knull+1:knull+kend
        coefhomog += coeffs[i+1]*coeffs[kcoef+knull-i+1]
    end
    @inbounds aux = ac[kcoef+knull+1] - two*coefhomog

    if kodd == 0
        @inbounds aux = aux - (coeffs[kend+knull+2])^2
    end
    @inbounds coefhomog = aux / (two*coeffs[knull+1])

    coefhomog
end

## Exp ##
"""
    exp(a)

For `a::Taylor1`, computes `e^a` of order `a.order` as an `Taylor1` object.
"""
function exp(a::Taylor1)
    @inbounds aux = exp( a.coeffs[1] )
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    coeffs = similar(v)
    @inbounds coeffs[1] = aux
    @inbounds for k = 1:a.order
        coeffs[k+1] = expHomogCoef(k, v, coeffs)
    end
    Taylor1( coeffs, a.order )
end

# Homogeneous coefficients for exp
"""
Let be ``f(x)`` an analitical function, then, the **`r-th` expansion coefficient** of ``e(x) = exp(f(x))`` is

<center>
``pᵣ = 1/k ∑ᵏ⁻¹(k-j)fᵣ₋ⱼ eⱼ``.
"""
function expHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, coeffs::Array{T,1})
    kcoef == 0 && return exp(ac[1])
    coefhomog = zero(T)
    @inbounds for i = 0:kcoef-1
        coefhomog += (kcoef-i) * ac[kcoef-i+1] * coeffs[i+1]
    end
    coefhomog = coefhomog/kcoef
    coefhomog
end

## Log ##
"""
    log(a)

For `a::Taylor1`, computes the natural logarithm's expansion of `a` of order `a.order` as an `Taylor1` object.
"""
function log(a::Taylor1)
    ( firstnonzero(a)>0 ) && throw(
        ArgumentError("Impossible to expand `log` around 0."))
    @inbounds aux = log( a.coeffs[1] )
    T = typeof(aux)
    ac = convert(Array{T,1}, a.coeffs)
    coeffs = similar(ac)
    @inbounds coeffs[1] = aux
    @inbounds for k = 1:a.order
        coeffs[k+1] = logHomogCoef(k, ac, coeffs)
    end
    Taylor1( coeffs, a.order )
end

# Homogeneous coefficients for log
"""
Let be ``f(x)`` an analitical function, then, the **`r-th` expansion coefficient** of ``l(x) = log(f(x))`` is

<center>
``lᵣ = 1/f₀ ( fᵣ - 1/k ∑ᵏ⁻¹ j fᵣ₋ⱼ lⱼ )``.
"""
function logHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, coeffs::Array{T,1})
    kcoef == 0 && return log( ac[1] )
    coefhomog = zero(T)
    @inbounds for i = 1:kcoef-1
        coefhomog += (kcoef-i) * ac[i+1] * coeffs[kcoef-i+1]
    end
    @inbounds coefhomog = (ac[kcoef+1] -coefhomog/kcoef) / ac[1]
    coefhomog
end

## Sin and cos ##
"""
    sin(a)

For `a::Taylor1`, computes the sine's expansion of `a` of order `a.order` as an `Taylor1` object.
"""
sin(a::Taylor1) = sincos(a)[1]
"""
    cos(a)

For `a::Taylor1`, computes the cosine's expansion of `a` of order `a.order` as an `Taylor1` object.
"""
cos(a::Taylor1) = sincos(a)[2]
function sincos(a::Taylor1)
    @inbounds aux = sin( a.coeffs[1] )
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    sincoeffs = similar(v)
    coscoeffs = similar(v)
    @inbounds sincoeffs[1] = aux
    @inbounds coscoeffs[1] = cos( a.coeffs[1] )
    @inbounds for k = 1:a.order
        sincoeffs[k+1], coscoeffs[k+1] = sincosHomogCoef(k, v, sincoeffs, coscoeffs)
    end
    return Taylor1(sincoeffs, a.order), Taylor1(coscoeffs, a.order)
end

# Homogeneous coefficients for sincos
"""
Let be ``f(x)`` an analitical function, then, the **`r-th` expansion coefficients ** of ``s(x) = sin(f(x))`` and ``c(x) = cos(f(x))`` are

<center>
``sᵣ = 1/k ∑ᵏ⁻¹(r-j)fᵣ₋ⱼ cⱼ ``.
</center>

and

<center>
``cᵣ = -1/k ∑ᵏ⁻¹(r-j)fᵣ₋ⱼ sⱼ ``.
</center>
"""
function sincosHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1},
    scoeffs::Array{T,1}, ccoeffs::Array{T,1})

    kcoef == 0 && return sin( ac[1] ), cos( ac[1] )
    sincoefhom = zero(T)
    coscoefhom = zero(T)

    @inbounds for i = 1:kcoef
        x = i * ac[i+1]
        sincoefhom += x * ccoeffs[kcoef-i+1]
        coscoefhom -= x * scoeffs[kcoef-i+1]
    end

    sincoefhom = sincoefhom/kcoef
    coscoefhom = coscoefhom/kcoef
    return sincoefhom, coscoefhom
end

## Tan ##
"""
    tan(a)

For `a::Taylor1`, computes the tangent's expansion of `a` of order `a.order` as an `Taylor1` object.
"""
function tan(a::Taylor1)
    aux = tan( a.coeffs[1] )
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    coeffs = similar(v)
    coeffst2 = similar(v)
    @inbounds coeffs[1] = aux
    @inbounds coeffst2[1] = aux^2
    @inbounds for k = 1:a.order
        coeffs[k+1] = tanHomogCoef(k, v, coeffst2)
        coeffst2[k+1] = squareHomogCoef(k, coeffs)
    end
    Taylor1( coeffs, a.order )
end
# Homogeneous coefficients for tan
"""
Let be ``f(x)`` an analitical function, then, the **`r-th` expansion coefficient** of ``t(x) = tan(f(x))`` with ``p(x) = t(x)²`` is

<center>
``tᵣ = fᵣ + 1/k ∑ᵏ⁻¹(r-j)fᵣ₋ⱼ pⱼ )``.
"""
function tanHomogCoef{T<:Number}(kcoef::Int,ac::Array{T,1},coeffst2::Array{T,1})
    kcoef == 0 && return tan( ac[1] )
    coefhomog = zero(T)
    @inbounds for i = 0:kcoef-1
        coefhomog += (kcoef-i)*ac[kcoef-i+1]*coeffst2[i+1]
    end
    @inbounds coefhomog = ac[kcoef+1] + coefhomog/kcoef
    coefhomog
end

## Differentiating ##
"""
    diffTaylor(a::Taylor1)

Return the `Taylor1` polynomial of the differential of `a::Taylor1`; the last
coefficient is set to zero.
"""
function diffTaylor(a::Taylor1)
    coeffs = zero(a.coeffs)
    @inbounds coeffs[1] = a.coeffs[2]
    @inbounds for i = 1:a.order
        coeffs[i] = i*a.coeffs[i+1]
    end
    return Taylor1(coeffs, a.order)
end

## Integrating ##
"""
    integTaylor(a, x)
    integTaylor(a)

Return the integral of `a::Taylor1`. The constant of integration
(0th order coefficient) is set to `x`, which is zero if ommitted.
"""
function integTaylor{T<:Number, S<:Number}(a::Taylor1{T}, x::S)
    R = promote_type(T, typeof(a.coeffs[1] / 1), S)
    coeffs = zeros(R, a.order+1)
    @inbounds for i = 1:a.order
        coeffs[i+1] = a.coeffs[i] / i
    end
    @inbounds coeffs[1] = convert(R, x)
    return Taylor1(coeffs, a.order)
end
integTaylor{T<:Number}(a::Taylor1{T}) = integTaylor(a, zero(T))

## Evaluates a Taylor1 polynomial on a given point using Horner's rule ##
"""
    evaluate(a, dx)
    evaluate(a)

Evaluate a `Taylor1` polynomial using Horner's rule (hand coded).
"""
function evaluate{T<:Number,S<:Number}(a::Taylor1{T}, dx::S)
    R = promote_type(T,S)
    @inbounds suma = convert(R, a.coeffs[end])
    @inbounds for k = a.order:-1:1
        suma = suma*dx + a.coeffs[k]
    end
    suma
end
evaluate{T<:Number}(a::Taylor1{T}) = a.coeffs[1]

"""
    evaluate(a, x)

Return the substitution of `x::Taylor1` as independent variable in
`a::Taylor1`.
"""
function evaluate{T<:Number,S<:Number}(a::Taylor1{T}, x::Taylor1{S})
    a, x = fixshape(a, x)
    @inbounds suma = a.coeffs[end]
    @inbounds for k = a.order:-1:1
        suma = suma*x + a.coeffs[k]
    end
    suma
end

## Returns de n-th derivative of a series expansion
"""
    deriv(a, [n=1])

Return the value of the `n`-th derivative of `a`.
"""
function deriv{T<:Number}(a::Taylor1{T}, n::Int=1)
    @assert a.order >= n >= 0
    factorial( widen(n) ) * a.coeffs[n+1] :: T
end
