# utils_Taylor1.jl: 1-variable Taylor expansions
#
# Last modification: 2015.05.08
#
# Luis Benet & David P. Sanders
# UNAM
#


## Constructors ##
@doc """
    DataType for polynomial expansions in one independent variable

    Fieldnames:

    - `coeffs`: vector containing the expansion coefficients; the i-th
    component is the i-1 coefficient of the expansion

    - `order` : maximum order of the expansion considered
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

## get_coeff
get_coeff(a::Taylor1, n::Int) = (@assert 0 <= n <= a.order+1;
    return a.coeffs[n+1])

## Type, length ##
eltype{T<:Number}(::Taylor1{T}) = T
length{T<:Number}(a::Taylor1{T}) = a.order

## Conversion and promotion rules ##
convert{T<:Number}(::Type{Taylor1{T}}, a::Taylor1) =
    Taylor1(convert(Array{T,1}, a.coeffs), a.order)
convert{T<:Number, S<:Number}(::Type{Taylor1{T}}, b::Array{S,1}) =
    Taylor1(convert(Array{T,1},b))
convert{T<:Number}(::Type{Taylor1{T}}, b::Number) = Taylor1([convert(T,b)], 0)
convert{T<:Number}(::Type{Taylor1{T}}, a::Taylor1{T}) = a
convert{T<:Number}(::Type{Taylor1{T}}, b::Array{T,1}) = Taylor1(b)
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
        ($f)(a::Taylor1) = Taylor1(($f)(a.coeffs), a.order)
    end
end

## Multiplication ##
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
function mulHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, bc::Array{T,1})
    kcoef == 0 && return ac[1] * bc[1]
    coefhomog = zero(T)
    @inbounds for i = 0:kcoef
        coefhomog += ac[i+1] * bc[kcoef-i+1]
    end
    coefhomog
end

## Division ##
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
    aux = abs2(cdivfact)

    # Is the polynomial factorizable?
    if isinf(aux) || isnan(aux)
        # info("Order k=$(orddivfact) => coeff[$(orddivfact+1)]=$(cdivfact)")
        error("Division does not define a Taylor1 polynomial\n",
            " or its first non-zero coefficient is Inf/NaN.\n",
            "Order k=$(orddivfact) => coeff[$(orddivfact+1)]=$(cdivfact).")
    ##else orddivfact>0
    ##    warn("Factorizing the polynomial.\n",
    ##        "The last k=$(orddivfact) Taylor1 coefficients ARE SET to 0.\n")
    end

    return orddivfact, cdivfact
end

# Homogeneous coefficient for the division
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
        function ($op){T<:Real}(a::Taylor1{T}, x::Real)
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

## Int power ##
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
^(a::Taylor1,x::Rational) = a^(x.num/x.den)

^(a::Taylor1, b::Taylor1) = exp( b*log(a) )

## Real power ##
function ^(a::Taylor1, x::Real)
    x == zero(x) && return one(a)
    x == one(x)/2 && return sqrt(a)
    isinteger(x) && return a^round(Int,x)

    # First non-zero coefficient
    l0nz = firstnonzero(a)
    l0nz > a.order && return zero(a)

    # The first non-zero coefficient of the result; must be integer
    lnull = x*l0nz
    !isinteger(lnull) &&
        error("""The 0th order Taylor1 coefficient must be non-zero
        to raise the Taylor1 polynomial to a non-integer exponent""")

    # Reaching this point, it is possible to implement the power of the Taylor1
    # polynomial. The last l0nz coefficients are set to zero.
    lnull = trunc(Int,lnull)
    #l0nz > 0 && warn("The last k=$(l0nz) Taylor1 coefficients ARE SET to 0.\n")
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
^(a::Taylor1, x::Complex) = exp( x*log(a) )

# Homogeneous coefficients for real power
function powHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, x::Real,
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
    l0nz = firstnonzero(a)
    if l0nz > a.order
        return zero(a)
    elseif l0nz%2 == 1 # l0nz must be pair
        error("""First non-vanishing Taylor1 coefficient must correspond
        to an **even power** in order to expand `sqrt` around 0""")
    end

    # Reaching this point, it is possible to implement the sqrt of the Taylor1 polynomial.
    # The last l0nz coefficients are set to zero.
    ##l0nz > 0 && warn("The last k=$(l0nz) Taylor1 coefficients ARE SET to 0.\n")
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
function log(a::Taylor1)
    ( firstnonzero(a)>0 ) && error("Impossible to expand `log` around 0.")
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
sin(a::Taylor1) = sincos(a)[1]
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
function diffTaylor(a::Taylor1)
    coeffs = zero(a.coeffs)
    @inbounds coeffs[1] = a.coeffs[2]
    @inbounds for i = 1:a.order
        coeffs[i] = i*a.coeffs[i+1]
    end
    return Taylor1(coeffs, a.order)
end

## Integrating ##
function integTaylor{T<:Number}(a::Taylor1{T}, x::Number)
    R = promote_type(T, typeof(a.coeffs[1] / 1), typeof(x))
    coeffs = zeros(R, a.order+1)
    @inbounds for i = 1:a.order
        coeffs[i+1] = a.coeffs[i] / i
    end
    @inbounds coeffs[1] = convert(R, x)
    return Taylor1(coeffs, a.order)
end
integTaylor{T<:Number}(a::Taylor1{T}) = integTaylor(a, zero(T))

## Evaluates a Taylor1 polynomial on a given point using Horner's rule ##
function evalTaylor{T<:Number,S<:Number}(a::Taylor1{T}, dx::S)
    R = promote_type(T,S)
    @inbounds suma = convert(R, a.coeffs[end])
    @inbounds for k = a.order:-1:1
        suma = suma*dx + a.coeffs[k]
    end
    suma
end
evalTaylor{T<:Number}(a::Taylor1{T}) = evalTaylor(a, zero(T))

function evalTaylor{T<:Number,S<:Number}(a::Taylor1{T}, x::Taylor1{S})
    a, x = fixshape(a, x)
    @inbounds suma = a.coeffs[end]
    @inbounds for k = a.order:-1:1
        suma = suma*x + a.coeffs[k]
    end
    suma
end

## Returns de n-th derivative of a series expansion
function deriv{T}(a::Taylor1{T}, n::Int=1)
    @assert a.order >= n >= 0
    res::T = factorial( widen(n) ) * a.coeffs[n+1]
    return res
end
