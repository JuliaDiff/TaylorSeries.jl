# utils_Taylor1.jl: 1-variable Taylor expansions
#
# Last modification: 2014.04.12
#
# Luis Benet & David P. Sanders
# UNAM
#


## Constructors ##
immutable Taylor{T<:Number}
    coeffs :: Array{T,1}
    order :: Int
    ## Inner constructor ##
    function Taylor(coeffs::Array{T,1}, order::Int)
        lencoef = length(coeffs)
        order = max(order, lencoef-1)
        v = zeros(T, order+1)
        @inbounds v[1:lencoef] = coeffs[1:lencoef]
        new(v, order)
    end
end
## Outer constructors ##
Taylor{T<:Number}(x::Taylor{T}, order::Int) = Taylor{T}(x.coeffs, order)
Taylor{T<:Number}(x::Taylor{T}) = Taylor{T}(x.coeffs, x.order)
Taylor{T<:Number}(coeffs::Array{T,1}, order::Int) = Taylor{T}(coeffs, order)
Taylor{T<:Number}(coeffs::Array{T,1}) = Taylor{T}(coeffs, length(coeffs)-1)
Taylor{T<:Number}(x::T, order::Int) = Taylor{T}([x], order)
Taylor{T<:Number}(x::T) = Taylor{T}([x], 0)

## Type, length ##
eltype{T<:Number}(::Taylor{T}) = T
length{T<:Number}(a::Taylor{T}) = a.order

## Conversion and promotion rules ##
convert{T<:Number}(::Type{Taylor{T}}, a::Taylor) = Taylor(convert(Array{T,1}, a.coeffs), a.order)
convert{T<:Number, S<:Number}(::Type{Taylor{T}}, b::Array{S,1}) = Taylor(convert(Array{T,1},b))
convert{T<:Number, S<:Number}(::Type{Taylor{T}}, b::S) = Taylor([convert(T,b)], 0)
promote_rule{T<:Number, S<:Number}(::Type{Taylor{T}}, ::Type{Taylor{S}}) = Taylor{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{Taylor{T}}, ::Type{Array{S,1}}) = Taylor{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{Array{S,1}}, ::Type{Taylor{T}}) = Taylor{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{Taylor{T}}, ::Type{S}) = Taylor{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{S}, ::Type{Taylor{T}}) = Taylor{promote_type(T, S)}

## Auxiliary function ##
function firstnonzero{T<:Number}(a::Taylor{T})
    order = a.order
    nonzero::Int = order+1
    z = zero(T)
    for i = 1:order+1
        if a.coeffs[i] != z
            nonzero = i-1
            break
        end
    end
    nonzero
end
function fixshape{T<:Number, S<:Number}(a::Taylor{T}, b::Taylor{S})
    order = max(a.order, b.order)
    a1, b1 = promote(a, b)
    return Taylor(a1, order), Taylor(b1, order), order
end

## real, imag, conj and ctranspose ##
for f in (:real, :imag, :conj)
    @eval ($f)(a::Taylor) = Taylor(($f)(a.coeffs), a.order)
end
ctranspose(a::Taylor) = conj(a)

## zero and one ##
zero{T<:Number}(a::Taylor{T}) = Taylor(zero(T), a.order)
one{T<:Number}(a::Taylor{T}) = Taylor(one(T), a.order)

## Equality ##
function ==(a::Taylor, b::Taylor)
    a1, b1, order = fixshape(a, b)
    return a1.coeffs == b1.coeffs
end
==(a::Taylor, b::Number) = ==(a, Taylor(b, a.order))
==(a::Number, b::Taylor) = ==(b, Taylor(a, b.order))

## Addition and substraction ##
for f in (:+, :-)
    @eval begin
        function ($f)(a::Taylor, b::Taylor)
            a1, b1, order = fixshape(a, b)
            v = ($f)(a1.coeffs, b1.coeffs)
            return Taylor(v, order)
        end
        ($f)(a::Taylor, b::Number) = ($f)(a, Taylor(b, a.order))
        ($f)(a::Number, b::Taylor) = ($f)(Taylor(a, b.order), b)
        ($f)(a::Taylor) = Taylor(($f)(a.coeffs), a.order)
    end
end

## Multiplication ##
function *(a::Taylor, b::Taylor)
    a1, b1, order = fixshape(a, b)
    T = eltype(a1)
    coeffs = zeros(T,order+1)
    coeffs[1] = a1.coeffs[1] * b1.coeffs[1]
    for k = 1:order
        coeffs[k+1] = mulHomogCoef(k, a1.coeffs, b1.coeffs)
    end
    Taylor(coeffs, order)
end
# Homogeneous coefficient for the multiplication
function mulHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, bc::Array{T,1})
    kcoef == 0 && return ac[1] * bc[1]
    coefhomog = zero(T)
    for i = 0:kcoef
        @inbounds coefhomog += ac[i+1] * bc[kcoef-i+1]
    end
    coefhomog
end
*(a::Taylor, b::Number) = Taylor(b*a.coeffs, a.order)
*(a::Number, b::Taylor) = Taylor(a*b.coeffs, b.order)

## Division ##
function /(a::Taylor, b::Taylor)
    a1, b1, order = fixshape(a, b)
    ordLHopital, cLHopital = divlhopital(a1, b1) # L'Hôpital order and coefficient
    T = typeof(cLHopital)
    v1 = convert(Array{T,1}, a1.coeffs)
    v2 = convert(Array{T,1}, b1.coeffs)
    coeffs = zeros(T,order+1)
    coeffs[1] = cLHopital
    for k = ordLHopital+1:order
        coeffs[k-ordLHopital+1] = divHomogCoef(k, v1, v2, coeffs, ordLHopital)
    end
    Taylor(coeffs, order)
end
function divlhopital(a1::Taylor, b1::Taylor)
    # L'Hôpital order is calculated; a1 and b1 are assumed to be of the same order (length)
    a1nz = firstnonzero(a1)
    b1nz = firstnonzero(b1)
    ordLHopital = min(a1nz, b1nz)
    if ordLHopital > a1.order
        ordLHopital = a1.order
    end
    cLHopital = a1.coeffs[ordLHopital+1] / b1.coeffs[ordLHopital+1]
    aux = abs2(cLHopital)
    # Can L'Hôpital be applied?
    if isinf(aux)
        info("Order k=$(ordLHopital) => coeff[$(ordLHopital+1)]=$(cLHopital)")
        error("Division does not define a Taylor polynomial or its first coefficient is infinite.\n")
    elseif isnan(aux)
        info("Order k=$(ordLHopital) => coeff[$(ordLHopital+1)]=$(cLHopital)")
        error("Impossible to apply L'Hôpital...\n")
    elseif ordLHopital>0
        warn("Applying L'Hôpital. The last k=$(ordLHopital) Taylor coefficients ARE SET to 0.\n")
    end
    return ordLHopital, cLHopital
end
# Homogeneous coefficient for the division
function divHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, bc::Array{T,1}, 
    coeffs::Array{T,1}, ordLHopital::Int)
    #
    kcoef == ordLHopital && return ac[ordLHopital+1] / bc[ordLHopital+1]
    coefhomog = mulHomogCoef(kcoef, coeffs, bc)
    coefhomog = (ac[kcoef+1]-coefhomog) / bc[ordLHopital+1]
    coefhomog
end
/(a::Taylor,b::Number) = Taylor(a.coeffs/b, a.order)
/(a::Number,b::Taylor) = Taylor([a], b.order) / b

## Int power ##
function ^(a::Taylor, n::Integer)
    uno = one(a)
    n < 0 && return uno / a^(-n)
    n == 0 && return uno
    if n%2 == 0     # even power
        n == 2 && return square(a)
        pow = div(n, 2)
        return square( a^pow )
    else            # odd power
        n == 1 && return a
        pow = div(n-1, 2)
        return a*square( a^pow )
    end
end
## Real power ##
function ^(a::Taylor, x::Real)
    uno = one(a)
    if x == zero(x)
        return uno
    elseif x == 0.5
        return sqrt(a)
    end
    order = a.order
    # First non-zero coefficient
    l0nz = firstnonzero(a)
    if l0nz > order
        return zero(a)
    end
    # The first non-zero coefficient of the result; must be integer
    lnull = x*l0nz
    if !isinteger(lnull)
        error("Integer exponent REQUIRED if the Taylor polynomial is expanded around 0.\n")
    end
    # Reaching this point, it is possible to implement the power of the Taylor polynomial. 
    # The last l0nz coefficients are set to zero.
    lnull = itrunc(lnull)
    if l0nz > 0
        warn("The last k=$(l0nz) Taylor coefficients ARE SET to 0.\n")
    end
    aux = (a.coeffs[l0nz+1])^x
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    coeffs = zeros(T, order+1)
    coeffs[lnull+1] = aux
    k0 = lnull+l0nz
    for k = k0+1:order
        coeffs[k-l0nz+1] = powHomogCoef(k, v, x, coeffs, l0nz)
    end
    Taylor(coeffs,order)
end
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
^{T<:Number,S<:Number}(a::Taylor{T}, x::Complex{S}) = exp( x*log(a) )
^(a::Taylor, b::Taylor) = exp( b*log(a) )

## Square ##
function square{T<:Number}(a::Taylor{T})
    order = a.order
    coeffs = zeros(T,order+1)
    coeffs[1] = a.coeffs[1]^2
    for k = 1:order
        coeffs[k+1] = squareHomogCoef(k, a.coeffs)
    end
    Taylor(coeffs,order)
end
# Homogeneous coefficients for square
function squareHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1})
    kcoef == 0 && return ac[1]^2
    coefhomog = zero(T)
    kodd = kcoef%2
    kend = div(kcoef - 2 + kodd, 2)
    for i = 0:kend
        @inbounds coefhomog += ac[i+1]*ac[kcoef-i+1]
    end
    coefhomog = 2coefhomog
    if kodd == 0
        @inbounds coefhomog += ac[div(kcoef,2)+1]^2
    end
    coefhomog
end

## Square root ##
function sqrt(a::Taylor)
    order = a.order
    # First non-zero coefficient
    l0nz = firstnonzero(a)
    if l0nz > order
        return Taylor(zero(T),order)
    elseif l0nz%2 == 1 # l0nz must be pair
        error("First non-vanishing Taylor coefficient must be an EVEN POWER\n",
            "to expand SQRT around 0.\n")
    end
    # Reaching this point, it is possible to implement the sqrt of the Taylor polynomial. 
    # The last l0nz coefficients are set to zero.
    if l0nz > 0
        warn("The last k=$(l0nz) Taylor coefficients ARE SET to 0.\n")
    end
    lnull = div(l0nz, 2)
    aux = sqrt(a.coeffs[l0nz+1])
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    coeffs = zeros(T, order+1)
    coeffs[lnull+1] = aux
    for k = lnull+1:order-l0nz
        coeffs[k+1] = sqrtHomogCoef(k, v, coeffs, lnull)
    end
    Taylor(coeffs, order)
end
# Homogeneous coefficients for the square-root
function sqrtHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, coeffs::Array{T,1}, knull::Int)
    kcoef == knull && return sqrt(ac[2*knull+1])
    coefhomog = zero(T)
    kodd = (kcoef - knull)%2
    kend = div(kcoef - knull - 2 + kodd, 2)
    for i = knull+1:knull+kend
        @inbounds coefhomog += coeffs[i+1]*coeffs[kcoef+knull-i+1]
    end
    @inbounds aux = ac[kcoef+knull+1]-2coefhomog
    if kodd == 0
        @inbounds aux = aux - (coeffs[kend+knull+2])^2
    end
    coefhomog = aux / (2coeffs[knull+1])
    coefhomog
end

## Exp ##
function exp(a::Taylor)
    order = a.order
    aux = exp( a.coeffs[1] )
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    coeffs = zeros(T, order+1)
    coeffs[1] = aux
    for k = 1:order
        coeffs[k+1] = expHomogCoef(k, v, coeffs)
    end
    Taylor( coeffs, order )
end
# Homogeneous coefficients for exp
function expHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, coeffs::Array{T,1})
    kcoef == 0 && return exp(ac[1])
    coefhomog = zero(T)
    for i = 0:kcoef-1
        @inbounds coefhomog += (kcoef-i) * ac[kcoef-i+1] * coeffs[i+1]
    end
    coefhomog = coefhomog/kcoef
    coefhomog
end

## Log ##
function log(a::Taylor)
    order = a.order
    l0nz = firstnonzero(a)
    if firstnonzero(a)>0
        error("Not possible to expand LOG around 0.\n")
    end
    aux = log( a.coeffs[1] )
    T = typeof(aux)
    ac = convert(Array{T,1}, a.coeffs)
    coeffs = zeros(T, order+1)
    coeffs[1] = aux
    for k = 1:order
        coeffs[k+1] = logHomogCoef(k, ac, coeffs)
    end
    Taylor( coeffs, order )
end
# Homogeneous coefficients for log
function logHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, coeffs::Array{T,1})
    kcoef == 0 && return log( ac[1] )
    coefhomog = zero(T)
    for i = 1:kcoef-1
        @inbounds coefhomog += (kcoef-i) * ac[i+1] * coeffs[kcoef-i+1]
    end
    @inbounds coefhomog = (ac[kcoef+1] -coefhomog/kcoef) / ac[1]
    coefhomog
end

## Sin and cos ##
sin(a::Taylor) = sincos(a, "sin")
cos(a::Taylor) = sincos(a, "cos")
function sincos(a::Taylor, fun::String)
    order = a.order
    aux = sin( a.coeffs[1] )
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    sincoeffs = zeros(T,order+1)
    coscoeffs = zeros(T,order+1)
    sincoeffs[1] = aux
    coscoeffs[1] = cos( a.coeffs[1] )
    for k = 1:order
        sincoeffs[k+1], coscoeffs[k+1] = sincosHomogCoef(k, v, sincoeffs, coscoeffs)
    end
    if fun == "sin"
        return Taylor( sincoeffs, order )
    else
        return Taylor( coscoeffs, order )
    end
end
# Homogeneous coefficients for log
function sincosHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, 
        sincoeffs::Array{T,1}, coscoeffs::Array{T,1})
    kcoef == 0 && return sin( ac[1] ), cos( ac[1] )
    sincoefhom = zero(T)
    coscoefhom = zero(T)
    for i = 1:kcoef
        @inbounds begin
            number = i * ac[i+1]
            sincoefhom += number * coscoeffs[kcoef-i+1]
            coscoefhom -= number * sincoeffs[kcoef-i+1]
        end
    end
    sincoefhom = sincoefhom/kcoef
    coscoefhom = coscoefhom/kcoef
    return sincoefhom, coscoefhom
end

## Tan ##
function tan(a::Taylor)
    order = a.order
    aux = tan( a.coeffs[1] )
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    coeffs = zeros(T,order+1)
    coeffst2 = zeros(T,order+1)
    coeffs[1] = aux
    coeffst2[1] = aux^2
    for k = 1:order
        coeffs[k+1] = tanHomogCoef(k, v, coeffst2)
        coeffst2[k+1] = squareHomogCoef(k, coeffs)
    end
    Taylor( coeffs, order )
end
# Homogeneous coefficients for tan
function tanHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, coeffst2::Array{T,1})
    kcoef == 0 && return tan( ac[1] )
    coefhomog = zero(T)
    for i = 0:kcoef-1
        @inbounds coefhomog += (kcoef-i)*ac[kcoef-i+1]*coeffst2[i+1]
    end
    @inbounds coefhomog = ac[kcoef+1] + coefhomog/kcoef
    coefhomog
end

## Differentiating ##
function diffTaylor{T<:Number}(a::Taylor{T})
    order = a.order
    coeffs = zeros(T, order+1)
    coeffs[1] = a.coeffs[2]
    for i = 1:order
        @inbounds coeffs[i] = i*a.coeffs[i+1]
    end
    return Taylor(coeffs, order)
end

## Integrating ##
function integTaylor{T<:Number}(a::Taylor{T}, x::Number)
    order = a.order
    R = promote_type(T, typeof(a.coeffs[1] / 1), typeof(x))
    coeffs = zeros(R, order+1)
    for i = 1:order
        @inbounds coeffs[i+1] = a.coeffs[i] / i
    end
    coeffs[1] = convert(R, x)
    return Taylor(coeffs, order)
end
integTaylor{T<:Number}(a::Taylor{T}) = integTaylor(a, zero(T))

## Evaluates a Taylor polynomial on a given point using Horner's rule ##
function evalTaylor{T}(a::Taylor{T}, dx::Number)
    orden = a.order
    suma = a.coeffs[end]
    for k = orden:-1:1
        @inbounds suma = suma*dx + a.coeffs[k]
    end
    suma
end
evalTaylor{T<:Number}(a::Taylor{T}) = evalTaylor(a, zero(T))

## Returns de n-th derivative of a series expansion
function deriv{T}(a::Taylor{T}, n::Int=1)
    @assert n>= 0
    #n < 0 && error("Needs a non-negative value for the derivative!")
    n > a.order && error(
        "You need to increase the order of the Taylor series (currently ", a.order, ")\n",
        "to calculate its ", n,"-th derivative.")
    res::T = factorial(BigInt(n))*a.coeffs[n+1]
    return res
end

## showcompact ##
function showcompact{T<:Number}(io::IO, a::Taylor{T})
    z = zero(T)
    if a == z
        println(io, a.order, "-order Taylor{", T, "}:\n ", z, " ")
        return
    end
    space = " "
    ser = space
    ifirst = true
    print(io, a.order, "-order Taylor{", T, "}:\n ")
    for i = 0:a.order
        monom = i==0 ? "" : i==1 ? " * x" : string(" * x^", i)
        @inbounds c = a.coeffs[i+1]
        if c != z
            if ifirst
                plusmin = c > 0 ? "" : "-"
                print(io, ser, plusmin, abs(c), monom, space)
                ifirst = false
                continue
            end
            plusmin = c > 0 ? "+ " : "- "
            print(io, ser, plusmin, abs(c), monom, space)
            ser = ""
        end
    end
    print(io, "\n")
end
function showcompact{T<:Complex}(io::IO, a::Taylor{T})
    z = zero(T)
    zre = real(z)
    if a == z
        println(io, a.order, "-order Taylor{", T, "}\n ", zre)
        return
    end
    space = " "
    ifirst = true
    print(io, a.order, "-order Taylor{", T, "}\n ")
    for i = 0:a.order
        monom = i==0 ? "" : i==1 ? " * x" : string(" * x^", i)
        @inbounds c = a.coeffs[i+1]
        if c == z
            continue
        end
        cadena = compactCmplx(c, ifirst)
        print(io, cadena * monom * space )
        ifirst = false
    end
    print(io, ") \n")
end
function compactCmplx{T}(zz::Complex{T}, ifirst::Bool)
    zre = zero(T)
    if zz == zero(Complex{T})
        return zre
    end
    cre, cim = reim(zz)
    if cre > zre
        if ifirst
            cadena = string("( ", abs(cre)," ")
        else
            cadena = string(" + ( ", abs(cre)," ")
        end
        if cim > zre
            cadena = string(cadena, "+ ", abs(cim), " im )")
        elseif cim < zre
            cadena = string(cadena, "- ", abs(cim), " im )")
        else
            cadena = string(cadena, ")")
        end
    elseif cre < zre
        cadena = string(" - ( ", abs(cre), " ")
        if cim > zre
            cadena = string(cadena, "- ", abs(cim), " im )")
        elseif cim < zre
            cadena = string(cadena, "+ ", abs(cim), " im )")
        else
            cadena = string(cadena, ")")
        end
    else
        if cim > zre
            if ifirst
                cadena = string("( ", abs(cim), " im )")
            else
                cadena = string(" + ( ", abs(cim), " im )")
            end
        else
            cadena = string(" - ( ", abs(cim), " im )")
        end
    end
    return cadena
end
