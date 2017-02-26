# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

## Exp ## 
doc"""
    exp(a)

Return the Taylor expansion of $e^a$, of order `a.order`, for
`a::Taylor1` polynomial.

For details on making the Taylor expansion, see
[`TaylorSeries.expHomogCoef`](@ref).
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
doc"""
    expHomogCoef(kcoef, ac, coeffs)

Compute the `k-th` expansion coefficient of $c = \exp(a)$ given by

\begin{equation*}
c_k = \frac{1}{k} \sum_{j=0}^{k-1} (k-j) a_{k-j} c_j,
\end{equation*}

with $a$ a `Taylor1` polynomial.

Inputs are the `kcoef`-th coefficient, the vector of the expansion coefficients
`ac` of $a$, and the already calculated expansion coefficients `coeffs` of `c`.
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
doc"""
    log(a)

Return the Taylor expansion of $\log(a)$, of order `a.order`, for `a::Taylor1`
polynomial.

For details on making the Taylor expansion, see
[`TaylorSeries.logHomogCoef`](@ref).
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
doc"""
    logHomogCoef(kcoef, ac, coeffs)

Compute the `k-th` expansion coefficient of $c = \log(a)$, given by

\begin{equation*}
c_k = \frac{1}{a_0} (a_k - \frac{1}{k} \sum_{j=0}^{k-1} j a_{k-j} c_j ),
\end{equation*}

with $a$ a `Taylor1` polynomial.

Inputs are the `kcoef`-th coefficient, the vector of the expansion coefficients
`ac` of `a`, and the already calculated expansion coefficients `coeffs` of `c`.
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


### TRIGONOMETRIC FUNCTIONS ###

## Sin ## 
doc"""
    sin(a)

Return the Taylor expansion of $\sin(a)$, of order `a.order`, for
`a::Taylor1` polynomial.

For details on making the Taylor expansion, see
[`TaylorSeries.sincosHomogCoef`](@ref).
"""
sin(a::Taylor1) = sincos(a)[1]

## Cos ## 
doc"""
    cos(a)

Return the Taylor expansion of $\cos(a)$, of order `a.order`,
for `a::Taylor1` polynomial

For details on making the Taylor expansion, see
[`TaylorSeries.sincosHomogCoef`](@ref).
"""
cos(a::Taylor1) = sincos(a)[2]

## Sin and Cos ## 
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
doc"""
    sincosHomogCoef(kcoef, ac, scoeffs, ccoeffs)

Compute the `k-th` expansion coefficient of $s = \sin(a)$ and $c=\cos(a)$
simultaneously given by

\begin{eqnarray*}
s_k &=& \frac{1}{k} \sum_{j=0}^{k-1} (k-j) a_{k-j} c_j \\\ \\

c_k &=& -\frac{1}{k}\sum_{j=0}^{k-1} (k-j) a_{k-j} s_j
\end{eqnarray*}

with $a$ a `Taylor1` polynomial.

Inputs are the `kcoef`-th coefficient, the vector of the expansion coefficients
`ac` of `a`, and the already calculated expansion coefficients `scoeffs`
and `ccoeffs` of `sin(a)` and `cos(a)`, respectvely.
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

## tan ##
doc"""
    tan(a)

Return the Taylor expansion of $\tan(a)$, of order `a.order`, for
`a::Taylor1` polynomial.

For details on making the Taylor expansion, see
[`TaylorSeries.tanHomogCoef`](@ref).
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
doc"""
    tanHomogCoef(kcoef, ac, coeffst2)

Compute the `k-th` expansion coefficient of $c = \tan(a)$ given by

\begin{equation*}
c_k = a_k + \frac{1}{k} \sum_{j=0}^{k-1} (k-j) a_{k-j} p_j,
\end{equation*}

with $a$ a `Taylor1` polynomial and $p = c^2$.

Inputs are the `kcoef`-th coefficient, the vector of the expansion coefficients
`ac` of `a`, and the already calculated expansion coefficients `coeffst2`
of `c^2`.
"""
function tanHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, coeffst2::Array{T,1})
    kcoef == 0 && return tan( ac[1] )
    coefhomog = zero(T)
    @inbounds for i = 0:kcoef-1
        coefhomog += (kcoef-i)*ac[kcoef-i+1]*coeffst2[i+1]
    end
    @inbounds coefhomog = ac[kcoef+1] + coefhomog/kcoef
    coefhomog
end

### INVERSE TRIGONOMETRIC FUNCTIONS ### 

## Arcsin ##
doc"""
    asin(a)

Return the Taylor expansion of $\arcsin(a)$, of order `a.order`, for
`a::Taylor1` polynomial.

For details on making the Taylor expansion, see
[`TaylorSeries.asinHomogCoef`](@ref).
"""
function asin(a::Taylor1)
    @inbounds aux = asin( a.coeffs[1] )
    T = typeof(aux)
    ac = convert(Array{T,1}, a.coeffs)
    ac[1]^2 == one(T) && throw(ArgumentError("""
        Recursion formula diverges due to vanishing `sqrt`."""))
    rc = sqrt(one(T) - a^2).coeffs
    coeffs = zeros(ac)
    @inbounds coeffs[1] = aux
    @inbounds for k in 1:a.order
        coeffs[k+1] = asinHomogCoef(k, ac, rc, coeffs)
    end
    Taylor1( coeffs, a.order )
end

# Homogeneous coefficients for arcsin
doc"""
    asinHomogCoef(kcoef, ac, rc, coeffs)

Compute the `k-th` expansion coefficient of $s = \arcsin(a)$ given by

\begin{equation*}
s_k = \frac{1}{ \sqrt{r_0} } \big( a_k - \frac{1}{k}
    \sum_{j=1}^{k-1} j r_{k-j} s_j \big),
\end{equation*}

with $a$ a `Taylor1` polynomial and $r = \sqrt{1 - a^2}$.

Inputs are the `kcoef`-th coefficient, the vector of the expansion coefficients
`ac` of `a`, the already calculated expansion coefficients `rc`
of $r$ (see above), and the already calculated expansion coefficients
`coeffs` of `asin(a)`.
"""
function asinHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, rc::Array{T,1},
        coeffs::Array{T,1})
    kcoef == 0 && return asin( ac[1] )
    coefhomog = zero(T)
    @inbounds for i in 1:kcoef-1
        coefhomog += (kcoef-i) * rc[i+1] * coeffs[kcoef-i+1]
    end
    @inbounds coefhomog = (ac[kcoef+1] - coefhomog/kcoef) / rc[1]
    coefhomog
end

## Arccos ## 
doc"""
    acos(a)

Return the Taylor expansion of $\arccos(a)$, of order `a.order`, for
`a::Taylor1` polynomial.

For details on making the Taylor expansion, see
[`TaylorSeries.acosHomogCoef`](@ref).
"""
function acos(a::Taylor1)
    @inbounds aux = asin( a.coeffs[1] )
    T = typeof(aux)
    ac = convert(Array{T,1}, a.coeffs)
    ac[1]^2 == one(T) && throw(ArgumentError("""
        Recursion formula diverges due to vanishing `sqrt`."""))
    rc = sqrt(one(T) - a^2).coeffs
    coeffs = zeros(ac)
    @inbounds coeffs[1] = aux
    @inbounds for k in 1:a.order
        coeffs[k+1] = acosHomogCoef(k, ac, rc, coeffs)
    end
    @inbounds coeffs[1] = acos( a.coeffs[1] )
    Taylor1( coeffs, a.order )
end

# Homogeneous coefficients for arccos
doc"""
    acosHomogCoef(kcoef, ac, rc, coeffs)

Compute the `k-th` expansion coefficient of $c = \arccos(a)$ given by

\begin{equation*}
c_k = - \frac{1}{ r_0 } \big( a_k - \frac{1}{k}
    \sum_{j=1}^{k-1} j r_{k-j} c_j \big),
\end{equation*}

with $a$ a `Taylor1` polynomial and $r = \sqrt{1 - a^2}$.

Inputs are the `kcoef`-th coefficient, the vector of the expansion coefficients
`ac` of `a`, the already calculated expansion coefficients `rc`
of $r$ (see above), and the already calculated expansion coefficients
`coeffs` of `acos(a)`.
"""
function acosHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, rc::Array{T,1},
        coeffs::Array{T,1})
    kcoef == 0 && return acos( ac[1] )
    asinHomogCoef(kcoef, -ac, rc, coeffs)
end


## Arctan
doc"""
    atan(a)

Return the Taylor expansion of $\arctan(a)$, of order `a.order`, for
`a::Taylor1` polynomial.

For details on making the Taylor expansion, see
[`TaylorSeries.atanHomogCoef`](@ref).
"""
function atan(a::Taylor1)
    @inbounds aux = atan( a.coeffs[1] )
    T = typeof(aux)
    rc = (one(T) + a^2).coeffs
    rc[1] == zero(T) && throw(ArgumentError("""
        Recursion formula has a pole."""))
    ac = convert(Array{T,1}, a.coeffs)
    coeffs = zeros(ac)
    @inbounds coeffs[1] = aux
    @inbounds for k in 1:a.order
        coeffs[k+1] = atanHomogCoef(k , ac, rc, coeffs)
    end
    Taylor1( coeffs, a.order )
end

# Homogeneous coefficients for arctan
doc"""
    atanHomogCoef(kcoef, ac, rc, coeffs)

Compute the `k-th` expansion coefficient of $c = \arctan(a)$ given by

\begin{equation*}
t_k = \frac{1}{r_0}(a_k - \frac{1}{k} \sum_{j=1}^{k-1} j r_{k-j} t_j) ,
\end{equation*}

with $a$ a `Taylor1` polynomial and $r = 1 + a^2$.

Inputs are the `kcoef`-th coefficient, the vector of the expansion coefficients
`ac` of `a`, the already calculated expansion coefficients `rc`
of $r$ (see above), and the already calculated expansion coefficients
`coeffs` of `asin(a)`.
"""
function atanHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, rc::Array{T,1},
        coeffs::Array{T,1})
    kcoef == 0 && return atan( ac[1] )
    coefhomog = zero(T)
    @inbounds for i in 1:kcoef-1
        coefhomog += (kcoef-i) * rc[i+1] * coeffs[kcoef-i+1]
    end
    @inbounds coefhomog = (ac[kcoef+1] - coefhomog/kcoef) / rc[1]
    coefhomog
end

### HYPERBOLIC FUNCTIONS ###

## Sinh ## 
doc"""
    sinh(a)

Return the Taylor expansion of $\sinh(a)$, of order `a.order`, for
`a::Taylor1` polynomial.

For details on making the Taylor expansion, see
[`TaylorSeries.sinhcoshHomogCoef`](@ref).
"""
sinh(a::Taylor1) = sinhcosh(a)[1]

## Cosh ## 
doc"""
    cosh(a)

Return the Taylor expansion of $\cosh(a)$, of order `a.order`,
for `a::Taylor1` polynomial

For details on making the Taylor expansion, see
[`TaylorSeries.sinhcoshHomogCoef`](@ref).
"""
cosh(a::Taylor1) = sinhcosh(a)[2]

## Sinh and Cosh ## 
function sinhcosh(a::Taylor1)
    @inbounds aux = sinh( a.coeffs[1] )
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    sinhcoeffs = similar(v)
    coshcoeffs = similar(v)
    @inbounds sinhcoeffs[1] = aux
    @inbounds coshcoeffs[1] = cosh( a.coeffs[1] )
    @inbounds for k = 1:a.order
        sinhcoeffs[k+1], coshcoeffs[k+1] = sinhcoshHomogCoef(k, v, sinhcoeffs, coshcoeffs)
    end
    return Taylor1(sinhcoeffs, a.order), Taylor1(coshcoeffs, a.order)
end

# Homogeneous coefficients for sinhcosh
doc"""
    sinhcoshHomogCoef(kcoef, ac, scoeffs, ccoeffs)

Compute the `k-th` expansion coefficient of $sh = \sinh(a)$ and $ch=\cosh(a)$
simultaneously given by

\begin{eqnarray*}
sh_k &=& \frac{1}{k} \sum_{j=0}^{k-1} (k-j) a_{k-j} ch_j \\\ \\

ch_k &=& \frac{1}{k}\sum_{j=0}^{k-1} (k-j) a_{k-j} sh_j
\end{eqnarray*}

with $a$ a `Taylor1` polynomial.

Inputs are the `kcoef`-th coefficient, the vector of the expansion coefficients
`ac` of `a`, and the already calculated expansion coefficients `shcoeffs`
and `chcoeffs` of `sinh(a)` and `cosh(a)`, respectvely.
"""
function sinhcoshHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1},
        shcoeffs::Array{T,1}, chcoeffs::Array{T,1})

    kcoef == 0 && return sinh( ac[1] ), cosh( ac[1] )
    sinhcoefhom = zero(T)
    coshcoefhom = zero(T)

    @inbounds for i = 1:kcoef
        x = i * ac[i+1]
        sinhcoefhom += x * chcoeffs[kcoef-i+1]
        coshcoefhom += x * shcoeffs[kcoef-i+1]
    end

    sinhcoefhom = sinhcoefhom/kcoef
    coshcoefhom = coshcoefhom/kcoef
    return sinhcoefhom, coshcoefhom
end


## tanh ##
doc"""
    tanh(a)

Return the Taylor expansion of $\tanh(a)$, of order `a.order`, for
`a::Taylor1` polynomial.

For details on making the Taylor expansion, see
[`TaylorSeries.tanhHomogCoef`](@ref).
"""
function tanh(a::Taylor1)
    aux = tanh( a.coeffs[1] )
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    coeffs = similar(v)
    coeffst2 = similar(v)
    @inbounds coeffs[1] = aux
    @inbounds coeffst2[1] = aux^2
    @inbounds for k = 1:a.order
        coeffs[k+1] = tanhHomogCoef(k, v, coeffst2)
        coeffst2[k+1] = squareHomogCoef(k, coeffs)
    end
    Taylor1( coeffs, a.order )
end

# Homogeneous coefficients for tanh
doc"""
    tanhHomogCoef(kcoef, ac, coeffst2)

Compute the `k-th` expansion coefficient of $c = \tanh(a)$ given by

\begin{equation*}
th_k = a_k - \frac{1}{k} \sum_{j=0}^{k-1} (k-j) a_{k-j} p_j,
\end{equation*}

with $a$ a `Taylor1` polynomial and $p = th^2$.

Inputs are the `kcoef`-th coefficient, the vector of the expansion coefficients
`ac` of `a`, and the already calculated expansion coefficients `coeffst2`
of `th^2`.
"""
function tanhHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, coeffst2::Array{T,1})
    kcoef == 0 && return tanh( ac[1] )
    coefhomog = zero(T)
    @inbounds for i = 0:kcoef-1
        coefhomog += (kcoef-i)*ac[kcoef-i+1]*coeffst2[i+1]
    end
    @inbounds coefhomog = ac[kcoef+1] - coefhomog/kcoef
    coefhomog
end


doc"""
    reverse(f)

Return the Taylor expansion of $f^{-1}(t)$, of order `N = f.order`,
for `f::Taylor1` polynomial if the first coefficient of `f` is zero.
Otherwise, an `ArgumentError` is thrown.

The algorithm implements Lagrange inversion at $t=0$ if $f(0)=0$:
\begin{equation*}
f^{-1}(t) = \sum_{n=1}^{N} \frac{t^n}{n!} \left.
    \frac{{\rm d}^{n-1}}{{\rm d} z^{n-1}}\left(\frac{z}{f(z)}\right)^n
    \right\vert_{z=0}.
\end{equation*}
"""
function reverse{T<:Number}(f::Taylor1{T})
    if f.coeffs[1] != zero(T)
        throw(ArgumentError(
        """Evaluation of Taylor1 series at 0 is non-zero. For high accuracy, revert
        a Taylor1 series with first coefficient 0 and re-expand about f(0)."""))
    end
    z = Taylor1(T,f.order)
    zdivf = z/f
    zdivfpown = zdivf
    S = eltype(zdivf)
    coeffs = zeros(S,f.order+1)

    coeffs[1] = zero(S)
    @inbounds for n = 1:f.order
        coeffs[n+1] = zdivfpown.coeffs[n]/n
        zdivfpown *= zdivf
    end
    Taylor1(coeffs, f.order)
end



## exp ##
function exp(a::TaylorN)
    order = a.order
    @inbounds aux = exp( a.coeffs[1].coeffs[1] )
    T = typeof(aux)
    coeffs = zeros(HomogeneousPolynomial{T}, order)
    @inbounds coeffs[1] = HomogeneousPolynomial(aux, 0)

    @inbounds for ord in eachindex(coeffs)
        ord == order+1 && continue
        @inbounds for j = 0:ord-1
            coeffs[ord+1] += (ord-j) * a.coeffs[ord-j+1] * coeffs[j+1]
        end
        @inbounds coeffs[ord+1] = coeffs[ord+1] / ord
    end

    return TaylorN{T}(coeffs, order)
end

## log ##
function log(a::TaylorN)
    order = a.order
    @inbounds a0 = a.coeffs[1].coeffs[1]
    if a0 == zero(a0)
        throw(ArgumentError(
        """The 0-th order `TaylorN` coefficient must be non-zero
        in order to expand `log` around 0."""))
    end
    l0 = log( a0 )
    T = typeof(l0)
    coeffs = zeros(HomogeneousPolynomial{T}, order)
    @inbounds coeffs[1] = HomogeneousPolynomial(l0)

    @inbounds for ord in eachindex(coeffs)
        ord == order+1 && continue
        @inbounds for j = 1:ord-1
            coeffs[ord+1] += j * a.coeffs[ord-j+1] * coeffs[j+1]
        end
        coeffs[ord+1] = (a.coeffs[ord+1] - coeffs[ord+1] / ord ) / a0
    end

    return TaylorN{T}(coeffs, order)
end

## sin and cos ##
sin(a::TaylorN) = imag( exp(im*a) )
cos(a::TaylorN) = real( exp(im*a) )
# sin(a::TaylorN) = sincos(a)[1]
# cos(a::TaylorN) = sincos(a)[2]
# function sincos(a::TaylorN)
#     order = a.order
#     @inbounds a0 = a.coeffs[1].coeffs[1]
#     s0 = sin( a0 )
#     c0 = cos( a0 )
#     T = typeof(s0)
#     coeffsSin = zeros(HomogeneousPolynomial{T}, order)
#     coeffsCos = zeros(HomogeneousPolynomial{T}, order)
#     @inbounds coeffsSin[1] = HomogeneousPolynomial(s0)
#     @inbounds coeffsCos[1] = HomogeneousPolynomial(c0)
#
#     @inbounds for ord in eachindex(coeffsSin)
#         ord == order+1 && continue
#         @inbounds for j = 0:ord-1
#             coeffsSin[ord+1] += (ord-j) * a.coeffs[ord-j+1] * coeffsCos[j+1]
#             coeffsCos[ord+1] += (ord-j) * a.coeffs[ord-j+1] * coeffsSin[j+1]
#         end
#         @inbounds coeffsSin[ord+1] =  coeffsSin[ord+1] / ord
#         @inbounds coeffsCos[ord+1] = -coeffsCos[ord+1] / ord
#     end
#     return TaylorN{T}(coeffsSin, order), TaylorN{T}(coeffsCos, order)
# end

## tan ##
function tan(a::TaylorN)
    order = a.order
    @inbounds a0 = a.coeffs[1].coeffs[1]
    t0 = tan(a0)
    T = typeof(t0)
    coeffsTan = zeros(HomogeneousPolynomial{T}, order)
    @inbounds coeffsTan[1] = HomogeneousPolynomial(t0)

    @inbounds for ord in eachindex(coeffsTan)
        ord == order+1 && continue
        v = coeffsTan[1:ord]
        tAux = (TaylorN(v, ord))^2
        @inbounds for j = 0:ord-1
            coeffsTan[ord+1] += (ord-j) * a.coeffs[ord-j+1] * tAux.coeffs[j+1]
        end
        coeffsTan[ord+1] = a.coeffs[ord+1] + coeffsTan[ord+1] / ord
    end

    return TaylorN{T}(coeffsTan, order)
end
