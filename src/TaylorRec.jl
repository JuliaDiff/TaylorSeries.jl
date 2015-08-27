# utils_TaylorRecRec.jl: TaylorRec polynomials in recursive dense representation
#
# Last modification: 2015.04.21
#
# Luis Benet & David P. Sanders
# UNAM
#

abstract AbstractSeries{T<:Number,N}# <: Number
export AbstractSeries

export TaylorRec, taylorRec_variable, @taylorRec


@doc doc"""Basic (immutable) type structure for polynomial expansions in N
independent variables in the recursive dense representation

`coeffs`: vector containing the expansion coefficients, the i-th component is
the i-1 coefficient of the expansion. The vector entries may be TaylorRec
polynomials (in another variable), which is the basis of the TaylorRec
recursive representation
""" ->
immutable TaylorRec{T<:Number, N} <: AbstractSeries{T,N}
    coeffs :: Array

    ## Inner constructors ##
    function TaylorRec(coeffs::Array)
        nn = @compat(Int(N))
        @assert nn > 0
        @assert nn == get_numvars(coeffs)+1
        if !isa(eltype(coeffs), Type{T})
            if nn == 1
                coeffs = convert(Array{T,1}, coeffs)
            else
                coeffs = convert(Array{TaylorRec{T,N-1},1}, coeffs)
            end
        end
        new(coeffs)
    end

    function TaylorRec(c::T)
        nn = @compat(Int(N))
        @assert nn > 0
        nn == 1 && return TaylorRec{T,1}([c])
        c = TaylorRec{T,nn-1}(c)
        return TaylorRec{T,nn}([c])
    end
end


## Outer constructors
TaylorRec{T<:Number,N}(x::TaylorRec{T,N}) = x
TaylorRec{T<:Number,N}(coeffs::Array{TaylorRec{T,N},1}) =
    TaylorRec{T,N+1}(coeffs)
TaylorRec{T<:Number}(coeffs::Array{T,1}) = TaylorRec{T,1}(coeffs)
TaylorRec{T<:Number}(c::T) = TaylorRec{T,1}([c])


## eltype, get_numvars, get_order ##
eltype{T<:Number,N}(::TaylorRec{T,N}) = T
eltype{T<:Number,N}(::Array{TaylorRec{T,N},1}) = T

get_numvars{T<:Number,N}(::TaylorRec{T,N}) = N
get_numvars{T<:Number,N}(v::Array{TaylorRec{T,N},1}) = get_numvars(v[1])
get_numvars{T<:Number}(::Array{T,1}) = 0
get_numvars{T<:Number}(::T) = 0

function get_order{T<:Number,N}(a::TaylorRec{T,N}, var::Int)
    nn = @compat Int(N)
    (var > nn) && return 0
    (nn == 1 || var == 1) && return length(a.coeffs)-1

    ord = 0
    for i in eachindex(a.coeffs)
        tmp = get_order(a.coeffs[i],var-1)
        ord = max(ord, tmp)
    end

    return ord
end
function get_order{T<:Number,N}(a::TaylorRec{T,N})
    nn = @compat Int(N)
    ords = Array(Int,nn)
    for i = 1:nn
        ords[i] = get_order(a, i)
    end
    return (ords...)
end


## taylorRec_variable: creates an independent TaylorRec{T,N}
@doc doc"""
Creates a TaylorRec{T,1} structure corresponding to the 1-independent
variable of type T. `ord` specifies the maximum degree.""" ->
function taylorRec_variable{T<:Number}(::Type{TaylorRec{T,1}}, ord::Int)
    v = zeros(T,ord+1)
    v[2] = one(T)
    TaylorRec{T,1}(v)
end
@doc doc"""
Creates a TaylorRec{T,N} structure corresponding to the N-independent
variable of type T. `ord` specifies the maximum degree of *that*
variable.""" ->
function taylorRec_variable{T<:Number,N}(::Type{TaylorRec{T,N}}, ord::Int)
    @compat nn = Int(N)
    x = taylorRec_variable(TaylorRec{T,1}, ord)
    for n = 2:nn
        x = TaylorRec{T,n}([x])
    end
    x
end

@doc doc"""
Creates a TaylorRec structure corresponding to the i-independent
variable. The `type` (default is `Float64`) and the `order` of the variable
(default is 1) can be specified.""" ->
macro taylorRec(nv, args...)
    local ord = 1 :: Int
    local T = Float64
    if !isempty(args)
        ll = length(args)
        ord = args[1] < 1 ? 1 : args[1] :: Int
        T = ll >=2 ? args[2] : Float64
    end
    :(taylorRec_variable(TaylorRec{$T,$nv}, $ord))
end


## fixorder!
function fixorder!{T<:Number,N}(x::TaylorRec{T,N}, order::Int64)
    ll = length(x.coeffs)
    order+1 <= ll && return TaylorRec{T,N}(x.coeffs)

    resize!(x.coeffs,order+1)
    @inbounds z = zero(typeof(x.coeffs[1]))
    @simd for i=ll+1:order+1
        @inbounds x.coeffs[i] = z
    end

    return x
end


## zero, zeros, one, ones ##
zero{T<:Number}(::Type{TaylorRec{T,1}}) = TaylorRec{T,1}([zero(T)])
function zero{T<:Number, N}(::Type{TaylorRec{T,N}})
    z = zero(TaylorRec{T,N-1})
    TaylorRec{T,N}([z])
end
zero{T<:Number}(a::TaylorRec{T,1}) = TaylorRec{T,1}(zeros(T, length(a.coeffs)))
function zero{T<:Number,N}(a::TaylorRec{T,N})
    ll = length(a.coeffs)
    z = Array(TaylorRec{T,N-1}, ll)
    for i in eachindex(a.coeffs)
        z[i] = zero(a.coeffs[i])
    end
    TaylorRec{T,N}(z)
end

function zeros{T<:Number, N}(TT::Type{TaylorRec{T,N}}, dim::Int64)
    v = Array(TT, dim)
    z = zero(TT)
    @inbounds for i = 1:dim
        v[i] = z
    end
    return v
end
function zeros{T<:Number, N}(a::TaylorRec{T,N}, dim::Int64)
    v = Array(TT, dim)
    @inbounds for i = 1:dim
        v[i] = zero(a)
    end
    return v
end

one{T<:Number}(::Type{TaylorRec{T,1}}) = TaylorRec{T,1}([one(T)])
function one{T<:Number, N}(::Type{TaylorRec{T,N}})
    z = zero(TaylorRec{T,N-1})
    z.coeffs[1] = one(z.coeffs[1])
    TaylorRec{T,N}([z])
end
function one{T<:Number}(a::TaylorRec{T,1})
    z = zero(a)
    z.coeffs[1] = one(T)
    z
end
function one{T<:Number,N}(a::TaylorRec{T,N})
    z = zero(a)
    z.coeffs[1] = one(z.coeffs[1])
    z
end

function ones{T<:Number, N}(TT::Type{TaylorRec{T,N}}, dim::Int64)
    v = Array(TT, dim)
    z = one(TT)
    @inbounds for i = 1:dim
        v[i] = z
    end
    return v
end
function ones{T<:Number, N}(a::TaylorRec{T,N}, dim::Int64)
    v = Array(TT, dim)
    @inbounds for i = 1:dim
        v[i] = one(a)
    end
    return v
end


## Equality ##
function ==(a::TaylorRec, b::TaylorRec)
    a, b = fixshape(a, b)
    test = true
    for i in eachindex(a.coeffs)
        @inbounds test = a.coeffs[i] == b.coeffs[i]
        ~test && return test
    end
    return test
end
==(a::TaylorRec, b::Number) = ==(promote(a,b)...)
==(b::Number, a::TaylorRec) = ==(promote(a,b)...)


## Auxiliary functions: fixshape, firstnonzero ##
function fixshape{T<:Number, S<:Number}(aa::TaylorRec{T,1}, bb::TaylorRec{S,1})
    a = deepcopy(aa)
    b = deepcopy(bb)
    a, b = promote(a, b)
    la = length(a.coeffs)
    lb = length(b.coeffs)
    la == lb && return a, b
    if la < lb
        fixorder!(a, lb-1)
    elseif lb < la
        fixorder!(b, la-1)
    end
    return a, b
end

function fixshape{T<:Number,S<:Number,N,M}(aa::TaylorRec{T,N}, bb::TaylorRec{S,M})
    a = deepcopy(aa)
    b = deepcopy(bb)
    a, b = promote(a, b)
    R = eltype(a)
    nn = get_numvars(a)

    la = length(a.coeffs)
    lb = length(b.coeffs)
    if la < lb
        fixorder!(a, lb-1)
    elseif lb < la
        fixorder!(b, la-1)
    end

    la = length(a.coeffs)
    for i in eachindex(a.coeffs)
        a.coeffs[i], b.coeffs[i] = fixshape(a.coeffs[i], b.coeffs[i])
    end
    return a, b
end

function firstnonzero(a::TaylorRec)
    nonzero = length(a.coeffs)
    z = zero(a.coeffs[1])
    for i in eachindex(a.coeffs)
        if a.coeffs[i] != z
            nonzero = i-1
            break
        end
    end
    return nonzero
end


# ## Conversion and promotion rules ##
convert{T<:Number,N}(::Type{TaylorRec{T,N}}, a::TaylorRec{T,N}) = a
function convert{T<:Number,S<:Number,N,M}(::Type{TaylorRec{T,N}}, a::TaylorRec{S,M})
    nn = @compat(Int(N))
    mm = @compat(Int(M))
    @assert nn >= mm
    a = TaylorRec{T,M}(a.coeffs)
    nn == mm && return a

    for nvar = mm+1:nn
        ll = length(a.coeffs)
        mm = get_numvars(a)
        coeffs = Array(TaylorRec{T,mm}, ll)
        for i in eachindex(a.coeffs)
            coeffs[i] = convert(TaylorRec{T,mm}, a.coeffs[i])
        end
        a = TaylorRec{T,nvar}(coeffs)
    end

    return a
end
convert{T<:Number, S<:Number, N}(::Type{TaylorRec{T,N}}, a::S) =
    convert(TaylorRec{T,N},TaylorRec{T,1}([convert(T,a)]))
convert{T<:Number}(::Type{TaylorRec}, a::T) = TaylorRec{T,1}([a])

function promote_rule{T<:Number,S<:Number,N,M}(::Type{TaylorRec{T,N}},
        ::Type{TaylorRec{S,M}})
    TT = promote_type(T,S)
    nn = max( @compat(Int(N)), @compat(Int(M)) )
    TaylorRec{TT, nn}
end
promote_rule{T<:Number, S<:Number, N}(::Type{TaylorRec{T, N}}, ::Type{S}) =
    TaylorRec{promote_type(T,S),N}


## Addition and substraction ##
for f in (:+, :-)
    @eval begin
        function ($f)(a::TaylorRec, b::TaylorRec)
            a, b = fixshape(a, b)
            v = zero(a.coeffs)
            @inbounds for i in eachindex(a.coeffs)
                v[i] = ($f)(a.coeffs[i], b.coeffs[i])
            end
            return TaylorRec(v)
        end
        function ($f)(a::TaylorRec)
            v = zero(a.coeffs)
            @inbounds for i in eachindex(a.coeffs)
                v[i] = ($f)(a.coeffs[i])
            end
            return TaylorRec(v)
        end
        ($f){T<:Number,S<:Number,N}(a::TaylorRec{T,N}, b::S) = ($f)(promote(a,b)...)
        ($f){T<:Number,S<:Number,N}(b::S, a::TaylorRec{T,N}) = ($f)(promote(b,a)...)
    end
end


## Multiplication ##
function *{T<:Number,S<:Number}(b::S, a::TaylorRec{T,1})
    R = promote_type(T,S)
    coeffs = Array(R, length(a.coeffs))
    @inbounds for k in eachindex(a.coeffs)
        coeffs[k] = b * a.coeffs[k]
    end
    TaylorRec(coeffs)
end
*{T<:Number,S<:Number}(a::TaylorRec{T,1}, b::S) = *(b,a)

function *{T<:Number,S<:Number,N}(b::S, a::TaylorRec{T,N})
    R = promote_type(T,S)
    coeffs = Array(TaylorRec{R,N-1}, length(a.coeffs))
    @inbounds for k in eachindex(a.coeffs)
        coeffs[k] = b * a.coeffs[k]
    end
    TaylorRec(coeffs)
end
*{T<:Number,S<:Number,N}(a::TaylorRec{T,N}, b::S) = *(b,a)

function *{T<:Number, S<:Number, N, M}(aa::TaylorRec{T,N}, bb::TaylorRec{S,M})
    # a, b = fixshape(a, b)
    a, b = fixshape(aa, bb)
    la = length(a.coeffs)
    fixorder!(a, 2*(la-1))
    fixorder!(b, 2*(la-1))
    coeffs = Array(typeof(a.coeffs[1]), 2*la-1)
    @inbounds for k in eachindex(a.coeffs)
        coeffs[k] = mulHomogCoef(k-1, a.coeffs, b.coeffs)
    end
    TaylorRec(coeffs)
end

function mulHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, bc::Array{T,1})
    kcoef == 0 && return ac[1] * bc[1]
    coefhomog = zero(T)
    @inbounds for i = 0:kcoef
        coefhomog += ac[i+1] * bc[kcoef-i+1]
    end
    coefhomog
end

function mulHomogCoef{T<:Number, N}(kcoef::Int, ac::Array{TaylorRec{T,N},1},
    bc::Array{TaylorRec{T,N},1})

    kcoef == 0 && return ac[1] * bc[1]
    coefhomog = zero(ac[1])
    @inbounds for i = 0:kcoef
        coefhomog += ac[i+1] * bc[kcoef-i+1]
    end
    coefhomog
end


## Division ##
function /{T<:Number, S<:Number}(a::TaylorRec{T,1}, b::TaylorRec{S,1})
    a, b = fixshape(a, b)
    la = length(a.coeffs)
    order_ini, c = divfactorization(a, b) # order and coefficient of first factorized term
    TT = typeof(c)
    v1 = convert(Array{TT,1}, a.coeffs)
    v2 = convert(Array{TT,1}, b.coeffs)
    coeffs = zeros(TT, la)
    @inbounds coeffs[1] = c
    @inbounds for k = order_ini+1:la-1
        coeffs[k-order_ini+1] = divHomogCoefRec(k, v1, v2, coeffs, order_ini)
    end
    TaylorRec(coeffs)
end

function /{T<:Number,S<:Number, N, M}(a::TaylorRec{T,N}, b::TaylorRec{S,M})
    b0 = b.coeffs[1]
    @assert b0 != zero(b0)
    a, b = fixshape(a, b)
    la = length(a.coeffs)
    # ord_ini, c = divfactorization(a, b) # order and coefficient of first factorized term
    ord_ini = 0
    @inbounds c = a.coeffs[1] / b0
    TT = typeof(c)
    v1 = convert(Array{TT,1}, a.coeffs)
    v2 = convert(Array{TT,1}, b.coeffs)
    coeffs = zeros(TT, la)
    @inbounds coeffs[1] = c
    @inbounds for k = ord_ini+1:la-1
        coeffs[k-ord_ini+1] = divHomogCoefRec(k, v1, v2, coeffs, ord_ini)
    end
    TaylorRec(coeffs)
end

function /{T<:Number,S<:Union(Real,Complex),N}(a::TaylorRec{T,N}, b::S)
    R = eltype( one(T)/one(S) )
    coeffs = Array(TaylorRec{R,N}, length(a.coeffs))
    # @inbounds coeffs[1] = a.coeffs[1] * b.coeffs[1]
    @inbounds for k in eachindex(a.coeffs)
        coeffs[k] = a.coeffs[k] / b
    end
    TaylorRec(coeffs)
end

/{T<:Number,S<:Union(Real,Complex),N}(b::S, a::TaylorRec{T,N}) = /(TaylorRec(b),a)

function divHomogCoefRec{T<:Number}(kcoef::Int, ac::Array{T,1}, bc::Array{T,1},
    coeffs::Array{T,1}, ordfact=0::Int)

    @inbounds kcoef == ordfact && return ac[ordfact+1] / bc[ordfact+1]
    coefhomog = mulHomogCoef(kcoef, coeffs, bc)
    @inbounds coefhomog = (ac[kcoef+1]-coefhomog) / bc[ordfact+1]
    coefhomog
end
function divHomogCoefRec{T<:Number,N}(kcoef::Int, ac::Array{TaylorRec{T,N},1},
    bc::Array{TaylorRec{T,N},1}, coeffs::Array{TaylorRec{T,N},1}, ordfact=0::Int)

    @inbounds kcoef == ordfact && return ac[ordfact+1] / bc[ordfact+1]
    coefhomog = mulHomogCoef(kcoef, coeffs, bc)
    @inbounds coefhomog = (ac[kcoef+1]-coefhomog) / bc[ordfact+1]
    coefhomog
end

function divfactorization{T<:Number}(a1::TaylorRec{T,1}, b1::TaylorRec{T,1})
    # order of first factorized term; a1 and b1 are assumed to be of the same order (length)
    order = length(a1.coeffs)-1
    a1nz = firstnonzero(a1)
    b1nz = firstnonzero(b1)
    ord_ini = min(a1nz, b1nz)
    if ord_ini > order
        ord_ini = order
    end
    c = a1.coeffs[ord_ini+1] / b1.coeffs[ord_ini+1]
    aux = abs2(c)

    # Is the polynomial factorizable?
    if isinf(aux) || isnan(aux)
        info("Order k=$(ord_ini) => coeff[$(ord_ini+1)]=$(c)")
        error("Division does not define a TaylorRec polynomial\n",
            " or its first non-zero coefficient is Inf/NaN.\n")
    ##else ord_ini>0
    ##    warn("Factorizing the polynomial.\n",
    ##        "The last k=$(ord_ini) TaylorRec coefficients ARE SET to 0.\n")
    end

    return ord_ini, c
end


## real, imag, conj and ctranspose ##
for f in (:real, :imag, :conj)
    @eval begin
        function ($f)(a::TaylorRec)
            r = ($f)(a.coeffs[1])
            v = Array(typeof(r), length(a.coeffs))
            @inbounds for i in eachindex(v)
                v[i] = ($f)(a.coeffs[i])
            end
            return TaylorRec(v)
        end
    end
end
ctranspose(a::TaylorRec) = conj(a)


## Division functions: rem and mod
for op in (:mod, :rem)
    @eval begin
        function ($op){T<:Real,N}(a::TaylorRec{T,N}, x::Real)
            v = deepcopy(a)
            @inbounds v.coeffs[1] = ($op)(v.coeffs[1], x)
            return v
        end
    end
end

function mod2pi{T<:Real,N}(a::TaylorRec{T,N})
    v = deepcopy(a)
    @inbounds v[1] = mod2pi( v[1] )
    return TaylorRec(v)
end



## Int power ##
function ^{T<:Number}(a::TaylorRec{T}, n::Integer)
    n == 0 && return one(a)
    n == 1 && return a
    n == 2 && return square(a)
    n < 0  && return inv( a^(-n) )
    return power_by_squaring(a, n)
end
## Rational power ##
^(a::TaylorRec,x::Rational) = a^(x.num/x.den)
## Real power ##
function ^(a::TaylorRec, x::Real)
    uno = one(a)
    x == zero(x) && return uno
    x == one(x)/2 && return sqrt(a)
    isinteger(x) && return a^(@compat Int(x))
    order = length(a.coeffs)-1

    # First non-zero coefficient
    l0nz = firstnonzero(a)
    l0nz > order && return zero(a)

    # The first non-zero coefficient of the result; must be integer
    lnull = x*l0nz
    !isinteger(lnull) &&
        error("Integer exponent REQUIRED if the TaylorRec polynomial is expanded around 0.\n")

    # Reaching this point, it is possible to implement the power of the TaylorRec polynomial.
    # The last l0nz coefficients are set to zero.
    lnull = trunc(Int,lnull)
    ##l0nz > 0 && warn("The last k=$(l0nz) TaylorRec coefficients ARE SET to 0.\n")
    @inbounds aux = (a.coeffs[l0nz+1])^x
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    coeffs = zeros(T, order+1)
    @inbounds coeffs[lnull+1] = aux
    k0 = lnull+l0nz
    @inbounds for k = k0+1:order
        coeffs[k-l0nz+1] = powHomogCoef(k, v, x, coeffs, l0nz)
    end

    TaylorRec(coeffs)
end
# ^(a::TaylorRec, x::Complex) = exp( x*log(a) )
# ^(a::TaylorRec, b::TaylorRec) = exp( b*log(a) )

## power_by_squaring; slightly modified from base/intfuncs.jl
function power_by_squaring(x::TaylorRec, p::Integer)
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
function powHomogCoef{T<:Number,N}(kcoef::Int, ac::Array{TaylorRec{T,N},1}, x::Real,
    coeffs::Array{TaylorRec{T,N},1}, knull::Int)

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
function square(a::TaylorRec)
    order = length(a.coeffs)-1
    T = typeof(a.coeffs[1])
    coeffs = Array(T,order+1)
    # coeffs = zeros(T,order+1)
    # coeffs[1] = a.coeffs[1]^2
    @inbounds for k in eachindex(coeffs)
        coeffs[k] = squareHomogCoef(k-1, a.coeffs)
    end
    TaylorRec(coeffs)
end

# Homogeneous coefficients for square
function squareHomogCoef{T<:Number,N}(kcoef::Int, ac::Union(Array{T,1},
    Array{TaylorRec{T,N},1}))

    kcoef == 0 && return ac[1]^2
    coefhomog = zero(T)
    kodd = kcoef%2
    kend = div(kcoef - 2 + kodd, 2)
    @inbounds for i = 0:kend
        coefhomog += ac[i+1]*ac[kcoef-i+1]
    end
    coefhomog = 2coefhomog
    if kodd == 0
        @inbounds coefhomog += ac[div(kcoef,2)+1]^2
    end
    coefhomog
end

## Square root ##
function sqrt(a::TaylorRec)
    order = length(a.coeffs)-1
    # First non-zero coefficient
    l0nz = firstnonzero(a)
    if l0nz > order
        return zero(a)
    elseif l0nz%2 == 1 # l0nz must be pair
        error("First non-vanishing TaylorRec coefficient must be an EVEN POWER\n",
            "to expand SQRT around 0.\n")
    end

    # Reaching this point, it is possible to implement the sqrt of the TaylorRec polynomial.
    # The last l0nz coefficients are set to zero.
    ##l0nz > 0 && warn("The last k=$(l0nz) TaylorRec coefficients ARE SET to 0.\n")
    lnull = div(l0nz, 2)
    @inbounds aux = sqrt(a.coeffs[l0nz+1])
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeffs)
    coeffs = zeros(T, order+1)
    @inbounds coeffs[lnull+1] = aux
    @inbounds for k = lnull+1:order-l0nz
        coeffs[k] = sqrtHomogCoef(k-1, v, coeffs, lnull)
    end
    TaylorRec(coeffs)
end

# Homogeneous coefficients for the square-root
function sqrtHomogCoef{T<:Number,N}(kcoef::Int, ac::Union(Array{T,1},
    Array{TaylorRec{T,N},1}),coeffs::Union(Array{T,1},Array{TaylorRec{T,N},1}),knull::Int)

    kcoef == knull && return sqrt(ac[2*knull+1])
    coefhomog = zero(T)
    kodd = (kcoef - knull)%2
    kend = div(kcoef - knull - 2 + kodd, 2)
    @inbounds for i = knull+1:knull+kend
        coefhomog += coeffs[i+1]*coeffs[kcoef+knull-i+1]
    end
    @inbounds aux = ac[kcoef+knull+1]-2coefhomog

    if kodd == 0
        @inbounds aux = aux - (coeffs[kend+knull+2])^2
    end
    @inbounds coefhomog = aux / (2coeffs[knull+1])

    coefhomog
end



#------- Functions to get niceR printing
function pretty_print(a::TaylorRec)
    z = zero(eltype(a))
    space = utf8(" ")
    a == zero(a) && return string(space, z)
    strout::UTF8String = ""
    ifirst = true
    for i in eachindex(a.coeffs)
        pol = a.coeffs[i]
        pol == zero(pol) && continue
        cadena::UTF8String = taylorrec2str( pol, i-1, 1 )
        strsgn = (ifirst || i == 1 || cadena[2] == '-') ? string("") : string(" +")
        strout = string( strout, strsgn, cadena)
        ifirst = false
    end
    strout
end

function monomial(ord::Int, nv::Int)
    monom::UTF8String = string("")
    if ord == 0
        monom = string(monom)
    elseif ord == 1
        monom = string(monom, name_taylorNvar(nv))
    else
        monom = string(monom, name_taylorNvar(nv), superscriptify(ord))
    end
    monom
end

function taylorrec2str{T<:Number}(c::T, ord::Int, nv::Int)
    z = zero(T)
    space = utf8(" ")
    c == z && return string( space, z)
    strout::UTF8String = space
    ifirst = true

    monom = monomial(ord,nv)
    cadena = numbr2str(c, ifirst)
    strout = string(strout, cadena, monom)

    return strout
end

function taylorrec2str{T<:Number,N}(a::TaylorRec{T,N}, ord::Int, nv::Int)
    z = zero(T)
    space = utf8(" ")
    a == zero(a) && return string( space, z)
    strout::UTF8String = string("")
    ifirst = true

    monom = monomial(ord,nv)
    for i in eachindex(a.coeffs)
        pol = a.coeffs[i]
        pol == zero(pol) && continue
        cadena::UTF8String = taylorrec2str( pol, i-1, nv+1 )
        strsgn = (ifirst || i == 1 || cadena[2] == '-') ? string("") : string(" +")
        strout = string( strout, strsgn, cadena)
        ifirst = false
    end
    strout = string(space,"(",strout," )", monom)
    strout
end
