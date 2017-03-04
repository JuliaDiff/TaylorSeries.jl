# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

# Arithmetic operations: +, -, *, /

## Equality ##
for T in (:Taylor1, :TaylorN)
    @eval begin
        =={T<:Number,S<:Number}(a::$T{T}, b::$T{S}) = ==(promote(a,b)...)
        function =={T<:Number}(a::$T{T}, b::$T{T})
            a, b = fixorder(a, b)
            @simd for i in eachindex(a.coeffs)
                @inbounds (a.coeffs[i] == b.coeffs[i]) || return false
            end
            true
        end
    end
end

==(a::HomogeneousPolynomial, b::HomogeneousPolynomial) = a.coeffs == b.coeffs


function iszero{T}(a::HomogeneousPolynomial{T})
    @simd for i in eachindex(a.coeffs)
        @inbounds (a.coeffs[i] == zero(T)) || return false
    end
    return true
end


## zero and one ##
zero{T<:Number}(a::Taylor1{T}) = Taylor1(zero(a.coeffs[1]), a.order)
one{T<:Number}(a::Taylor1{T}) = Taylor1(one(a.coeffs[1]), a.order)


function zero{T<:Number}(a::HomogeneousPolynomial{T})
    a.order == 0 && return HomogeneousPolynomial([zero(a.coeffs[1])], 0)
    v = Array{T}( size_table[a.order+1] )
    @simd for ind in eachindex(a.coeffs)
        @inbounds v[ind] = zero(a.coeffs[ind])
    end
    return HomogeneousPolynomial(v, a.order)
end

function zeros{T<:Number}(a::HomogeneousPolynomial{T}, order::Int)
    order == 0 && return [HomogeneousPolynomial([zero(a.coeffs[1])], 0)]
    v = Array{HomogeneousPolynomial{T}}(order+1)
    @simd for ord in eachindex(v)
        @inbounds v[ord] = HomogeneousPolynomial([zero(a.coeffs[1])], ord-1)
    end
    return v
end

zeros{T<:Number}(::Type{HomogeneousPolynomial{T}}, order::Int) =
    zeros( HomogeneousPolynomial([zero(T)], 0), order)

function one{T<:Number}(a::HomogeneousPolynomial{T})
    a.order == 0 && return HomogeneousPolynomial([one(a.coeffs[1])], 0)
    v = Array{T}( size_table[a.order+1] )
    @simd for ind in eachindex(a.coeffs)
        @inbounds v[ind] = one(a.coeffs[ind])
    end
    return HomogeneousPolynomial{T}(v, a.order)
end

function ones{T<:Number}(a::HomogeneousPolynomial{T}, order::Int)
    order == 0 && return [HomogeneousPolynomial([one(a.coeffs[1])], 0)]
    v = Array{HomogeneousPolynomial{T}}(order+1)
    @simd for ord in eachindex(v)
        @inbounds num_coeffs = size_table[ord]
        vT = Array{T}(num_coeffs)
        @simd for ind in 1:num_coeffs
            @inbounds vT[ind] = one(a.coeffs[1])
        end
        @inbounds v[ord] = HomogeneousPolynomial(vT, ord-1)
    end
    return v
end

ones{T<:Number}(::Type{HomogeneousPolynomial{T}}, order::Int) =
    ones( HomogeneousPolynomial([one(T)], 0), order)

zero{T<:Number}(a::TaylorN{T}) = TaylorN(zero(T), a.order)
one{T<:Number}(a::TaylorN{T}) = TaylorN(one(T), a.order)



## Addition and substraction ##
for f in (:+, :-)
    for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
        @eval begin
            ($f){T<:Number,S<:Number}(a::$T{T}, b::$T{S}) = $f(promote(a,b)...)
            function $f{T<:Number}(a::$T{T}, b::$T{T})
                a, b = fixorder(a, b)
                v = similar(a.coeffs)
                @simd for i in eachindex(v)
                    @inbounds v[i] = $f(a.coeffs[i], b.coeffs[i])
                end
                return $T(v, a.order)
            end
            function $f(a::$T)
                v = similar(a.coeffs)
                @simd for i in eachindex(v)
                    @inbounds v[i] = $f(a.coeffs[i])
                end
                return $T(v, a.order)
            end
            if $T != HomogeneousPolynomial
                function $f(a::$T, b::RealOrComplex)
                    R = promote_type(eltype(a.coeffs), eltype(b))
                    coeffs = Array{R}(a.order+1)
                    @simd for i in eachindex(coeffs)
                        @inbounds coeffs[i] = a.coeffs[i]
                    end
                    @inbounds coeffs[1] = $f(a.coeffs[1], b)
                    return $T(coeffs, a.order)
                end
                function $f(b::RealOrComplex, a::$T)
                    R = promote_type(eltype(a.coeffs), eltype(b))
                    coeffs = Array{R}(a.order+1)
                    @simd for i in eachindex(coeffs)
                        @inbounds coeffs[i] = $f(a.coeffs[i])
                    end
                    @inbounds coeffs[1] = $f(b, a.coeffs[1])
                    return $T(coeffs, a.order)
                end
            end
        end
    end
    @eval begin
        function ($f){T<:RealOrComplex,S<:RealOrComplex}(
                a::TaylorN{Taylor1{T}}, b::Taylor1{S})
            @inbounds aux = $f(a.coeffs[1].coeffs[1], b)
            R = typeof(aux)
            coeffs = Array{HomogeneousPolynomial{Taylor1{R}}}(length(a.coeffs))
            @simd for i in eachindex(coeffs)
                @inbounds coeffs[i] = a.coeffs[i]
            end
            @inbounds coeffs[1] = aux
            return TaylorN(coeffs, a.order)
        end
        function ($f){T<:RealOrComplex,S<:RealOrComplex}(
            b::Taylor1{S}, a::TaylorN{Taylor1{T}})
            @inbounds aux = $f(b, a.coeffs[1].coeffs[1])
            R = typeof(aux)
            coeffs = Array{HomogeneousPolynomial{Taylor1{R}}}(length(a.coeffs))
            @simd for i in eachindex(coeffs)
                @inbounds coeffs[i] = $f(a.coeffs[i])
            end
            @inbounds coeffs[1] = aux
            return TaylorN(coeffs, a.order)
        end
    end
end



## Multiplication ##
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval begin
        *(a::Bool, b::$T) = *(convert(Int, a), b)
        *(a::$T, b::Bool) = b * a
        function *(a::RealOrComplex, b::$T)
            @inbounds aux = a * b.coeffs[1]
            v = Array{typeof(aux)}(length(b.coeffs))
            @simd for i in eachindex(v)
                @inbounds v[i] = a * b.coeffs[i]
            end
            $T(v, b.order)
        end
        *(b::$T, a::RealOrComplex) = a * b
        if $T != Taylor1
            function *{T<:RealOrComplex,R<:RealOrComplex}(
                    a::Taylor1{R}, b::$T{Taylor1{T}})
                @inbounds aux = a * b.coeffs[1]
                S = typeof(aux)
                coeffs = Array{S}(length(a.coeffs))
                @simd for i in eachindex(coeffs)
                    @inbounds coeffs[i] = a * b.coeffs[i]
                end
                return $T(coeffs, a.order)
            end
            *{T<:RealOrComplex,R<:RealOrComplex}(
                    b::$T{Taylor1{T}}, a::Taylor1{R}) = a * b
        end
    end
end


doc"""
```
*(a, b)
```

Return the Taylor expansion of $a \cdot b$, of order `max(a.order,b.order)`, for
`a::Taylor1`, `b::Taylor1` polynomials.

For details on making the Taylor expansion, see [`TaylorSeries.mulHomogCoef`](@ref).
"""
*{T<:Number,S<:Number}(a::Taylor1{T}, b::Taylor1{S}) = *(promote(a,b)...)
function *{T<:Number}(a::Taylor1{T}, b::Taylor1{T})
    a, b = fixorder(a, b)
    coeffs = similar(a.coeffs)
    @inbounds coeffs[1] = a.coeffs[1] * b.coeffs[1]
    @inbounds for k = 1:a.order
        coeffs[k+1] = mulHomogCoef(k, a.coeffs, b.coeffs)
    end
    Taylor1(coeffs, a.order)
end

*{T<:Number,S<:Number}(a::TaylorN{T}, b::TaylorN{S}) = *(promote(a,b)...)
function *{T<:Number}(a::TaylorN{T}, b::TaylorN{T})
    a, b = fixorder(a, b)
    coeffs = zeros(HomogeneousPolynomial{T}, a.order)

    for ord in eachindex(coeffs)
        for i = 0:ord-1
            @inbounds mul!(coeffs[ord], a.coeffs[i+1], b.coeffs[ord-i])
        end
    end

    return TaylorN(coeffs, a.order)
end

## Multiplication ##
*{T<:Number,S<:Number}(a::HomogeneousPolynomial{T}, b::HomogeneousPolynomial{S}) =
    *(promote(a,b)...)
function *{T<:Number}(a::HomogeneousPolynomial{T}, b::HomogeneousPolynomial{T})
    order = a.order + b.order
    order > get_order() && return HomogeneousPolynomial([zero(T)], get_order())

    res = HomogeneousPolynomial([zero(T)], order)
    mul!(res, a, b)
    return res
end



# Homogeneous coefficient for the multiplication
doc"""
    mulHomogCoef(kcoef, ac, bc)

Compute the `k`-th expansion coefficient of $c = a\cdot b$ given by

\begin{equation*}
c_k = \sum_{j=0}^k a_j b_{k-j},
\end{equation*}

with $a$ and $b$ `Taylor1` polynomials.

Inputs are the `kcoef`-th coefficient, and the vectors of the expansion coefficients
`ac` and `bc`, corresponding respectively to `a` and `b`.
"""
function mulHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, bc::Array{T,1})
    kcoef == 0 && return ac[1] * bc[1]
    coefhomog = zero(T)
    @inbounds for i = 0:kcoef
        coefhomog += ac[i+1] * bc[kcoef-i+1]
    end
    coefhomog
end


"""
```
mul!(c, a, b)
```

Return `c = a*b` with no allocation; all parameters are `HomogeneousPolynomial`.
"""
function mul!(c::HomogeneousPolynomial, a::HomogeneousPolynomial, b::HomogeneousPolynomial)
    (iszero(b) || iszero(a)) && return nothing

    T = eltype(c)
    @inbounds num_coeffs_a = size_table[a.order+1]
    @inbounds num_coeffs_b = size_table[b.order+1]
    @inbounds num_coeffs  = size_table[c.order+1]

    @inbounds posTb = pos_table[c.order+1]

    @inbounds for na = 1:num_coeffs_a
        ca = a.coeffs[na]
        ca == zero(T) && continue
        inda = index_table[a.order+1][na]

        @inbounds for nb = 1:num_coeffs_b
            cb = b.coeffs[nb]
            cb == zero(T) && continue
            indb = index_table[b.order+1][nb]

            pos = posTb[inda + indb]
            c.coeffs[pos] += ca * cb
        end
    end

    return nothing
end

"""
    mul!(c, a)

Return `c = a*a` with no allocation; all parameters are `HomogeneousPolynomial`.
"""
function mul!(c::HomogeneousPolynomial, a::HomogeneousPolynomial)
    iszero(a) && return nothing

    T = eltype(c)
    @inbounds num_coeffs_a = size_table[a.order+1]
    @inbounds num_coeffs  = size_table[c.order+1]

    @inbounds posTb = pos_table[c.order+1]

    for na = 1:num_coeffs_a
        @inbounds ca = a.coeffs[na]
        ca == zero(T) && continue
        @inbounds inda = index_table[a.order+1][na]
        @inbounds pos = posTb[2*inda]
        @inbounds c.coeffs[pos] += ca * ca
        @inbounds for nb = na+1:num_coeffs_a
            cb = a.coeffs[nb]
            cb == zero(T) && continue
            indb = index_table[a.order+1][nb]
            pos = posTb[inda+indb]
            c.coeffs[pos] += 2 * ca * cb
        end
    end

    return nothing
end



## Division ##
function /{T<:Integer, S<:RealOrComplex}(a::Taylor1{Rational{T}}, b::S)
    R = typeof( a.coeffs[1] // b)
    v = Array{R}(length(a.coeffs))
    @simd for i in eachindex(v)
        @inbounds v[i] = a.coeffs[i] // b
    end
    Taylor1(v, a.order)
end

for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval begin
        /{T<:RealOrComplex}(a::$T{T}, b::T) = a * inv(b)
        function /{T<:RealOrComplex}(a::$T, b::T)
            R = promote_type(eltype(a), T)
            convert($T{R}, a) * inv(convert(R, b))
        end
        if $T != Taylor1
            /{T<:RealOrComplex,R<:RealOrComplex}(
                a::$T{Taylor1{T}}, x::Taylor1{R}) = a * inv(x)
        end
    end
end


doc"""
```
/(a, b)
```

Return the Taylor expansion of $a/b$, of order `max(a.order,b.order)`, for
`a::Taylor1`, `b::Taylor1` polynomials.

For details on making the Taylor expansion, see
[`TaylorSeries.divHomogCoef`](@ref).
"""
/{T<:Number,S<:Number}(a::Taylor1{T}, b::Taylor1{S}) = /(promote(a,b)...)
function /{R<:Number}(a::Taylor1{R}, b::Taylor1{R})
    a, b = fixorder(a, b)
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

/{T<:Number,S<:Number}(a::TaylorN{T}, b::TaylorN{S}) = /(promote(a,b)...)
function /{R<:Number}(a::TaylorN{R}, b::TaylorN{R})
    @inbounds b0 = b.coeffs[1].coeffs[1]
    @assert b0 != zero(b0)
    a, b = fixorder(a, b)
    # order and coefficient of first factorized term
    # orddivfact, cdivfact = divfactorization(a, b)
    b0 = inv(b0)
    @inbounds cdivfact = a.coeffs[1] * b0
    T = eltype(cdivfact)
    coeffs = zeros(HomogeneousPolynomial{T}, a.order)
    @inbounds coeffs[1] = cdivfact

    for ord in eachindex(coeffs)
        ord == a.order+1 && continue
        @inbounds for i = 0:ord-1
            mul!(coeffs[ord+1], coeffs[i+1], b.coeffs[ord-i+1])
        end
        @inbounds coeffs[ord+1] = (a.coeffs[ord+1] - coeffs[ord+1]) * b0
    end

    return TaylorN(coeffs, a.order)
end


function divfactorization(a1::Taylor1, b1::Taylor1)
    # order of first factorized term; a1 and b1 assumed to be of the same order
    a1nz = findfirst(a1)
    b1nz = findfirst(b1)
    orddivfact = min(a1nz, b1nz)
    if orddivfact < 0
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
doc"""
    divHomogCoef(kcoef, ac, bc, coeffs, ordfact)

Compute the `k-th` expansion coefficient of $c = a / b$ given by

\begin{equation*}
c_k =  \frac{1}{b_0} (a_k - \sum_{j=0}^{k-1} c_j b_{k-j}),
\end{equation*}

with $a$ and $b$ `Taylor1` polynomials.

Inputs are the `kcoef`-th coefficient, the vectors of the expansion coefficients
`ac` and `bc`, corresponding respectively to `a` and `b`, the
already calculated expansion coefficients `coeffs` of `c`, and `ordfact`
which is the order of the factorized term of the denominator,
whenever `b_0` is zero.
"""
function divHomogCoef{T<:Number}(kcoef::Int, ac::Array{T,1}, bc::Array{T,1},
    coeffs::Array{T,1}, ordfact::Int)
    #
    @inbounds kcoef == ordfact && return ac[ordfact+1] / bc[ordfact+1]
    coefhomog = mulHomogCoef(kcoef, coeffs, bc)
    @inbounds coefhomog = (ac[kcoef+1]-coefhomog) / bc[ordfact+1]
    coefhomog
end

## TODO: Implement factorization (divfactorization) for TaylorN polynomials






"""
    A_mul_B!(Y, A, B)

Multiply A*B and save the result in Y.
"""
function A_mul_B!{T<:Number}(y::Vector{Taylor1{T}},
        a::Union{Matrix{T},SparseMatrixCSC{T}},
        b::Vector{Taylor1{T}})

    n, k = size(a)
    @assert (length(y)== n && length(b)== k)

    # determine the maximal order of b
    order = maximum([b1.order for b1 in b])

    # Use matrices of coefficients (of proper size) and A_mul_B!
    B = zeros(T, k, order+1)
    for i = 1:k
        @inbounds ord = b[i].order
        @inbounds for j = 1:ord+1
            B[i,j] = b[i].coeffs[j]
        end
    end
    Y = Array{T}(n, order+1)
    A_mul_B!(Y, a, B)
    @inbounds for i = 1:n
        y[i] = Taylor1( collect(Y[i,:]), order)
    end

    return y
end
