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
            la = a.order+1
            lb = b.order+1
            if a.order == b.order
                return all( a.coeffs .== b.coeffs )
            elseif a.order < b.order
                res1 = all( a[1:end] .== b[1:la] )
                res2 = iszero(b[la+1:end])
            else a.order > b.order
                res1 = all( a[1:lb] .== b[1:end] )
                res2 = iszero(a[lb+1:end])
            end
            return res1 && res2
        end
    end
end

==(a::HomogeneousPolynomial, b::HomogeneousPolynomial) = a.coeffs == b.coeffs


for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval iszero(a::$T) = iszero(a.coeffs)
end

## zero and one ##
for T in (:Taylor1, :TaylorN), f in (:zero, :one)
    @eval begin
        ($f)(a::$T) = $T(($f)(a[1]), a.order)

        ($f)(a::$T, order) = $T(($f)(a[1]), order)
    end
end

function zero{T<:Number}(a::HomogeneousPolynomial{T})
    a.order == 0 && return HomogeneousPolynomial(zero(T), 0)
    v = Array{T}( size_table[a.order+1] )
    v .= zero.(a.coeffs)
    return HomogeneousPolynomial(v, a.order)
end

function zeros{T<:Number}(a::HomogeneousPolynomial{T}, order::Int)
    order == 0 && return [HomogeneousPolynomial(zero(T), 0)]
    v = Array{HomogeneousPolynomial{T}}(order+1)
    @simd for ord in eachindex(v)
        @inbounds v[ord] = HomogeneousPolynomial(zero(T), ord-1)
    end
    return v
end

zeros{T<:Number}(::Type{HomogeneousPolynomial{T}}, order::Int) =
    zeros( HomogeneousPolynomial([zero(T)], 0), order)

function one{T<:Number}(a::HomogeneousPolynomial{T})
    a.order == 0 && return HomogeneousPolynomial([one(a[1])], 0)
    v = Array{T}( size_table[a.order+1] )
    v .= one.(a.coeffs)
    return HomogeneousPolynomial{T}(v, a.order)
end

function ones{T<:Number}(a::HomogeneousPolynomial{T}, order::Int)
    order == 0 && return [HomogeneousPolynomial([one(a[1])], 0)]
    v = Array{HomogeneousPolynomial{T}}(order+1)
    @simd for ord in eachindex(v)
        @inbounds num_coeffs = size_table[ord]
        @inbounds v[ord] = HomogeneousPolynomial(ones(T, num_coeffs), ord-1)
    end
    return v
end

ones{T<:Number}(::Type{HomogeneousPolynomial{T}}, order::Int) =
    ones( HomogeneousPolynomial([one(T)], 0), order)



## Addition and substraction ##
for (f, fc) in ((:+, :(add!)), (:-, :(subst!)))
    dotf = Symbol(".",f)

    for T in (:Taylor1, :TaylorN)
        @eval begin
            ($f){T<:Number,S<:Number}(a::$T{T}, b::$T{S}) = $f(promote(a,b)...)

            function $f{T<:Number}(a::$T{T}, b::$T{T})
                la = a.order+1
                lb = b.order+1
                if a.order == b.order
                    v = similar(a.coeffs)
                    v .= $dotf(a.coeffs, b.coeffs)
                elseif a.order < b.order
                    v = similar(b.coeffs)
                    @inbounds v[1:la] .= $dotf(a[1:end], b[1:la])
                    @inbounds v[la+1:end] .= $f(b[la+1:end])
                else
                    v = similar(a.coeffs)
                    @inbounds v[1:lb] .= $dotf(a[1:lb], b[1:end])
                    @inbounds v[lb+1:end] .= a[lb+1:end]
                end
                return $T(v)
            end

            function $f(a::$T)
                v = similar(a.coeffs)
                broadcast!($f, v, a.coeffs)
                return $T(v, a.order)
            end

            ($f){T<:Number,S<:Number}(a::$T{T}, b::S) = $f(promote(a,b)...)

            function $f{T<:Number}(a::$T{T}, b::T)
                coeffs = similar(a.coeffs)
                coeffs .= a.coeffs
                @inbounds coeffs[1] = $f(a[1], b)
                return $T(coeffs, a.order)
            end

            ($f){T<:Number,S<:Number}(b::S, a::$T{T}) = $f(promote(b,a)...)

            function $f{T<:Number}(b::T, a::$T{T})
                coeffs = similar(a.coeffs)
                coeffs .= $f(a.coeffs)
                @inbounds coeffs[1] = $f(b, a[1])
                return $T(coeffs, a.order)
            end

            ## add! and subst! ##
            function ($fc)(v::$T, a::$T, b::$T, k::Int)
                @inbounds v[k+1] = ($f)(a[k+1], b[k+1])
                return nothing
            end
        end
    end

    @eval begin
        ($f){T<:NumberNotSeriesN,S<:NumberNotSeriesN}(a::HomogeneousPolynomial{T},
            b::HomogeneousPolynomial{S}) = $f(promote(a,b)...)

        function $f{T<:NumberNotSeriesN}(a::HomogeneousPolynomial{T},
                b::HomogeneousPolynomial{T})
            @assert a.order == b.order
            v = similar(a.coeffs)
            v .= $dotf(a.coeffs, b.coeffs)
            return HomogeneousPolynomial(v, a.order)
        end

        function $f(a::HomogeneousPolynomial)
            v = similar(a.coeffs)
            v .= $f(a.coeffs)
            return HomogeneousPolynomial(v, a.order)
        end

        function ($f){T<:NumberNotSeries,S<:NumberNotSeries}(
                a::TaylorN{Taylor1{T}}, b::Taylor1{S})
            @inbounds aux = $f(a[1][1], b)
            R = eltype(aux)
            coeffs = Array{HomogeneousPolynomial{Taylor1{R}}}(a.order+1)
            coeffs .= a.coeffs
            @inbounds coeffs[1] = aux
            return TaylorN(coeffs, a.order)
        end

        function ($f){T<:NumberNotSeries,S<:NumberNotSeries}(
                b::Taylor1{S}, a::TaylorN{Taylor1{T}})
            @inbounds aux = $f(b, a[1][1])
            R = eltype(aux)
            coeffs = Array{HomogeneousPolynomial{Taylor1{R}}}(a.order+1)
            coeffs .= $f(a.coeffs)
            @inbounds coeffs[1] = aux
            return TaylorN(coeffs, a.order)
        end

        function ($f){T<:NumberNotSeries,S<:NumberNotSeries}(
                a::Taylor1{TaylorN{T}}, b::TaylorN{S})
            @inbounds aux = $f(a[1], b)
            R = eltype(aux)
            coeffs = Array{TaylorN{R}}(a.order+1)
            coeffs .= a.coeffs
            @inbounds coeffs[1] = aux
            return Taylor1(coeffs, a.order)
        end

        function ($f){T<:NumberNotSeries,S<:NumberNotSeries}(
                b::TaylorN{S}, a::Taylor1{TaylorN{T}})
            @inbounds aux = $f(b, a[1])
            R = eltype(aux)
            coeffs = Array{TaylorN{R}}(a.order+1)
            coeffs .= $f(a.coeffs)
            @inbounds coeffs[1] = aux
            return Taylor1(coeffs, a.order)
        end
    end
end



## Multiplication ##
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)

    @eval begin
        *(a::Bool, b::$T) = *(convert(Int, a), b)

        *(a::$T, b::Bool) = b * a

        function *{T<:NumberNotSeries}(a::T, b::$T)
            @inbounds aux = a * b[1]
            v = Array{typeof(aux)}(length(b.coeffs))
            v .= a .* b.coeffs
            $T(v, b.order)
        end

        *{T<:NumberNotSeries}(b::$T, a::T) = a * b
    end
end

for T in (:HomogeneousPolynomial, :TaylorN)

    @eval begin
        function *{T<:NumberNotSeries,S<:NumberNotSeries}(a::Taylor1{T}, b::$T{Taylor1{S}})
            @inbounds aux = a * b[1]
            R = typeof(aux)
            coeffs = Array{R}(length(b.coeffs))
            coeffs .= a .* b.coeffs
            return $T(coeffs, b.order)
        end

        *{T<:NumberNotSeries,R<:NumberNotSeries}(b::$T{Taylor1{R}}, a::Taylor1{T}) = a * b

        function *{T<:NumberNotSeries,S<:NumberNotSeries}(a::$T{T}, b::Taylor1{$T{S}})
            @inbounds aux = a * b[1]
            R = typeof(aux)
            coeffs = Array{R}(length(b.coeffs))
            coeffs .= a .* b.coeffs
            return Taylor1(coeffs, b.order)
        end

        *{T<:NumberNotSeries,S<:NumberNotSeries}(b::Taylor1{$T{S}}, a::$T{T}) = a * b
    end
end

for (T, W) in ((:Taylor1, :Number), (:TaylorN, :NumberNotSeriesN))
    @eval *{T<:$W, S<:$W}(a::$T{T}, b::$T{S}) = *(promote(a,b)...)

    @eval function *{T<:$W}(a::$T{T}, b::$T{T})
        corder = max(a.order, b.order)
        c = zero(a, corder)
        @inbounds for ord = 0:corder
            mul!(c, a, b, ord) # updates c[ord+1]
        end
        return c
    end
end


*{T<:NumberNotSeriesN,S<:NumberNotSeriesN}(a::HomogeneousPolynomial{T},
    b::HomogeneousPolynomial{S}) = *(promote(a,b)...)

function *{T<:NumberNotSeriesN}(a::HomogeneousPolynomial{T}, b::HomogeneousPolynomial{T})

    order = a.order + b.order

    order > get_order() && return HomogeneousPolynomial(zero(T), get_order())

    res = HomogeneousPolynomial(zero(T), order)
    mul!(res, a, b)
    return res
end


# Homogeneous coefficient for the multiplication
for T in (:Taylor1, :TaylorN)
    @eval function mul!(c::$T, a::$T, b::$T, k::Int)
        a == b && return sqr!(c, a, k)

        if $T == Taylor1
            c[k+1] = zero(eltype(c))
        else
            c[k+1] = HomogeneousPolynomial(zero(eltype(c)), k)
        end

        @inbounds for i = 0:k
            c[k+1] += a[i+1] * b[k-i+1]
        end

        return nothing
    end
end

doc"""
    mul!(c, a, b, k::Int) --> nothing

Update the `k`-th expansion coefficient `c[k+1]` of `c = a * b`,
where all `c`, `a`, and `b` are either `Taylor1` or `TaylorN`.

The coefficients are given by

\begin{equation*}
c_k = \sum_{j=0}^k a_j b_{k-j}.
\end{equation*}

""" mul!


"""
    mul!(c, a, b) --> nothing

Return `c = a*b` with no allocation; all arguments are `HomogeneousPolynomial`.

"""
function mul!(c::HomogeneousPolynomial, a::HomogeneousPolynomial,
        b::HomogeneousPolynomial)

    (iszero(b) || iszero(a)) && return nothing

    T = eltype(c)
    @inbounds num_coeffs_a = size_table[a.order+1]
    @inbounds num_coeffs_b = size_table[b.order+1]
    @inbounds num_coeffs  = size_table[c.order+1]

    @inbounds posTb = pos_table[c.order+1]

    @inbounds for na = 1:num_coeffs_a
        ca = a[na]
        ca == zero(T) && continue
        inda = index_table[a.order+1][na]

        @inbounds for nb = 1:num_coeffs_b
            cb = b[nb]
            cb == zero(T) && continue
            indb = index_table[b.order+1][nb]

            pos = posTb[inda + indb]
            c[pos] += ca * cb
        end
    end

    return nothing
end



## Division ##
function /{T<:Integer, S<:NumberNotSeries}(a::Taylor1{Rational{T}}, b::S)
    R = typeof( a[1] // b)
    v = Array{R}(a.order+1)
    v .= a.coeffs .// b
    return Taylor1(v, a.order)
end

for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)

    @eval function /{T<:NumberNotSeries,S<:NumberNotSeries}(a::$T{T}, b::S)
        R = promote_type(T,S)
        return convert($T{R}, a) * inv(convert(R, b))
    end

    @eval /{T<:NumberNotSeries}(a::$T, b::T) = a * inv(b)
end

for T in (:HomogeneousPolynomial, :TaylorN)
    @eval function /{T<:NumberNotSeries,S<:NumberNotSeries}(
            b::$T{Taylor1{S}}, a::Taylor1{T})
        @inbounds aux = b[1] / a
        R = typeof(aux)
        coeffs = Array{R}(length(b.coeffs))
        coeffs .= b.coeffs ./ a
        return $T(coeffs, b.order)
    end

    @eval function /{T<:NumberNotSeries,S<:NumberNotSeries}(
            b::Taylor1{$T{S}}, a::$T{T})
        @inbounds aux = b[1] / a
        R = typeof(aux)
        coeffs = Array{R}(length(b.coeffs))
        coeffs .= b.coeffs ./ a
        return Taylor1(coeffs, b.order)
    end
end



/{T<:Number,S<:Number}(a::Taylor1{T}, b::Taylor1{S}) = /(promote(a,b)...)

function /{T<:Number}(a::Taylor1{T}, b::Taylor1{T})
    corder = max(a.order, b.order)

    # order and coefficient of first factorized term
    ordfact, cdivfact = divfactorization(a, b)

    c = Taylor1([cdivfact], corder)
    @inbounds for ord = 1:corder-ordfact
        div!(c, a, b, ord, ordfact) # updates c[ord+1]
    end

    return c
end

/{T<:NumberNotSeriesN,S<:NumberNotSeriesN}(a::TaylorN{T}, b::TaylorN{S}) =
    /(promote(a,b)...)

function /{T<:NumberNotSeriesN}(a::TaylorN{T}, b::TaylorN{T})
    @assert !iszero(constant_term(b))
    corder = max(a.order, b.order)

    # first coefficient
    @inbounds cdivfact = a[1] / constant_term(b)
    c = TaylorN(cdivfact, corder)
    for ord in 1:corder
        div!(c, a, b, ord) # updates c[ord+1]
    end

    return c
end


function divfactorization(a1::Taylor1, b1::Taylor1)
    # order of first factorized term; a1 and b1 assumed to be of the same order
    a1nz = findfirst(a1)
    b1nz = findfirst(b1)
    a1nz = a1nz ≥ 0 ? a1nz : a1.order
    b1nz = b1nz ≥ 0 ? b1nz : a1.order
    ordfact = min(a1nz, b1nz)
    cdivfact = a1[ordfact+1] / b1[ordfact+1]

    # Is the polynomial factorizable?
    if isinf(cdivfact) || isnan(cdivfact)
        throw(ArgumentError(
        """Division does not define a Taylor1 polynomial
        or its first non-zero coefficient is Inf/NaN.
        Order k=$(ordfact) => coeff[$(ordfact+1)]=$(cdivfact)."""))
    end

    return ordfact, cdivfact
end


## TODO: Implement factorization (divfactorization) for TaylorN polynomials


# Homogeneous coefficient for the division
doc"""
    div!(c, a, b, k::Int, ordfact::Int=0)

Compute the `k-th` expansion coefficient `c[k+1]` of `c = a / b`,
where all `c`, `a` and `b` are either `Taylor1` or `TaylorN`.

The coefficients are given by

\begin{equation*}
c_k =  \frac{1}{b_0} (a_k - \sum_{j=0}^{k-1} c_j b_{k-j}).
\end{equation*}

For `Taylor1` polynomials, `ordfact` is the order of the factorized
term of the denominator.
"""
function div!(c::Taylor1, a::Taylor1, b::Taylor1, k::Int, ordfact::Int=0)
    if k == 0
        @inbounds c[1] = a[ordfact+1] / b[ordfact+1]
        return nothing
    end

    coef = zero(constant_term(c))
    @inbounds for i = 0:k
        coef += c[i+1] * b[k+ordfact-i+1]
    end
    @inbounds c[k+1] = (a[k+ordfact+1]-coef) / b[ordfact+1]
    return nothing
end

function div!(c::TaylorN, a::TaylorN, b::TaylorN, k::Int)
    if k==0
        @inbounds c[1] = a[1] / constant_term(b)
        return nothing
    end

    c[k+1] = HomogeneousPolynomial(zero(constant_term(c)), k)
    @inbounds for i = 0:k-1
        mul!(c[k+1], c[i+1], b[k-i+1])
    end
    @inbounds c[k+1] = (a[k+1] - c[k+1]) / constant_term(b)
    return nothing
end







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
            B[i,j] = b[i][j]
        end
    end
    Y = Array{T}(n, order+1)
    A_mul_B!(Y, a, B)
    @inbounds for i = 1:n
        y[i] = Taylor1( collect(Y[i,:]), order)
    end

    return y
end
