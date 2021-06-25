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
        ==(a::$T{T}, b::$T{S}) where {T<:Number,S<:Number} = ==(promote(a,b)...)

        function ==(a::$T{T}, b::$T{T}) where {T<:Number}
            if a.order != b.order
                a, b = fixorder(a, b)
            end
            return a.coeffs == b.coeffs
        end
    end
end

function ==(a::HomogeneousPolynomial, b::HomogeneousPolynomial)
    a.order == b.order && return a.coeffs == b.coeffs
    return iszero(a.coeffs) && iszero(b.coeffs)
end

for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval iszero(a::$T) = iszero(a.coeffs)
end

## zero and one ##
for T in (:Taylor1, :TaylorN)
    @eval zero(a::$T) = $T(zero.(a.coeffs))
    @eval function one(a::$T)
        b = zero(a)
        b[0] = one(b[0])
        return b
    end
end

function zero(a::HomogeneousPolynomial{T}) where {T<:Number}
    v = zero.(a.coeffs)
    return HomogeneousPolynomial(v, a.order)
end

function zeros(a::HomogeneousPolynomial{T}, order::Int) where {T<:Number}
    order == 0 && return [HomogeneousPolynomial([zero(a[1])], 0)]
    v = Array{HomogeneousPolynomial{T}}(undef, order+1)
    @simd for ord in eachindex(v)
        @inbounds v[ord] = HomogeneousPolynomial(zero(a[1]), ord-1)
    end
    return v
end

zeros(::Type{HomogeneousPolynomial{T}}, order::Int) where {T<:Number} =
    zeros( HomogeneousPolynomial([zero(T)], 0), order)

function one(a::HomogeneousPolynomial{T}) where {T<:Number}
    v = one.(a.coeffs)
    return HomogeneousPolynomial(v, a.order)
end

function ones(a::HomogeneousPolynomial{T}, order::Int) where {T<:Number}
    order == 0 && return [HomogeneousPolynomial([one(a[1])], 0)]
    v = Array{HomogeneousPolynomial{T}}(undef, order+1)
    @simd for ord in eachindex(v)
        @inbounds num_coeffs = size_table[ord]
        @inbounds v[ord] = HomogeneousPolynomial(ones(T, num_coeffs), ord-1)
    end
    return v
end

ones(::Type{HomogeneousPolynomial{T}}, order::Int) where {T<:Number} =
    ones( HomogeneousPolynomial([one(T)], 0), order)



## Addition and substraction ##
for (f, fc) in ((:+, :(add!)), (:-, :(subst!)))

    for T in (:Taylor1, :TaylorN)
        @eval begin
            ($f)(a::$T{T}, b::$T{S}) where {T<:Number,S<:Number} =
                $f(promote(a,b)...)

            function $f(a::$T{T}, b::$T{T}) where {T<:Number}
                if a.order != b.order
                    a, b = fixorder(a, b)
                end
                v = similar(a.coeffs)
                @__dot__ v = $f(a.coeffs, b.coeffs)
                return $T(v, a.order)
            end

            function $f(a::$T)
                v = similar(a.coeffs)
                @__dot__ v = $f(a.coeffs)
                return $T(v, a.order)
            end

            ($f)(a::$T{T}, b::S) where {T<:Number,S<:Number} =
                $f(promote(a,b)...)

            function $f(a::$T{T}, b::T) where {T<:Number}
                coeffs = copy(a.coeffs)
                @inbounds coeffs[1] = $f(a[0], b)
                return $T(coeffs, a.order)
            end

            ($f)(b::S, a::$T{T}) where {T<:Number,S<:Number} =
                $f(promote(b,a)...)

            function $f(b::T, a::$T{T}) where {T<:Number}
                coeffs = similar(a.coeffs)
                @__dot__ coeffs = ($f)(a.coeffs)
                @inbounds coeffs[1] = $f(b, a[0])
                return $T(coeffs, a.order)
            end

            ## add! and subst! ##
            function ($fc)(v::$T{T}, a::$T{T}, k::Int) where {T<:Number}
                @inbounds v[k] = ($f)(a[k])
                return nothing
            end
            function ($fc)(v::$T{T}, a::T, k::Int) where {T<:Number}
                @inbounds v[k] = k==0 ? ($f)(a) : zero(a)
                return nothing
            end
            function ($fc)(v::$T, a::$T, b::$T, k::Int)
                @inbounds v[k] = ($f)(a[k], b[k])
                return nothing
            end
            function ($fc)(v::$T, a::$T, b::Number, k::Int)
                @inbounds v[k] = k==0 ? ($f)(a[0], b) : ($f)(a[k], zero(b))
                return nothing
            end
            function ($fc)(v::$T, a::Number, b::$T, k::Int)
                @inbounds v[k] = k==0 ? ($f)(a, b[0]) : ($f)(zero(a), b[k])
                return nothing
            end
        end
    end

    @eval begin
        ($f)(a::HomogeneousPolynomial{T}, b::HomogeneousPolynomial{S}) where
            {T<:NumberNotSeriesN,S<:NumberNotSeriesN} = $f(promote(a,b)...)

        function $f(a::HomogeneousPolynomial{T}, b::HomogeneousPolynomial{T}) where
                {T<:NumberNotSeriesN}
            @assert a.order == b.order
            v = similar(a.coeffs)
            @__dot__ v = $f(a.coeffs, b.coeffs)
            return HomogeneousPolynomial(v, a.order)
        end

        function $f(a::HomogeneousPolynomial)
            v = similar(a.coeffs)
            @__dot__ v = $f(a.coeffs)
            return HomogeneousPolynomial(v, a.order)
        end

        function ($f)(a::TaylorN{Taylor1{T}}, b::Taylor1{S}) where
                {T<:NumberNotSeries,S<:NumberNotSeries}
            @inbounds aux = $f(a[0][1], b)
            R = eltype(aux)
            coeffs = Array{HomogeneousPolynomial{Taylor1{R}}}(undef, a.order+1)
            coeffs .= a.coeffs
            @inbounds coeffs[1] = aux
            return TaylorN(coeffs, a.order)
        end

        function ($f)(b::Taylor1{S}, a::TaylorN{Taylor1{T}}) where
                {T<:NumberNotSeries,S<:NumberNotSeries}
            @inbounds aux = $f(b, a[0][1])
            R = eltype(aux)
            coeffs = Array{HomogeneousPolynomial{Taylor1{R}}}(undef, a.order+1)
            @__dot__ coeffs = $f(a.coeffs)
            @inbounds coeffs[1] = aux
            return TaylorN(coeffs, a.order)
        end

        function ($f)(a::Taylor1{TaylorN{T}}, b::TaylorN{S}) where
                {T<:NumberNotSeries,S<:NumberNotSeries}
            @inbounds aux = $f(a[0], b)
            R = eltype(aux)
            coeffs = Array{TaylorN{R}}(undef, a.order+1)
            coeffs .= a.coeffs
            @inbounds coeffs[1] = aux
            return Taylor1(coeffs, a.order)
        end

        function ($f)(b::TaylorN{S}, a::Taylor1{TaylorN{T}}) where
                {T<:NumberNotSeries,S<:NumberNotSeries}
            @inbounds aux = $f(b, a[0])
            R = eltype(aux)
            coeffs = Array{TaylorN{R}}(undef, a.order+1)
            @__dot__ coeffs = $f(a.coeffs)
            @inbounds coeffs[1] = aux
            return Taylor1(coeffs, a.order)
        end
    end
end

+(a::Taylor1{T}, b::TaylorN{S}) where {T<:NumberNotSeries,S<:NumberNotSeries} =
    +(promote(a,b)...)
-(a::Taylor1{T}, b::TaylorN{S}) where {T<:NumberNotSeries,S<:NumberNotSeries} =
    -(promote(a,b)...)


## Multiplication ##
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)

    @eval begin
        function *(a::T, b::$T{S}) where {T<:NumberNotSeries,S<:NumberNotSeries}
            @inbounds aux = a * b.coeffs[1]
            v = Array{typeof(aux)}(undef, length(b.coeffs))
            @__dot__ v = a * b.coeffs
            return $T(v, b.order)
        end

        *(b::$T{S}, a::T) where {T<:NumberNotSeries,S<:NumberNotSeries} = a * b

        function *(a::T, b::$T{T}) where {T<:Number}
            v = Array{T}(undef, length(b.coeffs))
            @__dot__ v = a * b.coeffs
            return $T(v, b.order)
        end

        *(b::$T{T}, a::T) where {T<:Number} = a * b
    end
end

for T in (:HomogeneousPolynomial, :TaylorN)

    @eval begin
        function *(a::Taylor1{T}, b::$T{Taylor1{S}}) where
                {T<:NumberNotSeries,S<:NumberNotSeries}

            @inbounds aux = a * b.coeffs[1]
            R = typeof(aux)
            coeffs = Array{R}(undef, length(b.coeffs))
            @__dot__ coeffs = a * b.coeffs
            return $T(coeffs, b.order)
        end

        *(b::$T{Taylor1{R}}, a::Taylor1{T}) where
            {T<:NumberNotSeries,R<:NumberNotSeries} = a * b

        function *(a::$T{T}, b::Taylor1{$T{S}}) where {T<:NumberNotSeries,S<:NumberNotSeries}
            @inbounds aux = a * b[0]
            R = typeof(aux)
            coeffs = Array{R}(undef, length(b.coeffs))
            @__dot__ coeffs = a * b.coeffs
            return Taylor1(coeffs, b.order)
        end

        *(b::Taylor1{$T{S}}, a::$T{T}) where {T<:NumberNotSeries,S<:NumberNotSeries} = a * b
    end
end

for (T, W) in ((:Taylor1, :Number), (:TaylorN, :NumberNotSeriesN))
    @eval function *(a::$T{T}, b::$T{T}) where {T<:$W}
        if a.order != b.order
            a, b = fixorder(a, b)
        end
        c = $T(zero(a[0]), a.order)
        for ord in eachindex(c)
            mul!(c, a, b, ord) # updates c[ord]
        end
        return c
    end
end


*(a::HomogeneousPolynomial{T}, b::HomogeneousPolynomial{S}) where
    {T<:NumberNotSeriesN,S<:NumberNotSeriesN} = *(promote(a,b)...)

function *(a::HomogeneousPolynomial{T}, b::HomogeneousPolynomial{T}) where
        {T<:NumberNotSeriesN}

    order = a.order + b.order

    order > get_order() && return HomogeneousPolynomial(zero(a[1]), get_order())

    res = HomogeneousPolynomial(zero(a[1]), order)
    mul!(res, a, b)
    return res
end


# Internal multiplication functions
for T in (:Taylor1, :TaylorN)
    @eval @inline function mul!(c::$T{T}, a::$T{T}, b::$T{T}, k::Int) where {T<:Number}

        if $T == Taylor1
            @inbounds c[k] = a[0] * b[k]
        else
            @inbounds mul!(c[k], a[0], b[k])
        end
        @inbounds for i = 1:k
            if $T == Taylor1
                c[k] += a[i] * b[k-i]
            else
                mul!(c[k], a[i], b[k-i])
            end
        end

        return nothing
    end

    @eval @inline function mul!(v::$T, a::$T, b::NumberNotSeries, k::Int)
        @inbounds v[k] = a[k] * b
        return nothing
    end
    @eval @inline function mul!(v::$T, a::NumberNotSeries, b::$T, k::Int)
        @inbounds v[k] = a * b[k]
        return nothing
    end
end

@doc doc"""
    mul!(c, a, b, k::Int) --> nothing

Update the `k`-th expansion coefficient `c[k]` of `c = a * b`,
where all `c`, `a`, and `b` are either `Taylor1` or `TaylorN`.

The coefficients are given by

```math
c_k = \sum_{j=0}^k a_j b_{k-j}.
```

""" mul!


"""
    mul!(c, a, b) --> nothing

Return `c = a*b` with no allocation; all arguments are `HomogeneousPolynomial`.

"""
@inline function mul!(c::HomogeneousPolynomial, a::HomogeneousPolynomial,
        b::HomogeneousPolynomial)

    (iszero(b) || iszero(a)) && return nothing

    @inbounds num_coeffs_a = size_table[a.order+1]
    @inbounds num_coeffs_b = size_table[b.order+1]

    @inbounds posTb = pos_table[c.order+1]

    @inbounds indTa = index_table[a.order+1]
    @inbounds indTb = index_table[b.order+1]

    @inbounds for na in 1:num_coeffs_a
        ca = a[na]
        # iszero(ca) && continue
        inda = indTa[na]

        @inbounds for nb in 1:num_coeffs_b
            cb = b[nb]
            # iszero(cb) && continue
            indb = indTb[nb]

            pos = posTb[inda + indb]
            c[pos] += ca * cb
        end
    end

    return nothing
end



## Division ##
function /(a::Taylor1{Rational{T}}, b::S) where {T<:Integer,S<:NumberNotSeries}
    R = typeof( a[0] // b)
    v = Array{R}(undef, a.order+1)
    @__dot__ v = a.coeffs // b
    return Taylor1(v, a.order)
end

for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval function /(a::$T{T}, b::S) where {T<:NumberNotSeries,S<:NumberNotSeries}
        @inbounds aux = a.coeffs[1] / b
        v = Array{typeof(aux)}(undef, length(a.coeffs))
        @__dot__ v = a.coeffs / b
        return $T(v, a.order)
    end

    @eval function /(a::$T{T}, b::T) where {T<:Number}
        @inbounds aux = a.coeffs[1] / b
        v = Array{typeof(aux)}(undef, length(a.coeffs))
        @__dot__ v = a.coeffs / b
        return $T(v, a.order)
    end
end

for T in (:HomogeneousPolynomial, :TaylorN)
    @eval function /(b::$T{Taylor1{S}}, a::Taylor1{T}) where
            {T<:NumberNotSeries,S<:NumberNotSeries}
        @inbounds aux = b.coeffs[1] / a
        R = typeof(aux)
        coeffs = Array{R}(undef, length(b.coeffs))
        @__dot__ coeffs = b.coeffs / a
        return $T(coeffs, b.order)
    end

    @eval function /(b::Taylor1{$T{S}}, a::$T{T}) where
            {T<:NumberNotSeries,S<:NumberNotSeries}
        @inbounds aux = b[0] / a
        R = typeof(aux)
        coeffs = Array{R}(undef, length(b.coeffs))
        @__dot__ coeffs = b.coeffs / a
        return Taylor1(coeffs, b.order)
    end
end

/(a::Taylor1{T}, b::Taylor1{S}) where {T<:Number,S<:Number} = /(promote(a,b)...)

function /(a::Taylor1{T}, b::Taylor1{T}) where {T<:Number}
    iszero(a) && !iszero(b) && return zero(a)
    if a.order != b.order
        a, b = fixorder(a, b)
    end

    # order and coefficient of first factorized term
    ordfact, cdivfact = divfactorization(a, b)

    c = Taylor1(cdivfact, a.order-ordfact)
    for ord in eachindex(c)
        div!(c, a, b, ord) # updates c[ord]
    end

    return c
end

/(a::TaylorN{T}, b::TaylorN{S}) where
    {T<:NumberNotSeriesN,S<:NumberNotSeriesN} = /(promote(a,b)...)

function /(a::TaylorN{T}, b::TaylorN{T}) where {T<:NumberNotSeriesN}
    @assert !iszero(constant_term(b))

    if a.order != b.order
        a, b = fixorder(a, b)
    end

    # first coefficient
    @inbounds cdivfact = a[0] / constant_term(b)
    c = TaylorN(cdivfact, a.order)
    for ord in eachindex(c)
        div!(c, a, b, ord) # updates c[ord]
    end

    return c
end


@inline function divfactorization(a1::Taylor1, b1::Taylor1)
    # order of first factorized term; a1 and b1 assumed to be of the same order
    a1nz = findfirst(a1)
    b1nz = findfirst(b1)
    a1nz = a1nz ≥ 0 ? a1nz : a1.order
    b1nz = b1nz ≥ 0 ? b1nz : a1.order
    ordfact = min(a1nz, b1nz)
    cdivfact = a1[ordfact] / b1[ordfact]

    # Is the polynomial factorizable?
    iszero(b1[ordfact]) && throw( ArgumentError(
        """Division does not define a Taylor1 polynomial;
        order k=$(ordfact) => coeff[$(ordfact)]=$(cdivfact).""") )

    return ordfact, cdivfact
end


## TODO: Implement factorization (divfactorization) for TaylorN polynomials


# Homogeneous coefficient for the division
@doc doc"""
    div!(c, a, b, k::Int)

Compute the `k-th` expansion coefficient `c[k]` of `c = a / b`,
where all `c`, `a` and `b` are either `Taylor1` or `TaylorN`.

The coefficients are given by

```math
c_k =  \frac{1}{b_0} \big(a_k - \sum_{j=0}^{k-1} c_j b_{k-j}\big).
```

For `Taylor1` polynomials, a similar formula is implemented which
exploits `k_0`, the order of the first non-zero coefficient of `a`.
""" div!

@inline function div!(c::Taylor1, a::Taylor1, b::Taylor1, k::Int)

    # order and coefficient of first factorized term
    ordfact, cdivfact = divfactorization(a, b)
    if k == 0
        @inbounds c[0] = cdivfact
        return nothing
    end

    imin = max(0, k+ordfact-b.order)
    @inbounds c[k] = c[imin] * b[k+ordfact-imin]
    @inbounds for i = imin+1:k-1
        c[k] += c[i] * b[k+ordfact-i]
    end
    if k+ordfact ≤ b.order
        @inbounds c[k] = (a[k+ordfact]-c[k]) / b[ordfact]
    else
        @inbounds c[k] = - c[k] / b[ordfact]
    end
    return nothing
end

@inline function div!(v::Taylor1, a::Taylor1, b::NumberNotSeries, k::Int)
    @inbounds v[k] = a[k] / b
    return nothing
end

div!(v::Taylor1, b::NumberNotSeries, a::Taylor1, k::Int) =
    div!(v::Taylor1, Taylor1(b, a.order), a, k)

@inline function div!(c::TaylorN, a::TaylorN, b::TaylorN, k::Int)
    if k==0
        @inbounds c[0] = a[0] / constant_term(b)
        return nothing
    end

    @inbounds for i = 0:k-1
        mul!(c[k], c[i], b[k-i])
    end
    @inbounds c[k] = (a[k] - c[k]) / constant_term(b)
    return nothing
end





"""
    mul!(Y, A, B)

Multiply A*B and save the result in Y.
"""
function mul!(y::Vector{Taylor1{T}},
        a::Union{Matrix{T},SparseMatrixCSC{T}},
        b::Vector{Taylor1{T}}) where {T<:Number}

    n, k = size(a)
    @assert (length(y)== n && length(b)== k)

    # determine the maximal order of b
    # order = maximum([b1.order for b1 in b])
    order = maximum(get_order.(b))

    # Use matrices of coefficients (of proper size) and mul!
    # B = zeros(T, k, order+1)
    B = Array{T}(undef, k, order+1)
    B = zero.(B)
    for i = 1:k
        @inbounds ord = b[i].order
        @inbounds for j = 1:ord+1
            B[i,j] = b[i][j-1]
        end
    end
    Y = Array{T}(undef, n, order+1)
    mul!(Y, a, B)
    @inbounds for i = 1:n
        # y[i] = Taylor1( collect(Y[i,:]), order)
        y[i] = Taylor1( Y[i,:], order)
    end

    return y
end


# Adapted from (Julia v1.2) stdlib/v1.2/LinearAlgebra/src/dense.jl#721-734,
# licensed under MIT "Expat".
# Specialize a method of `inv` for Matrix{Taylor1{T}}. Simply, avoid pivoting,
# since the polynomial field is not an ordered one.
# function Base.inv(A::StridedMatrix{Taylor1{T}}) where T
#     checksquare(A)
#     S = Taylor1{typeof((one(T)*zero(T) + one(T)*zero(T))/one(T))}
#     AA = convert(AbstractArray{S}, A)
#     if istriu(AA)
#         Ai = triu!(parent(inv(UpperTriangular(AA))))
#     elseif istril(AA)
#         Ai = tril!(parent(inv(LowerTriangular(AA))))
#     else
#         # Do not use pivoting !!
#         Ai = inv!(lu(AA, Val(false)))
#         Ai = convert(typeof(parent(Ai)), Ai)
#     end
#     return Ai
# end

# see https://github.com/JuliaLang/julia/pull/40623
const LU_RowMaximum = VERSION >= v"1.7.0-DEV.1188" ? RowMaximum() : Val(true)
const LU_NoPivot = VERSION >= v"1.7.0-DEV.1188" ? NoPivot() : Val(false)

# Adapted from (Julia v1.2) stdlib/v1.2/LinearAlgebra/src/lu.jl#240-253
# and (Julia v1.4.0-dev) stdlib/LinearAlgebra/v1.4/src/lu.jl#270-274,
# licensed under MIT "Expat".
# Specialize a method of `lu` for Matrix{Taylor1{T}}, which avoids pivoting,
# since the polynomial field is not an ordered one.
# We can't assume an ordered field so we first try without pivoting
function lu(A::AbstractMatrix{Taylor1{T}}; check::Bool = true) where {T<:Number}
    S = Taylor1{lutype(T)}
    F = lu!(copy_oftype(A, S), LU_NoPivot; check = false)
    if issuccess(F)
        return F
    else
        return lu!(copy_oftype(A, S), LU_RowMaximum; check = check)
    end
end
