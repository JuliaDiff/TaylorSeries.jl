# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

## Auxiliary function ##

"""
    _coeffsHP(x::T, order::Int) where {T<:Number}
    _coeffsHP(coeffs::DenseVector{T}, order::Int) where {T<:Number}

Returns a `Memory{T}` of size `size_table[order+1]`
to be used in the construction of an HomogeneousPolynomial of order
`order`. The returned vector has the first entries of `coeffs`,
and then it's filled with zeros.
"""
function _coeffsHP(x::T, order::Int) where {T<:NumberNotSeries}
    @assert order ≤ get_order()
    num_coeffs = size_table[order+1]
    v = Memory{T}(undef, num_coeffs)
    v .= zero.(x)
    v[1] = x
    return v
end
function _coeffsHP(x::Taylor1{T}, order::Int) where {T<:NumberNotSeries}
    @assert order ≤ get_order()
    v = Memory{Taylor1{T}}(undef, size_table[order+1])
    v .= zero.(x)
    v[1].coeffs .= x.coeffs
    return v
end
function _coeffsHP(coeffs::DenseVector{T}, order::Int) where {T<:Number}
    @assert order ≤ get_order()
    ll = length( coeffs )
    num_coeffs = size_table[order+1]
    # @assert ll ≤ num_coeffs
    v = Memory{T}(undef, num_coeffs)
    for ord in eachindex(coeffs)
        v[ord] = coeffs[ord]
    end
    for ord in ll+1:num_coeffs
        v[ord] = zero(v[1])
    end
    return v
end

"""
    _coeffsTN(v::DenseVector{T}, order::Int) where {T<:Number}

Returns a `Memory{HomogeneousPolynomial{T}}` of
size `order+1`, to be used in the construction of a TaylorN{T} of order
`order`. The returned vector has the entries of `v` at the proper
location according to their `order`, and otherwise it is filled with
the corresponding zeros.
"""
function _coeffsTN(v::DenseVector{HomogeneousPolynomial{T}},
        order::Int) where {T}
    coeffs = zeros(v[1], order)
    vord = get_order.(v)
    max_order = maximum(vord)
    if allunique(vord) && (max_order ≤ order)
        for i in eachindex(v)
            coeffs[vord[i]+1].coeffs .= v[i].coeffs
        end
    elseif max_order ≤ order
        for i in eachindex(v)
            coeffs[vord[i]+1].coeffs .+= v[i].coeffs
        end
    else
        for i in eachindex(v)
            ord = vord[i]
            ord > order && continue
            coeffs[ord+1].coeffs .+= v[i].coeffs
        end
    end
    return coeffs
end


## Minimum order of an HomogeneousPolynomial compatible with the vector's length
function orderH(coeffs::DenseVector{T}) where {T<:Number}
    ord = 0
    ll = length(coeffs)
    for i = 1:get_order()+1
        ll ≤ size_table[i] && return ord
        ord += 1
    end
    return ord
end

## Maximum order of a HomogeneousPolynomial vector; used by TaylorN constructor
maxorderH(v::DenseVector{HomogeneousPolynomial{T}}) where {T<:Number} =
    isempty(v) ? 0 : maximum(get_order.(v))


## getcoeff ##
"""
    getcoeff(a, n)

Return the coefficient of order `n::Int` of a `a::Taylor1` polynomial; the constant
term corresponds to n=0.
"""
getcoeff(a::Taylor1, n::Int) = (@assert 0 ≤ n ≤ get_order(a); return a[n])

getindex(a::Taylor1, n::Int) = a.coeffs[n+1]
getindex(a::Taylor1, u::UnitRange{Int}) = view(a.coeffs, u .+ 1 )
getindex(a::Taylor1, c::Colon) = view(a.coeffs, c)
getindex(a::Taylor1{T}, u::StepRange{Int,Int}) where {T<:Number} =
    view(a.coeffs, u .+ 1)

setindex!(a::Taylor1{T}, x::T, n::Int) where {T<:NumberNotSeries} =
    a.coeffs[n+1] = x
# setindex!(a::Taylor1{T}, x::T, n::Int) where {T<:AbstractSeries} =
#     setindex!(a.coeffs, deepcopy(x), n+1)
setindex!(a::Taylor1{TaylorN{T}}, x::TaylorN{T}, n::Int) where
    {T<:NumberNotSeries} = a.coeffs[n+1] = TaylorN(x.coeffs, get_order(x))
setindex!(a::TaylorN{Taylor1{T}}, x::Taylor1{T}, n::Int) where
    {T<:NumberNotSeries} = a.coeffs[n+1] = Taylor1{T}(x.coeffs[:])
function setindex!(a::Taylor1{Taylor1{T}}, x::Taylor1{T}, n::Int) where
        {T<:Taylor1{<:Number}}
    a.coeffs[n+1] = zero(x)
    for i in eachindex(x)
        a.coeffs[n+1].coeffs[i+1] = x.coeffs[i+1]
    end
    return a.coeffs[n+1]
end
setindex!(a::Taylor1{Taylor1{T}}, x::Taylor1{T}, n::Int) where
    {T<:NumberNotSeries} = a.coeffs[n+1] = Taylor1(x.coeffs[:], get_order(x))
setindex!(a::Taylor1{T}, x::T, u::UnitRange{Int}) where {T<:Number} =
    a.coeffs[u .+ 1] .= x
function setindex!(a::Taylor1{T}, x::DenseVector{T},
        u::UnitRange{Int}) where {T<:Number}
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a.coeffs[u[ind]+1] = x[ind]
    end
end
setindex!(a::Taylor1{T}, x::T, c::Colon) where {T<:Number} = a.coeffs[c] .= x
setindex!(a::Taylor1{T}, x::DenseVector{T}, c::Colon) where {T<:Number} =
    a.coeffs[c] .= x
setindex!(a::Taylor1{T}, x::T, u::StepRange{Int,Int}) where {T<:Number} =
    a.coeffs[u[:] .+ 1] .= x
function setindex!(a::Taylor1{T}, x::Vector{T}, u::StepRange{Int,Int}) where
        {T<:Number}
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a.coeffs[u[ind]+1] = x[ind]
    end
end


"""
    getcoeff(a, v)

Return the coefficient of `a::HomogeneousPolynomial`, specified by `v`,
which is a tuple (or vector) with the indices of the specific
monomial.
"""
function getcoeff(a::HomogeneousPolynomial, v::NTuple{N,Int}) where {N}
    @assert N == get_numvars() && all(v .>= 0)
    kdic = in_base(get_order(), v)
    @inbounds n = pos_table[get_order(a)+1][kdic]
    a[n]
end
getcoeff(a::HomogeneousPolynomial, v::DenseVector{Int}) =
    getcoeff(a, (v...,))

getindex(a::HomogeneousPolynomial, n::Int) = a.coeffs[n]
getindex(a::HomogeneousPolynomial, n::UnitRange{Int}) = view(a.coeffs, n)
getindex(a::HomogeneousPolynomial, c::Colon) = view(a.coeffs, c)
getindex(a::HomogeneousPolynomial, u::StepRange{Int,Int}) = view(a.coeffs, u[:])

setindex!(a::HomogeneousPolynomial{T}, x::T, n::Int) where {T<:Number} =
    a.coeffs[n] = x
setindex!(a::HomogeneousPolynomial{T}, x::T, n::UnitRange{Int}) where
    {T<:Number} = a.coeffs[n] .= x
setindex!(a::HomogeneousPolynomial{T}, x::DenseVector{T},
    n::UnitRange{Int}) where {T<:Number} = a.coeffs[n] .= x
setindex!(a::HomogeneousPolynomial{T}, x::T, c::Colon) where {T<:Number} =
    a.coeffs[c] .= x
setindex!(a::HomogeneousPolynomial{T}, x::DenseVector{T},
    c::Colon) where {T<:Number} = a.coeffs[c] .= x
setindex!(a::HomogeneousPolynomial{T}, x::T,
    u::StepRange{Int,Int}) where {T<:Number} = a.coeffs[u[:]] .= x
setindex!(a::HomogeneousPolynomial{T}, x::DenseVector{T},
    u::StepRange{Int,Int}) where {T<:Number} = a.coeffs[u[:]] .= x[:]


"""
    getcoeff(a, v)

Return the coefficient of `a::TaylorN`, specified by `v`,
which is a tuple (or vector) with the indices of the specific
monomial.
"""
function getcoeff(a::TaylorN, v::NTuple{N,Int}) where {N}
    order = sum(v)
    @assert order ≤ get_order(a)
    getcoeff(a[order], v)
end
getcoeff(a::TaylorN, v::DenseVector{Int}) = getcoeff(a, (v...,))

getindex(a::TaylorN, n::Int) = a.coeffs[n+1]
getindex(a::TaylorN, u::UnitRange{Int}) = view(a.coeffs, u .+ 1)
getindex(a::TaylorN, c::Colon) = view(a.coeffs, c)
getindex(a::TaylorN, u::StepRange{Int,Int}) = view(a.coeffs, u[:] .+ 1)

function setindex!(a::TaylorN{T}, x::HomogeneousPolynomial{T}, n::Int) where
        {T<:Number}
    @assert get_order(x) == n
    return a.coeffs[n+1] = x
end
setindex!(a::TaylorN{T}, x::T, n::Int) where {T<:Number} =
    a.coeffs[n+1] = HomogeneousPolynomial(x, n)
function setindex!(a::TaylorN{T}, x::T, u::UnitRange{Int}) where {T<:Number}
    for ind in u
        a[ind] = x
    end
    return a[u]
end
function setindex!(a::TaylorN{T},
        x::DenseVector{HomogeneousPolynomial{T}},
        u::UnitRange{Int}) where {T<:Number}
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a[u[ind]] = x[ind]
    end
    return a[u]
end
function setindex!(a::TaylorN{T}, x::DenseVector{T},
        u::UnitRange{Int}) where {T<:Number}
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a[u[ind]] = x[ind]
    end
    return a[u]
end
setindex!(a::TaylorN{T}, x::T, ::Colon) where {T<:Number} =
    (a[0:end] = x; a[:])
setindex!(a::TaylorN{T}, x::DenseVector{HomogeneousPolynomial{T}},
    ::Colon) where {T<:Number} = (a[0:end] = x; a[:])
setindex!(a::TaylorN{T}, x::DenseVector{T}, ::Colon) where {T<:Number} =
    (a[0:end] = x; a[:])
function setindex!(a::TaylorN{T}, x::T, u::StepRange{Int,Int}) where {T<:Number}
    for ind in u
        a[ind] = x
    end
    return a[u]
end
function setindex!(a::TaylorN{T}, x::DenseVector{HomogeneousPolynomial{T}},
        u::StepRange{Int,Int}) where {T<:Number}
    # a[u[:]] .= x[:]
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a[u[ind]] = x[ind]
    end
    return a[u]
end
function setindex!(a::TaylorN{T}, x::Vector{T},
        u::StepRange{Int,Int}) where {T<:Number}
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a[u[ind]] = x[ind]
    end
    return a[u]
end


## eltype, length, get_order, etc ##
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval begin
        if $T == HomogeneousPolynomial
            @inline get_order(a::$T) = a.order
            @inline iterate(a::$T, state=1) =
                state > length(a) ? nothing : (a.coeffs[state], state+1)
            # Base.iterate(rS::Iterators.Reverse{$T}, state=rS.itr.order) = state < 0 ? nothing : (a.coeffs[state], state-1)
            @inline length(a::$T) = size_table[get_order(a)+1]
            @inline firstindex(a::$T) = 1
            @inline lastindex(a::$T) = length(a)
        else
            @inline get_order(a::$T) = size(a.coeffs, 1)-1
            @inline iterate(a::$T, state=0) =
                state > get_order(a) ? nothing : (a.coeffs[state+1], state+1)
            # Base.iterate(rS::Iterators.Reverse{$T}, state=rS.itr.order) = state < 0 ? nothing : (a.coeffs[state], state-1)
            @inline length(a::$T) = length(a.coeffs)
            @inline firstindex(a::$T) = 0
            @inline lastindex(a::$T) = get_order(a)
        end
        @inline eachindex(a::$T) = firstindex(a):lastindex(a)
        @inline numtype(::$T{S}) where {S<:Number} = S
        @inline size(a::$T) = size(a.coeffs)
        @inline axes(a::$T) = ()
    end
end
numtype(a) = eltype(a)

@doc doc"""
    numtype(a::AbstractSeries)

Returns the type of the elements of the coefficients of `a`.
""" numtype

# Dumb methods included to properly export normalize_taylor
# (if IntervalArithmetic is loaded)
@inline normalize_taylor(a::AbstractSeries) = a
@inline aff_normalize(a::AbstractSeries) = a


## _minorder
function _minorder(a, b)
    minorder, maxorder = minmax(get_order(a), get_order(b))
    if minorder ≤ 0
        minorder = maxorder
    end
    return minorder
end


## fixorder ##
for T in (:Taylor1, :TaylorN)
    @eval begin
        @inline function fixorder(a::$T, b::$T)
            get_order(a) == get_order(b) && return a, b
            minorder = _minorder(a, b)
            return $T(a.coeffs, minorder), $T(b.coeffs, minorder)
        end
    end
end

function fixorder(a::HomogeneousPolynomial, b::HomogeneousPolynomial)
    @assert get_order(a) == get_order(b)
    return a, b
end

for T in (:HomogeneousPolynomial, :TaylorN)
    @eval function fixorder(a::Taylor1{$T{T}}, b::Taylor1{$T{S}}) where
            {T<:NumberNotSeries, S<:NumberNotSeries}
        (get_order(a) == get_order(b)) && (all(get_order.(a.coeffs) .==
            get_order.(b.coeffs))) && return a, b
        minordT = _minorder(a, b)
        aa = Taylor1(a.coeffs, minordT)
        bb = Taylor1(b.coeffs, minordT)
        for ind in eachindex(aa)
            get_order(aa[ind]) == get_order(bb[ind]) && continue
            minordQ = _minorder(aa[ind], bb[ind])
            aa[ind] = $T(aa[ind].coeffs, minordQ)
            bb[ind] = $T(bb[ind].coeffs, minordQ)
        end
        return aa, bb
    end
end


## minlength
@inline function minlength(a::Taylor1, b::Taylor1)
    length(eachindex(a)) < length(eachindex(b)) && return a
    return b
end
@inline function minlength(a::Taylor1, b::Taylor1, c::Taylor1)
    length(minlength(a, c)) < length(eachindex(b)) && return minlength(a, c)
    return b
end


## _isthinzero
"""
    _isthinzero(x)

Generic wrapper to function `iszero`, which allows using the correct
function for `Interval`s
"""
_isthinzero(x) = iszero(x)


## findfirst, findlast
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    # Finds the first non zero entry
    @eval function Base.findfirst(a::$T{T}) where {T<:Number}
        first = findfirst(!_isthinzero, a.coeffs)
        isnothing(first) && return -1
        if $T == HomogeneousPolynomial
            return first
        else
            return first-1
        end
    end
    # Finds the last non-zero entry
    @eval function Base.findlast(a::$T{T}) where {T<:Number}
        last = findlast(!_isthinzero, a.coeffs)
        isnothing(last) && return -1
        if $T == HomogeneousPolynomial
            return last
        else
            return last-1
        end
    end
end


## copyto! ##
# Inspired from base/abstractarray.jl, line 665
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval function copyto!(dst::$T{T}, src::$T{T}) where {T<:Number}
        length(dst) < length(src) && throw(
            ArgumentError(string(
                "Destination has fewer elements than required; no copy performed")))
        destiter = eachindex(dst)
        y = iterate(destiter)
        for x in src
            dst[y[1]] = x
            y = iterate(destiter, y[2])
        end
        return dst
    end
end


"""
    constant_term(a)

Return the constant value (zero order coefficient) for `Taylor1`
and `TaylorN`. The fallback behavior is to return `a` itself if
`a::Number`, or `a[1]` when `a::Vector`.
"""
constant_term(a::Taylor1) = a[0]

constant_term(a::TaylorN) = a[0][1]

constant_term(a::DenseVector{T}) where {T<:Number} = constant_term.(a)

constant_term(a::Number) = a

"""
    linear_polynomial(a)

Returns the linear part of `a` as a polynomial (`Taylor1` or `TaylorN`),
*without* the constant term. The fallback behavior is to return `a` itself.
"""
linear_polynomial(a::Taylor1) = Taylor1([zero(a[1]), a[1]], get_order(a))

linear_polynomial(a::HomogeneousPolynomial) = HomogeneousPolynomial(a[1],
    get_order(a))

linear_polynomial(a::TaylorN) = TaylorN(a[1], get_order(a))

linear_polynomial(a::DenseVector{T}) where {T<:Number} = linear_polynomial.(a)

linear_polynomial(a::Number) = a

"""
    nonlinear_polynomial(a)

Returns the nonlinear part of `a`. The fallback behavior is to return `zero(a)`.
"""
nonlinear_polynomial(a::AbstractSeries) = a - constant_term(a) - linear_polynomial(a)

nonlinear_polynomial(a::DenseVector{T}) where {T<:Number} = nonlinear_polynomial.(a)

nonlinear_polynomial(a::Number) = zero(a)


"""
    @isonethread (expr)

Internal macro used to check the number of threads in use, to prevent a data race
that modifies `coeff_table` when using `differentiate` or `integrate`; see
https://github.com/JuliaDiff/TaylorSeries.jl/issues/318.

This macro is inspired by the macro `@threaded`; see
https://github.com/trixi-framework/Trixi.jl/blob/main/src/auxiliary/auxiliary.jl;
and https://github.com/trixi-framework/Trixi.jl/pull/426/files.
"""
macro isonethread(expr)
    return esc(quote
        if Threads.nthreads() == 1
            $(expr)
        else
            copy($(expr))
        end
    end)
end
