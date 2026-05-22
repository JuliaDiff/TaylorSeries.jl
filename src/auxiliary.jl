# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

## Auxiliary function ##

"""
    space(a::Union{HomogeneousPolynomial,TaylorN})

Return the `JetSpace` associated with the multivariate Taylor object `a`.
"""
@inline space(a::HomogeneousPolynomial) = a.space
@inline space(a::TaylorN) = a.space

function _space_mismatch_error(space_a::JetSpace, space_b::JetSpace)
    throw(ArgumentError(
        "JetSpace mismatch: operands belong to different spaces. " *
        "Use an explicit projection or conversion before combining them."))
end

"""
    _check_same_space(space_a::JetSpace, space_b::JetSpace)
    _check_same_space(a::Taylor1, b::Taylor1)
    _check_same_space(a::Taylor1, b::Taylor1, c::Taylor1)
    _check_same_space(a::Union{HomogeneousPolynomial,TaylorN},
        b::Union{HomogeneousPolynomial,TaylorN})
    _check_same_space(a::Union{HomogeneousPolynomial,TaylorN},
        b::Union{HomogeneousPolynomial,TaylorN},
        c::Union{HomogeneousPolynomial,TaylorN})
    _check_same_space(space::JetSpace, v::AbstractVector{<:HomogeneousPolynomial})

Throw an `ArgumentError` unless all arguments belong to the same `JetSpace`
by object identity.
"""
@inline function _check_same_space(space_a::JetSpace, space_b::JetSpace)
    space_a === space_b || _space_mismatch_error(space_a, space_b)
    return nothing
end
@inline _check_same_space(a::Taylor1, b::Taylor1) = nothing
@inline _check_same_space(a::Taylor1, b::Taylor1, c::Taylor1) = nothing
@inline _check_same_space(a::Union{HomogeneousPolynomial,TaylorN},
    b::Union{HomogeneousPolynomial,TaylorN}) =
        _check_same_space(space(a), space(b))
@inline function _check_same_space(a::Union{HomogeneousPolynomial,TaylorN},
        b::Union{HomogeneousPolynomial,TaylorN},
        c::Union{HomogeneousPolynomial,TaylorN})
    _check_same_space(a, b)
    _check_same_space(a, c)
    return nothing
end

function _check_same_space(space::JetSpace,
        v::AbstractVector{<:HomogeneousPolynomial})
    for pol in v
        _check_same_space(space, pol.space)
    end
    return nothing
end

function _space_from_homogeneous_vector(v::AbstractVector{<:HomogeneousPolynomial},
        fallback::JetSpace)
    isempty(v) && return fallback
    space = v[1].space
    _check_same_space(space, v)
    return space
end

_constant_series_like(a::Taylor1, x, order::Int) = Taylor1(x, order)
_constant_series_like(a::TaylorN, x, order::Int) = TaylorN(a.space, x, order)

"""
    _coeffsHP(x::T, order::Int) where {T<:Number}
    _coeffsHP(coeffs::AbstractArray{T,1}, order::Int) where {T<:Number}
    _coeffsHP(space::JetSpace, x::T, order::Int) where {T<:Number}
    _coeffsHP(space::JetSpace, coeffs::AbstractArray{T,1}, order::Int) where {T<:Number}

Returns a `FixedSizeVectorDefault` of size `space.size_table[order+1]`
to be used in the construction of a `HomogeneousPolynomial` of order
`order`. The returned vector has the first entries of `coeffs`,
and then is filled with zeros.
"""
function _coeffsHP(space::JetSpace, x::T, order::Int) where {T<:NumberNotSeries}
    @assert order ≤ TS.order(space)
    num_coeffs = space.size_table[order+1]
    v = FixedSizeVectorDefault{T}(undef, num_coeffs)
    v .= zero.(x)
    v[1] = x
    return v
end
_coeffsHP(x::T, order::Int) where {T<:NumberNotSeries} =
    _coeffsHP(default_space[], x, order)
function _coeffsHP(space::JetSpace, x::Taylor1{T}, order::Int) where
        {T<:NumberNotSeries}
    @assert order ≤ TS.order(space)
    v = FixedSizeVectorDefault{Taylor1{T}}(undef, space.size_table[order+1])
    v .= zero.(x)
    v[1].coeffs .= x.coeffs
    return v
end
_coeffsHP(x::Taylor1{T}, order::Int) where {T<:NumberNotSeries} =
    _coeffsHP(default_space[], x, order)
function _coeffsHP(space::JetSpace, coeffs::AbstractArray{T,1},
        order::Int) where {T<:Number}
    @assert order ≤ TS.order(space)
    ll = length( coeffs )
    num_coeffs = space.size_table[order+1]
    num_coeffs == ll && return FixedSizeVectorDefault(coeffs)
    # @assert ll ≤ num_coeffs
    v = FixedSizeVectorDefault{T}(undef, num_coeffs)
    for ord in eachindex(coeffs)
        v[ord] = coeffs[ord]
    end
    v[ll+1:num_coeffs] .= zero.(v[1])
    return v
end
_coeffsHP(coeffs::AbstractArray{T,1}, order::Int) where {T<:Number} =
    _coeffsHP(default_space[], coeffs, order)

"""
    _coeffsTN(v::AbstractArray{T,1}, order::Int) where {T<:Number}

Returns a `FixedSizeVectorDefault{HomogeneousPolynomial{T}}` of
size `order+1`, to be used in the construction of a TaylorN{T} of order
`order`. The returned vector has the entries of `v` at the proper
location according to their `order`, and otherwise it is filled with
the corresponding zeros.
"""
function _coeffsTN(space::JetSpace, v::AbstractVector{HomogeneousPolynomial{T}},
        order::Int) where {T}
    coeffs = zeros(HomogeneousPolynomial(space, v[1][1], TS.order(v[1])), order)
    vord = TS.order.(v)
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
_coeffsTN(v::AbstractVector{HomogeneousPolynomial{T}}, order::Int) where {T} =
    _coeffsTN(_space_from_homogeneous_vector(v, default_space[]), v, order)


## Minimum order of an HomogeneousPolynomial compatible with the vector's length
function orderH(space::JetSpace, coeffs::AbstractArray{T,1}) where {T<:Number}
    ord = 0
    ll = length(coeffs)
    for i = 1:TS.order(space)+1
        ll ≤ space.size_table[i] && return ord
        ord += 1
    end
    return ord
end
orderH(coeffs::AbstractArray{T,1}) where {T<:Number} =
    orderH(default_space[], coeffs)

## Maximum order of a HomogeneousPolynomial vector; used by TaylorN constructor
maxorderH(v::AbstractArray{HomogeneousPolynomial{T},1}) where {T<:Number} =
    isempty(v) ? 0 : maximum(TS.order.(v))


## getcoeff ##
"""
    getcoeff(a, n)

Return the coefficient of order `n::Int` of a `a::Taylor1` polynomial; the constant
term corresponds to n=0.
"""
getcoeff(a::Taylor1, n::Int) = (@assert 0 ≤ n ≤ TS.order(a); return a[n])

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
    {T<:NumberNotSeries} = a.coeffs[n+1] = TaylorN(x.coeffs, TS.order(x))
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
    {T<:NumberNotSeries} = a.coeffs[n+1] = Taylor1(x.coeffs[:], TS.order(x))
setindex!(a::Taylor1{T}, x::T, u::UnitRange{Int}) where {T<:Number} =
    a.coeffs[u .+ 1] .= x
function setindex!(a::Taylor1{T}, x::AbstractArray{T,1},
        u::UnitRange{Int}) where {T<:Number}
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a.coeffs[u[ind]+1] = x[ind]
    end
end
setindex!(a::Taylor1{T}, x::T, c::Colon) where {T<:Number} = a.coeffs[c] .= x
setindex!(a::Taylor1{T}, x::AbstractArray{T,1}, c::Colon) where {T<:Number} =
    a.coeffs[c] .= x
setindex!(a::Taylor1{T}, x::T, u::StepRange{Int,Int}) where {T<:Number} =
    a.coeffs[u[:] .+ 1] .= x
function setindex!(a::Taylor1{T}, x::Array{T,1}, u::StepRange{Int,Int}) where {T<:Number}
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
    sp = space(a)
    @assert N == get_numvars(sp) && all(v .>= 0)
    kdic = in_base(TS.order(sp), v)
    @inbounds n = sp.pos_table[TS.order(a)+1][kdic]
    a[n]
end
getcoeff(a::HomogeneousPolynomial, v::AbstractArray{Int,1}) =
    getcoeff(a, (v...,))

getindex(a::HomogeneousPolynomial, n::Int) = a.coeffs[n]
getindex(a::HomogeneousPolynomial, n::UnitRange{Int}) = view(a.coeffs, n)
getindex(a::HomogeneousPolynomial, c::Colon) = view(a.coeffs, c)
getindex(a::HomogeneousPolynomial, u::StepRange{Int,Int}) = view(a.coeffs, u[:])

setindex!(a::HomogeneousPolynomial{T}, x::T, n::Int) where {T<:Number} =
    a.coeffs[n] = x
setindex!(a::HomogeneousPolynomial{T}, x::T, n::UnitRange{Int}) where
    {T<:Number} = a.coeffs[n] .= x
setindex!(a::HomogeneousPolynomial{T}, x::AbstractArray{T,1},
    n::UnitRange{Int}) where {T<:Number} = a.coeffs[n] .= x
setindex!(a::HomogeneousPolynomial{T}, x::T, c::Colon) where {T<:Number} =
    a.coeffs[c] .= x
setindex!(a::HomogeneousPolynomial{T}, x::AbstractArray{T,1},
    c::Colon) where {T<:Number} = a.coeffs[c] .= x
setindex!(a::HomogeneousPolynomial{T}, x::T,
    u::StepRange{Int,Int}) where {T<:Number} = a.coeffs[u[:]] .= x
setindex!(a::HomogeneousPolynomial{T}, x::AbstractArray{T,1},
    u::StepRange{Int,Int}) where {T<:Number} = a.coeffs[u[:]] .= x[:]


"""
    getcoeff(a, v)

Return the coefficient of `a::TaylorN`, specified by `v`,
which is a tuple (or vector) with the indices of the specific
monomial.
"""
function getcoeff(a::TaylorN, v::NTuple{N,Int}) where {N}
    order = sum(v)
    @assert order ≤ TS.order(a)
    getcoeff(a[order], v)
end
getcoeff(a::TaylorN, v::AbstractArray{Int,1}) = getcoeff(a, (v...,))

getindex(a::TaylorN, n::Int) = a.coeffs[n+1]
getindex(a::TaylorN, u::UnitRange{Int}) = view(a.coeffs, u .+ 1)
getindex(a::TaylorN, c::Colon) = view(a.coeffs, c)
getindex(a::TaylorN, u::StepRange{Int,Int}) = view(a.coeffs, u[:] .+ 1)

function setindex!(a::TaylorN{T}, x::HomogeneousPolynomial{T}, n::Int) where
        {T<:Number}
    @assert TS.order(x) == n
    _check_same_space(a, x)
    return a.coeffs[n+1] = x
end
setindex!(a::TaylorN{T}, x::T, n::Int) where {T<:Number} =
    a.coeffs[n+1] = HomogeneousPolynomial(a.space, x, n)
function setindex!(a::TaylorN{T}, x::T, u::UnitRange{Int}) where {T<:Number}
    for ind in u
        a[ind] = x
    end
    return a[u]
end
function setindex!(a::TaylorN{T},
        x::AbstractArray{HomogeneousPolynomial{T},1},
        u::UnitRange{Int}) where {T<:Number}
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a[u[ind]] = x[ind]
    end
    return a[u]
end
function setindex!(a::TaylorN{T}, x::AbstractArray{T,1},
        u::UnitRange{Int}) where {T<:Number}
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a[u[ind]] = x[ind]
    end
    return a[u]
end
setindex!(a::TaylorN{T}, x::T, ::Colon) where {T<:Number} =
    (a[0:end] = x; a[:])
setindex!(a::TaylorN{T}, x::AbstractArray{HomogeneousPolynomial{T},1},
    ::Colon) where {T<:Number} = (a[0:end] = x; a[:])
setindex!(a::TaylorN{T}, x::AbstractArray{T,1}, ::Colon) where {T<:Number} =
    (a[0:end] = x; a[:])
function setindex!(a::TaylorN{T}, x::T, u::StepRange{Int,Int}) where {T<:Number}
    for ind in u
        a[ind] = x
    end
    return a[u]
end
function setindex!(a::TaylorN{T}, x::AbstractArray{HomogeneousPolynomial{T},1},
        u::StepRange{Int,Int}) where {T<:Number}
    # a[u[:]] .= x[:]
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a[u[ind]] = x[ind]
    end
    return a[u]
end
function setindex!(a::TaylorN{T}, x::Array{T,1},
        u::StepRange{Int,Int}) where {T<:Number}
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a[u[ind]] = x[ind]
    end
    return a[u]
end


## eltype, length, order, etc ##
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval begin
        if $T == HomogeneousPolynomial
            @inline order(a::$T) = a.order
            @inline iterate(a::$T, state=1) =
                state > length(a) ? nothing : (a.coeffs[state], state+1)
            # Base.iterate(rS::Iterators.Reverse{$T}, state=rS.itr.order) = state < 0 ? nothing : (a.coeffs[state], state-1)
            @inline length(a::$T) = a.space.size_table[TS.order(a)+1]
            @inline firstindex(a::$T) = 1
            @inline lastindex(a::$T) = length(a)
        else
            @inline order(a::$T) = size(a.coeffs, 1)-1
            @inline iterate(a::$T, state=0) =
                state > TS.order(a) ? nothing : (a.coeffs[state+1], state+1)
            # Base.iterate(rS::Iterators.Reverse{$T}, state=rS.itr.order) = state < 0 ? nothing : (a.coeffs[state], state-1)
            @inline length(a::$T) = length(a.coeffs)
            @inline firstindex(a::$T) = 0
            @inline lastindex(a::$T) = TS.order(a)
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

# Dumb methods included to properly export normalize_taylor (if IntervalArithmetic is loaded)
@inline normalize_taylor(a::AbstractSeries) = a
@inline aff_normalize(a::AbstractSeries) = a


## _minorder
function _minorder(a, b)
    minorder, maxorder = minmax(TS.order(a), TS.order(b))
    if minorder ≤ 0
        minorder = maxorder
    end
    return minorder
end


## fixorder ##
for T in (:Taylor1, :TaylorN)
    @eval begin
        @inline function fixorder(a::$T, b::$T)
            TS.order(a) == TS.order(b) && return a, b
            minorder = _minorder(a, b)
            return $T(a.coeffs, minorder), $T(b.coeffs, minorder)
        end
    end
end

function fixorder(a::HomogeneousPolynomial, b::HomogeneousPolynomial)
    @assert TS.order(a) == TS.order(b)
    return a, b
end

for T in (:HomogeneousPolynomial, :TaylorN)
    @eval function fixorder(a::Taylor1{$T{T}}, b::Taylor1{$T{S}}) where
            {T<:NumberNotSeries, S<:NumberNotSeries}
        (TS.order(a) == TS.order(b)) && (all(TS.order.(a.coeffs) .== TS.order.(b.coeffs))) && return a, b
        minordT = _minorder(a, b)
        aa = Taylor1(a.coeffs, minordT)
        bb = Taylor1(b.coeffs, minordT)
        for ind in eachindex(aa)
            TS.order(aa[ind]) == TS.order(bb[ind]) && continue
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
# Finds the first non zero entry
function Base.findfirst(a::HomogeneousPolynomial{T}) where {T<:Number}
    first = findfirst(!_isthinzero, a.coeffs)
    isnothing(first) && return -1
    return first
end

# Finds the last non-zero entry
function Base.findlast(a::HomogeneousPolynomial{T}) where {T<:Number}
    last = findlast(!_isthinzero, a.coeffs)
    isnothing(last) && return -1
    return last
end

for T in (:Taylor1, :TaylorN)
    # Finds the first non zero entry
    @eval function Base.findfirst(a::$T{T}) where {T<:Number}
        first = findfirst(!_isthinzero, a.coeffs)
        isnothing(first) && return -1
        return first-1
    end

    # Finds the last non-zero entry
    @eval function Base.findlast(a::$T{T}) where {T<:Number}
        last = findlast(!_isthinzero, a.coeffs)
        isnothing(last) && return -1
        return last-1
    end
end


## copyto! ##
# Inspired from base/abstractarray.jl, line 665
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval function copyto!(dst::$T{T}, src::$T{T}) where {T<:Number}
        length(dst) < length(src) && throw(ArgumentError(string("Destination has fewer elements than required; no copy performed")))
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

constant_term(a::Vector{T}) where {T<:Number} = constant_term.(a)

constant_term(a::Number) = a

"""
    linear_polynomial(a)

Returns the linear part of `a` as a polynomial (`Taylor1` or `TaylorN`),
*without* the constant term. The fallback behavior is to return `a` itself.
"""
linear_polynomial(a::Taylor1) = Taylor1([zero(a[1]), a[1]], TS.order(a))

linear_polynomial(a::HomogeneousPolynomial) =
    HomogeneousPolynomial(a.space, a[1], TS.order(a))

linear_polynomial(a::TaylorN) = TaylorN(a[1], TS.order(a))

linear_polynomial(a::Vector{T}) where {T<:Number} = linear_polynomial.(a)

linear_polynomial(a::Number) = a

"""
    nonlinear_polynomial(a)

Returns the nonlinear part of `a`. The fallback behavior is to return `zero(a)`.
"""
nonlinear_polynomial(a::AbstractSeries) = a - constant_term(a) - linear_polynomial(a)

nonlinear_polynomial(a::Vector{T}) where {T<:Number} = nonlinear_polynomial.(a)

nonlinear_polynomial(a::Number) = zero(a)


"""
    @isonethread (expr)

Internal macro used to check the number of threads in use, to prevent a data race
that modifies `coeff_table` when using `differentiate` or `integrate`; see
https://github.com/JuliaDiff/TaylorSeries.jl/issues/318.

This macro is inspired by the macro `@threaded`; see https://github.com/trixi-framework/Trixi.jl/blob/main/src/auxiliary/auxiliary.jl;
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
