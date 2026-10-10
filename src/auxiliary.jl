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

@noinline function _space_mismatch_error(space_a::JetSpace, space_b::JetSpace)
    throw(ArgumentError(
        "JetSpace mismatch: operands belong to different spaces. " *
        "Use an explicit projection or conversion before combining them."))
end


# Does the element type carry a JetSpace? (compile-time)
_has_space(::Type) = false
_has_space(::Type{<:Union{HomogeneousPolynomial,TaylorN}}) = true
_has_space(::Type{Taylor1{T}}) where {T} = _has_space(T)


## Scalar-space helpers ------------------------------------------------------
# Plain numeric values created by conversion live in `_scalar_space[]` (order 0, 0
# variables). The *raw* `_check_same_space(::JetSpace, ::JetSpace)` stays strict
# on purpose: in-place kernels run `@inbounds` loops sized by one operand, so a
# scalar leaking into them must be an error. Scalars are instead embedded at the
# entry points (`_unify_space`) and by containers (`_adopt`, `_checked_coeffs`).
#
# Naming: `_jetspace(x)` gives the JetSpace of `x` (the scalar space for numbers and
# for series with no non-scalar elements), `_is_scalar_space(x)` asks whether `x` is
# space-agnostic, and `_embed_scalar(x, sp[, ord])` returns `x` in `sp` *only if*
# `_is_scalar_space(x)`; otherwise it returns `x` itself, so a non-scalar space is never
# silently changed.

# The space of an object (recursive through nested Taylor1s)
_jetspace(::Number) = _scalar_space[]                 # plain numbers
_jetspace(a::Union{HomogeneousPolynomial,TaylorN}) = space(a)
function _jetspace(a::Taylor1{T}) where {T}
    _has_space(T) || return _scalar_space[]
    for c in a.coeffs
        sp = _jetspace(c)
        _is_scalar_space(sp) || return sp             # first non-scalar coefficient
    end
    return _scalar_space[]
end

# Is the object space-agnostic? For a `Taylor1{TaylorN}`: are *all* coefficients
@inline _is_scalar_space(x::Number) = _is_scalar_space(_jetspace(x))


"""
    _embed_scalar(x, sp::JetSpace, ord::Int=0)

Return `x` rebuilt in `sp` if `x` is space-agnostic (`_is_scalar_space(x)`), and `x`
otherwise. `ord` is the inner order given to the embedded `TaylorN`s
(0 means "take the order of the other operand", as in `fixorder`); it is
ignored by the other types. It also acts on vectors and tuples of series, which are
returned unchanged if nothing has to be embedded.
"""
_embed_scalar(x, ::JetSpace, ::Int=0) = x                          # nothing to embed
_embed_scalar(a::HomogeneousPolynomial, sp::JetSpace, ::Int=0) =
    _is_scalar_space(a) ? HomogeneousPolynomial(sp, a.coeffs[1], 0) : a
_embed_scalar(a::TaylorN, sp::JetSpace, ord::Int=0) =
    _is_scalar_space(a) ? TaylorN(sp, a.coeffs[1].coeffs[1], ord) : a
function _embed_scalar(a::Taylor1{TaylorN{T}}, sp::JetSpace, ord::Int=0) where {T<:Number}
    _is_scalar_space(a) || return a          # only a Taylor1 made of constants
    v = FixedSizeVectorDefault{TaylorN{T}}(undef, length(a.coeffs))
    for (i, c) in enumerate(a.coeffs)
        v[i] = _embed_scalar(c, sp, ord)
    end
    return Taylor1{TaylorN{T}}(v)
end
function _embed_scalar(v::Union{AbstractVector{<:Union{HomogeneousPolynomial,TaylorN}},Tuple},
        sp::JetSpace, ord::Int=0)
    _is_scalar_space(sp) && return v
    any(_is_scalar_space, v) || return v
    return map(x -> _embed_scalar(x, sp, ord), v)
end


# Order of the first non-scalar `TaylorN` (used to embed evaluation points; containers
# keep embedded constants at order 0, which acts as a wildcard in `fixorder`)
_reference_order(v) = 0
function _reference_order(v::Union{AbstractVector{<:TaylorN},Tuple})
    for x in v
        x isa TaylorN && !_is_scalar_space(x) && return order(x)
    end
    return 0
end


"""
    _common_space(sp::JetSpace, v)

JetSpace shared by the non-scalar entries of `v` (and `sp`, if it is not the scalar
space); the scalar space if there is none. Throws if two non-scalar spaces differ.
"""
@inline function _common_space(sp::JetSpace, v)
    ssp = _scalar_space[]
    for x in v
        s = _jetspace(x)
        if s === ssp || s === sp
            continue
        elseif sp === ssp
            sp = s
        else
            _check_same_space(sp, s)     # two different non-scalar spaces: throws
        end
    end
    return sp
end


"""
    _adopt(sp::JetSpace, x[, ord])

Value to be *stored* in a container living in `sp`: space-agnostic `x` is embedded
(a new object), anything else is checked against `sp` and copied (no aliasing).
"""
function _adopt(sp::JetSpace, x::TaylorN, ord::Int=0)
    _is_scalar_space(x) && return _embed_scalar(x, sp, ord)
    _check_same_space(sp, x.space)
    return TaylorN(sp, x.coeffs, order(x))     # space known: skips the 2-argument constructor path
end
function _adopt(sp::JetSpace, x::HomogeneousPolynomial, ::Int=0)
    _is_scalar_space(x) && return _embed_scalar(x, sp)
    _check_same_space(sp, x.space)
    return _copy_series(x)
end
function _adopt(sp::JetSpace, x::Taylor1, ord::Int=0)
    _is_scalar_space(x) || _check_same_space(sp, _jetspace(x))
    y = _embed_scalar(x, sp, ord)
    return y === x ? deepcopy(x) : y
end


# Coefficients for a `Taylor1` built from a user-supplied vector. Numbers need nothing
# (the constructors copy the container). For series, scalar-space entries are embedded
# and an object that appears in several slots (`[ξ, ξ, ξ]`, `fill(ξ, n)`) is copied in
# all but its first occurrence, so that no two slots share an object. Objects that
# appear once are stored as they are (no copy), as before.
_own_coeffs(v::AbstractVector{<:Number}) = v
function _own_coeffs(v::AbstractVector{T}) where {T<:AbstractSeries}
    n = length(v)
    out = FixedSizeVectorDefault{T}(undef, n)
    hasspace = _has_space(T)
    sp = hasspace ? _common_space(_scalar_space[], v) : _scalar_space[]
    dups = _repeated_objects(v)
    for (i, x) in enumerate(v)
        y = (hasspace && _is_scalar_space(x)) ? _embed_scalar(x, sp) : x
        out[i] = (y === x && dups !== nothing && dups[i]) ?
            (hasspace ? _adopt(sp, x) : deepcopy(x)) : y
    end
    return out
end


# Returns `nothing` if no object appears twice in `v`; otherwise a `BitVector` marking the
# entries that repeat an earlier one. O(n²) pointer comparisons for short vectors (no allocation
# when there are no repeated objects); above that, a small open-addressing table keyed by
# the address of the coefficient storage (exact: a hit is confirmed with `===`), which is
# O(n). Cheap identity hash for series (`Taylor1`, `HomogeneousPolynomial`, `TaylorN`):
# the address of their coefficients (not dereferenced; objects are alive during the call)
@inline _storage_hash(x) = UInt(pointer(x.coeffs)) >> 4

function _repeated_objects(v)
    n = length(v)
    mask = nothing
    if n <= 24
        for i in 2:n
            x = v[i]
            for j in 1:i-1
                if @inbounds v[j] === x
                    mask === nothing && (mask = falses(n))
                    mask[i] = true
                    break
                end
            end
        end
        return mask
    end
    m = nextpow(2, 2n)
    table = zeros(Int, m)          # slot -> index (in `v`) of the object stored there
    for i in 1:n
        x = v[i]
        h = Int(_storage_hash(x) & UInt(m - 1)) + 1
        while true
            j = table[h]
            if j == 0
                table[h] = i
                break
            elseif v[j] === x
                mask === nothing && (mask = falses(n))
                mask[i] = true
                break
            end
            h = h == m ? 1 : h + 1
        end
    end
    return mask
end


# Used by the inner `TaylorN` constructor: one pass that throws on a non-scalar space
# different from `sp`, and embeds scalar-space polynomials (new vector) if there are any
@inline function _checked_hps(sp::JetSpace, v::AbstractVector{HomogeneousPolynomial{T}}) where {T}
    ssp = _scalar_space[]
    nscalar = 0
    for pol in v
        s = pol.space
        s === sp && continue
        s === ssp ? (nscalar += 1) : _check_same_space(sp, s)
    end
    (nscalar == 0 || sp === ssp) && return v
    return _embed_entries(sp, v)
end

# Rare path of the checks above; returns a copy of the same container type (type-stable)
@noinline function _embed_entries(sp::JetSpace, v::AbstractVector)
    out = copy(v)
    for i in eachindex(out)
        out[i] = _embed_scalar(out[i], sp)
    end
    return out
end

# Trusted constructor: `v` already holds objects nobody else references
_taylor1_owned(v::AbstractVector{T}) where {T<:Number} =
    Taylor1{T}(v isa FixedSizeVectorDefault{T} ? v : FixedSizeVectorDefault(v))

# Used by the inner `Taylor1` constructors. Throws on mixed non-scalar spaces; if scalar-space
# entries are mixed with non-scalar ones, returns a new vector with them embedded
# (the caller's vector is not modified). Otherwise returns `v` itself.
@inline _checked_coeffs(v::AbstractVector{T}) where {T} =
    _has_space(T) ? _checked_coeffs_space(v) : v
@inline function _checked_coeffs_space(v::AbstractVector{T}) where {T}
    isempty(v) && return v
    # single pass: common non-scalar space (strict) and number of scalar-space entries
    ssp = _scalar_space[]
    sp = ssp
    nscalar = 0
    for x in v
        s = _jetspace(x)
        if s === ssp
            nscalar += 1
        elseif sp === ssp
            sp = s
        elseif s !== sp
            _check_same_space(sp, s)     # two different non-scalar spaces: throws
        end
    end
    (nscalar == 0 || sp === ssp) && return v
    return _embed_entries(sp, v)
end

# Entry-point helper for binary operations: embed the scalar operand (if any) into
# the space of the other one, error on two different non-scalar spaces, and keep the full
# coefficient check for `Taylor1{<:TaylorN}` (kernels may write `coeffs` directly).
@inline _unify_space(a::AbstractSeries, b::AbstractSeries) =
    (_check_same_space(a, b); (a, b))
@inline function _unify_space(a::HomogeneousPolynomial, b::HomogeneousPolynomial)
    sa, sb = a.space, b.space
    sa === sb && return a, b
    _is_scalar_space(sa) && return _embed_scalar(a, sb), b
    _is_scalar_space(sb) && return a, _embed_scalar(b, sa)
    _space_mismatch_error(sa, sb)
end
@inline function _unify_space(a::TaylorN, b::TaylorN)
    sa, sb = a.space, b.space
    sa === sb && return a, b
    _is_scalar_space(sa) && return _embed_scalar(a, sb, order(b)), b
    _is_scalar_space(sb) && return a, _embed_scalar(b, sa, order(a))
    _space_mismatch_error(sa, sb)
end
# `Taylor1{<:TaylorN}` with `TaylorN`; a `Taylor1` made of constants (e.g. built from a
# type, `Taylor1(TaylorN{Float64}, 5)`) takes the space of the other operand
function _unify_space(a::Taylor1{<:TaylorN}, b::TaylorN)
    sa, sb = _jetspace(a), b.space
    if sa !== sb
        if _is_scalar_space(sa)
            a = _embed_scalar(a, sb, order(b))
        elseif _is_scalar_space(sb)
            b = _embed_scalar(b, sa, order(a.coeffs[1]))
        else
            _space_mismatch_error(sa, sb)
        end
    end
    _check_same_space(a.coeffs[1], b)
    return a, b
end
_unify_space(b::TaylorN, a::Taylor1{<:TaylorN}) = reverse(_unify_space(a, b))

function _unify_space(a::Taylor1{<:TaylorN}, b::Taylor1{<:TaylorN})
    sa, sb = _jetspace(a), _jetspace(b)
    if sa !== sb
        if _is_scalar_space(sa)
            a = _embed_scalar(a, sb, order(b.coeffs[1]))
        elseif _is_scalar_space(sb)
            b = _embed_scalar(b, sa, order(a.coeffs[1]))
        else
            _space_mismatch_error(sa, sb)
        end
    end
    _check_same_space(a, b)     # full-coefficient check (still needed, item 4)
    return a, b
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
    _check_same_space(a::Taylor1{<:TaylorN}[, b::Taylor1{<:TaylorN}[, c::Taylor1{<:TaylorN}]])
    _check_same_space(space::JetSpace,v::AbstractVector{<:HomogeneousPolynomial})
    _check_same_space(v::AbstractVector)

Throw an `ArgumentError` unless all arguments belong to the same `JetSpace`
by object identity. For `Taylor1{<:TaylorN}` arguments, *every* coefficient
is checked; for other `Taylor1`s the check is a no-op.
"""
@inline function _check_same_space(space_a::JetSpace, space_b::JetSpace)
    space_a === space_b || _space_mismatch_error(space_a, space_b)
    return nothing
end
@inline _check_same_space(::Taylor1, ::Taylor1) = nothing
@inline _check_same_space(::Taylor1, ::Taylor1, ::Taylor1) = nothing
function _check_coeffs_space(sp::JetSpace, a::Taylor1{<:TaylorN},
        inds=eachindex(a))
    for k in inds
        _check_same_space(sp, space(a[k]))
    end
    return nothing
end
_check_same_space(a::Taylor1{<:TaylorN}) = _check_same_space(a.coeffs)

function _check_same_space(a::Taylor1{<:TaylorN}, b::Taylor1{<:TaylorN})
    sp = space(a[0])
    _check_coeffs_space(sp, a)
    _check_coeffs_space(sp, b)
    return nothing
end
function _check_same_space(a::Taylor1{<:TaylorN}, b::Taylor1{<:TaylorN},
        c::Taylor1{<:TaylorN})
    sp = space(a[0])
    _check_coeffs_space(sp, a)
    _check_coeffs_space(sp, b)
    _check_coeffs_space(sp, c)
    return nothing
end
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

_check_same_space(v::AbstractVector) =
    _has_space(eltype(v)) ? _check_vector_space(v) : nothing


function _check_vector_space(v::AbstractVector)
    isempty(v) && return nothing
    sp = _jetspace(first(v))
    for x in v
        _check_same_space(sp, _jetspace(x))
    end
    return nothing
end


"""
    _check_same_space_all(a::AbstractSeries, vals)
    _check_same_space_all(::AbstractVector)

Check that every element of `vals` (a tuple or vector of `HomogeneousPolynomial`
or `TaylorN`) belongs to the same `JetSpace` as `a`.
"""
function _check_same_space_all(a::AbstractSeries, vals)
    sp = _jetspace(a)
    for v in vals
        _check_same_space(sp, _jetspace(v))
    end
    return nothing
end
_check_same_space_all(v::AbstractVector) = (_has_space(eltype(v)) &&
    !isempty(v)) ? _check_same_space_all(first(v), v) : nothing


function _space_from_homogeneous_vector(v::AbstractVector{<:HomogeneousPolynomial},
        fallback::JetSpace)
    isempty(v) && return fallback
    return _common_space(_scalar_space[], v)
end

_constant_series_like(a::Taylor1, x, order::Int) = Taylor1(x, order)
_constant_series_like(a::TaylorN, x, order::Int) = TaylorN(a.space, x, order)


_copy_series(x::Taylor1) = Taylor1(x.coeffs, order(x))
_copy_series(x::HomogeneousPolynomial) = HomogeneousPolynomial(x.space, x.coeffs[:], order(x))
_copy_series(x::TaylorN) = TaylorN(x.coeffs, order(x))


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
    for i = 1:order(space)+1
        ll ≤ space.size_table[i] && return ord
        ord += 1
    end
    return ord
end
orderH(coeffs::AbstractArray{T,1}) where {T<:Number} =
    orderH(default_space[], coeffs)

## Maximum order of a HomogeneousPolynomial vector; used by TaylorN constructor
maxorderH(v::AbstractArray{HomogeneousPolynomial{T},1}) where {T<:Number} =
    isempty(v) ? 0 : maximum(order.(v))


"""
    _evaluation_order(a, x)

Return the smallest positive order among `x` and the coefficients of `a`.
Return zero if all those orders are zero. Order-zero series are treated as
exact constants, so they do not limit the order of the evaluated result.
"""
function _evaluation_order(a::Taylor1{<:AbstractSeries}, x::AbstractSeries)
    ord = order(x)
    for coeff in a.coeffs
        coeff_order = order(coeff)
        # Only positive orders limit the information available in the result.
        if ord == 0
            ord = coeff_order
        elseif coeff_order > 0
            ord = min(ord, coeff_order)
        end
    end
    return ord
end


## getcoeff ##
"""
    getcoeff(a, n)

Return the coefficient of order `n::Int` of a `a::Taylor1` polynomial; the constant
term corresponds to n=0.
"""
getcoeff(a::Taylor1, n::Int) = (@assert 0 ≤ n ≤ order(a); return a[n])

@inline getindex(a::Taylor1, n::Int) = a.coeffs[n+1]
getindex(a::Taylor1, u::UnitRange{Int}) = view(a.coeffs, u .+ 1 )
getindex(a::Taylor1, c::Colon) = view(a.coeffs, c)
getindex(a::Taylor1{T}, u::StepRange{Int,Int}) where {T<:Number} =
    view(a.coeffs, u .+ 1)

@inline setindex!(a::Taylor1{T}, x::T, n::Int) where {T<:NumberNotSeries} =
    a.coeffs[n+1] = x
# setindex!(a::Taylor1{T}, x::T, n::Int) where {T<:AbstractSeries} =
#     setindex!(a.coeffs, deepcopy(x), n+1)
@inline function setindex!(a::Taylor1{TaylorN{T}}, x::TaylorN{T}, n::Int) where
        {T<:NumberNotSeries}
    sp = _jetspace(a)
    if _is_scalar_space(sp) && !_is_scalar_space(x)
        sp = x.space           # `a` was made of constants: its coefficients adopt that space
        for i in eachindex(a.coeffs)
            a.coeffs[i] = _embed_scalar(a.coeffs[i], sp)
        end
    end
    return a.coeffs[n+1] = _adopt(sp, x, order(a.coeffs[n+1]))
end
# Build the `HomogeneousPolynomial` explicitly in `a.space`: storing the `Taylor1` directly
# would go through `convert(HomogeneousPolynomial{...}, ::Taylor1)`, which only sees the
# target type (scalar space). `_coeffsHP` copies `x`.
@inline setindex!(a::TaylorN{Taylor1{T}}, x::Taylor1{T}, n::Int) where
    {T<:NumberNotSeries} = a.coeffs[n+1] = HomogeneousPolynomial(a.space, x, n)
@inline function setindex!(a::Taylor1{Taylor1{T}}, x::Taylor1{T}, n::Int) where
        {T<:Taylor1{<:Number}}
    a.coeffs[n+1] = zero(x)
    for i in eachindex(x)
        a.coeffs[n+1].coeffs[i+1] = x.coeffs[i+1]
    end
    return a.coeffs[n+1]
end
@inline setindex!(a::Taylor1{Taylor1{T}}, x::Taylor1{T}, n::Int) where
    {T<:NumberNotSeries} = a.coeffs[n+1] = Taylor1(x.coeffs[:], order(x))
@inline setindex!(a::Taylor1{T}, x::T, u::UnitRange{Int}) where {T<:Number} =
    a.coeffs[u .+ 1] .= x
@inline function setindex!(a::Taylor1{T}, x::AbstractArray{T,1},
        u::UnitRange{Int}) where {T<:Number}
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a.coeffs[u[ind]+1] = x[ind]
    end
end
@inline setindex!(a::Taylor1{T}, x::T, c::Colon) where {T<:Number} = a.coeffs[c] .= x
@inline setindex!(a::Taylor1{T}, x::AbstractArray{T,1}, c::Colon) where {T<:Number} =
    a.coeffs[c] .= x
@inline setindex!(a::Taylor1{T}, x::T, u::StepRange{Int,Int}) where {T<:Number} =
    a.coeffs[u[:] .+ 1] .= x
function setindex!(a::Taylor1{T}, x::Array{T,1}, u::StepRange{Int,Int}) where {T<:Number}
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a.coeffs[u[ind]+1] = x[ind]
    end
end
# Range and colon assignment for Taylor1{TaylorN}: route through the scalar
# method above, so each entry is space-checked and copied (no aliasing).
# The value types match the generic methods to avoid ambiguities.
_taylor1_indices(a::Taylor1, u::AbstractRange{Int}) = u
_taylor1_indices(a::Taylor1, ::Colon) = eachindex(a)
for I in (:(UnitRange{Int}), :(StepRange{Int,Int}), :Colon)
    @eval function setindex!(a::Taylor1{TaylorN{T}}, x::TaylorN{T}, u::$I) where
            {T<:NumberNotSeries}
        for k in _taylor1_indices(a, u)
            a[k] = x
        end
        return x
    end
end
for (I, V) in ((:(UnitRange{Int}), :AbstractArray), (:(StepRange{Int,Int}), :Array),
        (:Colon, :AbstractArray))
    @eval function setindex!(a::Taylor1{TaylorN{T}}, x::$V{TaylorN{T},1}, u::$I) where
            {T<:NumberNotSeries}
        idx = _taylor1_indices(a, u)
        @assert length(idx) == length(x)
        for (k, xk) in zip(idx, x)
            a[k] = xk
        end
        return x
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
    kdic = in_base(order(sp), v)
    @inbounds n = sp.pos_table[order(a)+1][kdic]
    a[n]
end
getcoeff(a::HomogeneousPolynomial, v::AbstractArray{Int,1}) =
    getcoeff(a, (v...,))

@inline getindex(a::HomogeneousPolynomial, n::Int) = a.coeffs[n]
getindex(a::HomogeneousPolynomial, n::UnitRange{Int}) = view(a.coeffs, n)
getindex(a::HomogeneousPolynomial, c::Colon) = view(a.coeffs, c)
getindex(a::HomogeneousPolynomial, u::StepRange{Int,Int}) = view(a.coeffs, u[:])

@inline setindex!(a::HomogeneousPolynomial{T}, x::T, n::Int) where {T<:Number} =
    a.coeffs[n] = x
@inline setindex!(a::HomogeneousPolynomial{T}, x::T, n::UnitRange{Int}) where
    {T<:Number} = a.coeffs[n] .= x
@inline setindex!(a::HomogeneousPolynomial{T}, x::AbstractArray{T,1},
    n::UnitRange{Int}) where {T<:Number} = a.coeffs[n] .= x
@inline setindex!(a::HomogeneousPolynomial{T}, x::T, c::Colon) where {T<:Number} =
    a.coeffs[c] .= x
@inline setindex!(a::HomogeneousPolynomial{T}, x::AbstractArray{T,1},
    c::Colon) where {T<:Number} = a.coeffs[c] .= x
@inline setindex!(a::HomogeneousPolynomial{T}, x::T,
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

@inline getindex(a::TaylorN, n::Int) = a.coeffs[n+1]
@inline getindex(a::TaylorN, u::UnitRange{Int}) = view(a.coeffs, u .+ 1)
@inline getindex(a::TaylorN, c::Colon) = view(a.coeffs, c)
@inline getindex(a::TaylorN, u::StepRange{Int,Int}) = view(a.coeffs, u[:] .+ 1)

@inline function setindex!(a::TaylorN{T}, x::HomogeneousPolynomial{T}, n::Int) where
        {T<:Number}
    @assert order(x) == n
    return a.coeffs[n+1] = _adopt(a.space, x)
end
@inline setindex!(a::TaylorN{T}, x::T, n::Int) where {T<:Number} =
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
            @inline length(a::$T) = a.space.size_table[order(a)+1]
            @inline firstindex(a::$T) = 1
            @inline lastindex(a::$T) = length(a)
        else
            @inline order(a::$T) = size(a.coeffs, 1)-1
            @inline iterate(a::$T, state=0) =
                state > order(a) ? nothing : (a.coeffs[state+1], state+1)
            # Base.iterate(rS::Iterators.Reverse{$T}, state=rS.itr.order) = state < 0 ? nothing : (a.coeffs[state], state-1)
            @inline length(a::$T) = length(a.coeffs)
            @inline firstindex(a::$T) = 0
            @inline lastindex(a::$T) = order(a)
        end
        @inline eachindex(a::$T) = firstindex(a):lastindex(a)
        @inline numtype(::$T{S}) where {S<:Number} = S
        @inline size(a::$T) = size(a.coeffs)
        @inline axes(a::$T) = ()
    end
end
@inline numtype(a) = eltype(a)

@doc doc"""
    numtype(a::AbstractSeries)

Returns the type of the elements of the coefficients of `a`.
""" numtype

# Dumb methods included to properly export normalize_taylor (if IntervalArithmetic is loaded)
@inline normalize_taylor(a::AbstractSeries) = a
@inline aff_normalize(a::AbstractSeries) = a


## _minorder
function _minorder(a, b)
    minorder, maxorder = minmax(order(a), order(b))
    if minorder ≤ 0
        minorder = maxorder
    end
    return minorder
end


## fixorder ##
for T in (:Taylor1, :TaylorN)
    @eval begin
        @inline function fixorder(a::$T, b::$T)
            order(a) == order(b) && return a, b
            minorder = _minorder(a, b)
            return $T(a.coeffs, minorder), $T(b.coeffs, minorder)
        end
    end
end

function fixorder(a::HomogeneousPolynomial, b::HomogeneousPolynomial)
    @assert order(a) == order(b)
    return a, b
end

for T in (:HomogeneousPolynomial, :TaylorN)
    @eval function fixorder(a::Taylor1{$T{T}}, b::Taylor1{$T{S}}) where
            {T<:NumberNotSeries, S<:NumberNotSeries}
        (order(a) == order(b)) && (all(order.(a.coeffs) .== order.(b.coeffs))) && return a, b
        minordT = _minorder(a, b)
        aa = Taylor1(a.coeffs, minordT)
        bb = Taylor1(b.coeffs, minordT)
        for ind in eachindex(aa)
            order(aa[ind]) == order(bb[ind]) && continue
            minordQ = _minorder(aa[ind], bb[ind])
            aa[ind] = $T(space(aa[ind]), aa[ind].coeffs, minordQ)
            bb[ind] = $T(space(bb[ind]), bb[ind].coeffs, minordQ)
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
    first = findfirst(!_isthinzero, view(a.coeffs, :))
    first = isnothing(first) ? -1 : first
    return first
end

# Finds the last non-zero entry
function Base.findlast(a::HomogeneousPolynomial{T}) where {T<:Number}
    last = findlast(!_isthinzero, view(a.coeffs, :))
    last = isnothing(last) ? -1 : last
    return last
end

for T in (:Taylor1, :TaylorN)
    # Finds the first non zero entry
    @eval function Base.findfirst(a::$T{T}) where {T<:Number}
        first = findfirst(!_isthinzero, view(a.coeffs, :))
        first = isnothing(first) ? 0 : first
        return first-1
    end

    # Finds the last non-zero entry
    @eval function Base.findlast(a::$T{T}) where {T<:Number}
        last = findlast(!_isthinzero, view(a.coeffs, :))
        last = isnothing(last) ? 0 : last
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
@inline constant_term(a::Taylor1) = a.coeffs[1]

@inline constant_term(a::TaylorN) = a.coeffs[1].coeffs[1]

constant_term(a::Vector{T}) where {T<:Number} = constant_term.(a)

@inline constant_term(a::Number) = a

"""
    constant_term!(a, c)

Update the constant term (zero order coefficient) of `a` to `c`,
leaving higher order coefficients unchanged.
"""
@inline function constant_term!(a::Taylor1{T}, c::T) where {T<:Number}
    a.coeffs[1] = c
    return a
end

@inline function constant_term!(a::TaylorN{T}, c::T) where {T<:Number}
    a.coeffs[1].coeffs[1] = c
    return a
end

@inline function constant_term!(a::HomogeneousPolynomial{T}, c::T) where {T<:Number}
    iszero(order(a)) ||
        throw(ArgumentError("only zero-order HomogeneousPolynomial has a constant term"))
    a.coeffs[1] = c
    return a
end

"""
    linear_polynomial(a)

Returns the linear part of `a` as a polynomial (`Taylor1` or `TaylorN`),
*without* the constant term. The fallback behavior is to return `a` itself.
"""
linear_polynomial(a::Taylor1) = Taylor1([zero(a[1]), a[1]], order(a))

linear_polynomial(a::HomogeneousPolynomial) =
    HomogeneousPolynomial(a.space, a[1], order(a))

linear_polynomial(a::TaylorN) = TaylorN(space(a), a[1], order(a))

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
that modifies coefficient tables when using `differentiate` or `integrate`; see
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