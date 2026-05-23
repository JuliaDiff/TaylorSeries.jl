module TaylorSeriesJLD2Ext

import Base: convert
using TaylorSeries

import JLD2: writeas

@doc raw"""
    JetSpaceSerialization

Lightweight serialization record for the `JetSpace` associated with a `TaylorN`.

The session and space ids are only used to preserve sharing between `TaylorN`s
that belonged to the same `JetSpace` when they were saved.
"""
struct JetSpaceSerialization
    session_id::UInt
    space_id::UInt
    order::Int
    variable_names::Vector{String}
    variable_symbols::Vector{Symbol}
end

@doc raw"""
    TaylorNSerialization{T}

Legacy serialization struct used by older `.jld2` files.

# Fields
- `vars::Vector{String}`: jet transport variables.
- `varorder::Int`: order of jet transport perturbations.
- `x::Vector{T}`: vector of coefficients.
"""
struct TaylorNSerialization{T}
    vars::Vector{String}
    varorder::Int
    x::Vector{T}
end

@doc raw"""
    TaylorNSerializationV2{T}

Custom serialization struct to save a `TaylorN{T}` and its explicit `JetSpace`
metadata to a `.jld2` file.

# Fields
- `space::JetSpaceSerialization`: lightweight description of the associated
  `JetSpace`.
- `varorder::Int`: order of the stored `TaylorN`.
- `x::Vector{T}`: flattened vector of homogeneous coefficients.
"""
struct TaylorNSerializationV2{T}
    space::JetSpaceSerialization
    varorder::Int
    x::Vector{T}
end

const _jld2_session_id = Ref{UInt}(0)
const _write_space_cache = Dict{Tuple{UInt,UInt},JetSpace}()
const _read_space_cache = Dict{Any,JetSpace}()
const _read_space_cache_lock = ReentrantLock()

function _serialization_session_id()
    session_id = _jld2_session_id[]
    iszero(session_id) || return session_id
    session_id = hash((time_ns(), objectid(_read_space_cache)))
    iszero(session_id) && (session_id = one(UInt))
    _jld2_session_id[] = session_id
    return session_id
end

function _space_spec(space::JetSpace)
    session_id = _serialization_session_id()
    space_id = objectid(space)
    _write_space_cache[(session_id, space_id)] = space
    return JetSpaceSerialization(session_id, space_id,
        TS.order(space), copy(TS.get_variable_names(space)),
        copy(TS.get_variable_symbols(space)))
end

function _space_cache_key(spec::JetSpaceSerialization)
    return (spec.session_id, spec.space_id, spec.order,
        Tuple(spec.variable_names), Tuple(spec.variable_symbols))
end

function _reconstruct_space(spec::JetSpaceSerialization)
    variable_names = copy(spec.variable_names)
    variable_symbols = copy(spec.variable_symbols)
    num_vars = length(variable_names)
    num_vars == length(variable_symbols) ||
        error("invalid serialized JetSpace: variable names and symbols differ in length")
    tables = TS.generate_tables(num_vars, spec.order)
    return JetSpace(spec.order, num_vars, variable_names, variable_symbols, tables...)
end

function _cached_space(spec::JetSpaceSerialization)
    saved_space = get(_write_space_cache, (spec.session_id, spec.space_id), nothing)
    saved_space === nothing || return saved_space

    key = _space_cache_key(spec)
    lock(_read_space_cache_lock)
    try
        return get!(_read_space_cache, key) do
            _reconstruct_space(spec)
        end
    finally
        unlock(_read_space_cache_lock)
    end
end

# Tell JLD2 to save TaylorN{T} as TaylorNSerializationV2{T}. The legacy
# TaylorNSerialization reader is kept below so older files remain loadable.
writeas(::Type{TaylorN{T}}) where {T} = TaylorNSerializationV2{T}

# Convert method to write .jld2 files
function convert(::Type{TaylorNSerialization{T}}, eph::TaylorN{T}) where {T}
    # Variables
    vars = TS.get_variable_names()
    # Number of variables
    n = length(vars)
    # TaylorN order
    varorder = TS.order(eph)
    # Number of coefficients in each TaylorN
    L = varorder + 1
    # Number of coefficients in each HomogeneousPolynomial
    M = binomial(n + varorder, varorder)

    # Vector of coefficients
    x = Vector{T}(undef, M)

    # Save coefficients
    i = 1
    for i_1 in 0:varorder
        # Iterate over i_1 order HomogeneousPolynomial
        for i_2 in 1:binomial(n + i_1 - 1, i_1)
            x[i] = eph.coeffs[i_1+1].coeffs[i_2]
            i += 1
        end
    end

    return TaylorNSerialization{T}(vars, varorder, x)
end

# Convert method to write .jld2 files using explicit JetSpace metadata.
function convert(::Type{TaylorNSerializationV2{T}}, eph::TaylorN{T}) where {T}
    varorder = TS.order(eph)
    num_coeffs = 0
    for degree in 0:varorder
        num_coeffs += length(eph.coeffs[degree+1].coeffs)
    end

    x = Vector{T}(undef, num_coeffs)
    i = 1
    for degree in 0:varorder
        coeffs = eph.coeffs[degree+1].coeffs
        for j in eachindex(coeffs)
            x[i] = coeffs[j]
            i += 1
        end
    end

    return TaylorNSerializationV2{T}(_space_spec(eph.space), varorder, x)
end

# Convert method to read .jld2 files
function convert(::Type{TaylorN{T}}, eph::TaylorNSerialization{T}) where {T}
    # Variables
    vars = eph.vars
    # Number of variables
    n = length(vars)
    # TaylorN order
    varorder = eph.varorder
    # Number of coefficients in each TaylorN
    L = varorder + 1
    # Number of coefficients in each HomogeneousPolynomial
    M = binomial(n + varorder, varorder)

    # Set variables
    if TS.get_variable_names() != vars
        TS.variables!(T, vars, order = varorder)
    end

    # Reconstruct TaylorN
    i = 1
    TaylorN_coeffs = Vector{HomogeneousPolynomial{T}}(undef, L)
    for i_1 in 0:varorder
        # Reconstruct HomogeneousPolynomials
        TaylorN_coeffs[i_1 + 1] = HomogeneousPolynomial(eph.x[i : i + binomial(n + i_1 - 1, i_1)-1], i_1)
        i += binomial(n + i_1 - 1, i_1)
    end
    x = TaylorN{T}(TaylorN_coeffs, varorder)

    return x
end

# Convert method to read .jld2 files written with TaylorNSerializationV2.
function convert(::Type{TaylorN{T}}, eph::TaylorNSerializationV2{T}) where {T}
    space = _cached_space(eph.space)
    varorder = eph.varorder
    varorder ≤ TS.order(space) ||
        error("invalid serialized TaylorN: order exceeds serialized JetSpace order")

    i = 1
    TaylorN_coeffs = Vector{HomogeneousPolynomial{T}}(undef, varorder + 1)
    for degree in 0:varorder
        num_coeffs = space.size_table[degree+1]
        last_i = i + num_coeffs - 1
        last_i ≤ length(eph.x) ||
            error("invalid serialized TaylorN: not enough coefficients")
        TaylorN_coeffs[degree+1] =
            HomogeneousPolynomial(space, @view(eph.x[i:last_i]), degree)
        i = last_i + 1
    end
    i == length(eph.x) + 1 ||
        error("invalid serialized TaylorN: too many coefficients")

    return TaylorN(space, TaylorN_coeffs, varorder)
end

end
