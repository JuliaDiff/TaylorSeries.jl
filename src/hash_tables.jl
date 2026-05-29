# This file is part of the TaylorSeries.jl Julia package, MIT license

# Hash tables for HomogeneousPolynomial and TaylorN

"""
    generate_tables(num_vars, order)

Return the lookup tables `coeff_table`, `index_table`, `size_table`
and `pos_table` used by a `JetSpace`.

# Hash tables

    coeff_table :: Array{Array{Array{Int,1},1},1}

The ``i+1``-th component contains a vector with the vectors of all the possible
combinations of monomials of a `HomogeneousPolynomial` of order ``i``.

    index_table :: Array{Array{Int,1},1}

The ``i+1``-th component contains a vector of (hashed) indices that represent
the distinct monomials of a `HomogeneousPolynomial` of order (degree) ``i``.

    size_table :: Array{Int,1}

The ``i+1``-th component contains the number of distinct monomials of the
`HomogeneousPolynomial` of order ``i``, equivalent to `length(coeff_table[i])`.

    pos_table :: Array{Dict{Int,Int},1}

The ``i+1``-th component maps the hash index to the (lexicographic) position
of the corresponding monomial in `coeffs_table`.
"""
function generate_tables(num_vars, order)
    coeff_table = [generate_index_vectors(num_vars, i) for i in 0:order]

    index_table = Vector{Int}[map(x->in_base(order, x), coeffs) for coeffs in coeff_table]

    # Check uniqueness of labels as "non-collision" test
    @assert all(allunique.(index_table))

    pos_table = map(make_inverse_dict, index_table)
    size_table = map(length, index_table)

    # The next line tests the consistency of the number of monomials,
    # but it's commented because it may not pass due to the `binomial`
    # @assert sum(size_table) == binomial(num_vars+order, min(num_vars,order))

    return (coeff_table, index_table, size_table, pos_table)
end

"""Return the input-position schedule for multiplying two homogeneous degrees."""
function _homogeneous_product_table(index_table, pos_table, order_a::Int,
        order_b::Int)
    order_c = order_a + order_b
    indTa = index_table[order_a+1]
    indTb = index_table[order_b+1]
    posTc = pos_table[order_c+1]

    num_coeffs_a = length(indTa)
    num_coeffs_b = length(indTb)
    num_pairs = num_coeffs_a * num_coeffs_b
    num_pairs ≤ typemax(UInt32) ||
        error("Product table is too large for UInt32 pair indices")

    input_positions = Vector{Int}(undef, num_pairs)

    pair = 1
    @inbounds for na in eachindex(indTa)
        inda = indTa[na]
        for nb in eachindex(indTb)
            pos = posTc[inda + indTb[nb]]
            input_positions[pair] = pos
            pair += 1
        end
    end

    return HomogeneousProductTable(input_positions, Int[], UInt32[], num_coeffs_b)
end

"""Initialize and return the output-major product schedule for two positive degrees."""
function _init_output_major_product_table!(space::JetSpace, degree_a::Int,
        degree_b::Int)
    table = _product_table(space, degree_a, degree_b)
    !isempty(table.output_pairs) && return table
    num_coeffs_c = space.size_table[degree_a + degree_b + 1]
    num_pairs = length(table.input_positions)
    lock(space.mul_table_lock)
    try
        table = _product_table(space, degree_a, degree_b)
        !isempty(table.output_pairs) && return table
        counts = zeros(Int, num_coeffs_c)
        @inbounds for pair in 1:num_pairs
            counts[table.input_positions[pair]] += 1
        end

        output_offsets = Vector{Int}(undef, num_coeffs_c+1)
        output_offsets[1] = 1
        @inbounds for pos in 1:num_coeffs_c
            output_offsets[pos+1] = output_offsets[pos] + counts[pos]
        end

        output_pairs = Vector{UInt32}(undef, num_pairs)
        cursors = copy(output_offsets)

        @inbounds for pair in 1:num_pairs
            pos = table.input_positions[pair]
            cursor = cursors[pos]
            output_pairs[cursor] = UInt32(pair)
            cursors[pos] = cursor + 1
        end
        table.output_offsets = output_offsets
        table.output_pairs = output_pairs
        return table
    finally
        unlock(space.mul_table_lock)
    end
end

"""Return empty valid positive-degree product-table placeholders for lazy initialization."""
function generate_multiplication_tables(order::Int)
    empty_table = HomogeneousProductTable(Int[], Int[], UInt32[], 0)
    return [[empty_table for _ in 1:(order - degree_a)] for degree_a in 1:(order - 1)]
end

"""Return the cached product table for two positive degrees, initializing it if needed."""
@inline function _product_table(space::JetSpace, degree_a::Int,
        degree_b::Int)
    @boundscheck degree_a > 0 && degree_b > 0 ||
        throw(ArgumentError("product tables are stored only for positive degrees"))
    @boundscheck degree_a + degree_b ≤ space.order ||
        throw(DimensionMismatch("homogeneous product order exceeds JetSpace order"))
    @inbounds table = space.mul_table[degree_a][degree_b]
    !isempty(table.input_positions) && return table
    return _init_product_table!(space, degree_a, degree_b)
end

"""Initialize and cache the input-position product table for two positive degrees."""
function _init_product_table!(space::JetSpace, degree_a::Int,
        degree_b::Int)
    degree_a > 0 && degree_b > 0 ||
        throw(ArgumentError("product tables are stored only for positive degrees"))
    degree_a + degree_b ≤ space.order ||
        throw(DimensionMismatch("homogeneous product order exceeds JetSpace order"))
    lock(space.mul_table_lock)
    try
        @inbounds table = space.mul_table[degree_a][degree_b]
        if isempty(table.input_positions)
            table = _homogeneous_product_table(space.index_table, space.pos_table,
                degree_a, degree_b)
            @inbounds space.mul_table[degree_a][degree_b] = table
        end
        return table
    finally
        unlock(space.mul_table_lock)
    end
end

"""Construct a `JetSpace` from precomputed homogeneous lookup tables."""
function JetSpace(order::Int, num_vars::Int, variable_names::Vector{String},
        variable_symbols::Vector{Symbol}, coeff_table::Vector{Vector{Vector{Int}}},
        index_table::Vector{Vector{Int}}, size_table::Vector{Int},
        pos_table::Vector{Dict{Int,Int}})
    mul_table = generate_multiplication_tables(order)
    return JetSpace(order, num_vars, variable_names, variable_symbols,
        coeff_table, index_table, size_table, pos_table, mul_table, ReentrantLock())
end

"""Construct a `JetSpace` and generate its lookup tables from variable names."""
function JetSpace(order::Int, variable_names::Vector{String})
    order ≥ 1 || error("Order must be at least 1")
    num_vars = length(variable_names)
    num_vars ≥ 1 || error("Number of variables must be at least 1")
    tables = generate_tables(num_vars, order)
    return JetSpace(order, num_vars, copy(variable_names),
        Symbol.(variable_names), tables...)
end

"""
    generate_index_vectors(num_vars, degree)

Return a vector of index vectors with `num_vars` (number of variables) and
degree.
"""
function generate_index_vectors(num_vars, degree)
    if num_vars == 1
        return Vector{Int}[ [degree] ]
    end

    indices = Vector{Int}[]

    for k in degree:-1:0

        new_indices = [ [k, x...] for x in generate_index_vectors(num_vars-1, degree-k) ]
        append!(indices, new_indices)

    end

    return indices
end


function make_forward_dict(v::Vector)
    Dict(Dict(i=>x for (i,x) in enumerate(v)))
end

"""
    make_inverse_dict(v)

Return a Dict with the enumeration of `v`: the elements of `v` point to
the corresponding index.

It is used to construct `pos_table` from `index_table`.
"""
make_inverse_dict(v::Vector) = Dict(Dict(x=>i for (i,x) in enumerate(v)))

"""
    in_base(order, v)

Convert vector `v` of non-negative integers to base `oorder`, where
`oorder` is the next odd integer of `order`.
"""
function in_base(order, v)
    oorder = iseven(order) ? order+1 : order+2 # `oorder` is the next odd integer to `order`

    result = 0

    all(iszero.(v)) && return result

    for i in v
        result = result*oorder + i
    end

    return result
end


"""
Set the compatibility default `JetSpace`.

The global `default_space` binding is a constant `Ref`, but its contents are
mutable. Module loading initializes `default_space[]`; later compatibility calls
such as `variables!` replace `default_space[]` with `space`. Existing `TaylorN`
and `HomogeneousPolynomial` objects keep their original spaces, while future
default-space constructors use the new default algebra.

Use `nowarn=true` to suppress the warning emitted when replacing the default
space.
"""
function set_default_space!(space::JetSpace; nowarn::Bool=false)
    old_space = default_space[]
    old_space === space && return space

    msg = "Updating TaylorSeries.default_space[]; existing TaylorN and " *
        "HomogeneousPolynomial objects keep their original JetSpace, while " *
        "future default-space constructors use the new default JetSpace."
    old_order = order(old_space)
    old_numvars = get_numvars(old_space)
    new_order = order(space)
    new_numvars = get_numvars(space)
    nowarn || @warn msg old_order old_numvars new_order new_numvars

    # Replace the active default space. Do not mutate the old object in place:
    # existing TaylorN/HomogeneousPolynomial objects may still depend on it.
    default_space[] = space

    # The previous default space may own large lookup tables. If no existing
    # series still references it, this gives the GC a chance to reclaim them.
    GC.gc()
    return space
end

default_space[] = JetSpace(DEFAULT_TAYLORN_ORDER, copy(DEFAULT_TAYLORN_VARIABLE_NAMES))

# Garbage-collect here to free memory
GC.gc();


"""
    show_monomials(ord::Int) --> nothing
    show_monomials(space::JetSpace, ord::Int) --> nothing

List the indices and corresponding monomials of a `HomogeneousPolynomial`
of degree `ord`, using either the default `JetSpace` or the given explicit
`space`.
"""
show_monomials(ord::Int) = show_monomials(default_space[], ord)

function show_monomials(space::JetSpace, ord::Int)
    z = zeros(Int, space.size_table[ord+1])
    for index in eachindex(space.coeff_table[ord+1])
        z[index] = 1
        pol = HomogeneousPolynomial(space, z, ord)
        println(" $index  -->  $(homogPol2str(pol)[4:end])")
        z[index] = 0
    end
    nothing
end
