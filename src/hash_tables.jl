# This file is part of the TaylorSeries.jl Julia package, MIT license

# Hash tables for HomogeneousPolynomial and TaylorN

"""
    generate_tables(num_vars, order)

Return the hash tables `coeff_table`, `index_table`, `size_table`
and `pos_table`. Internally, these are treated as `const`.

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

function _homogeneous_product_table(index_table, pos_table, order_a::Int,
        order_b::Int)
    order_c = order_a + order_b
    indTa = index_table[order_a+1]
    indTb = index_table[order_b+1]
    posTc = pos_table[order_c+1]
    positions = Matrix{Int}(undef, length(indTb), length(indTa))
    @inbounds for na in eachindex(indTa)
        inda = indTa[na]
        for nb in eachindex(indTb)
            positions[nb, na] = posTc[inda + indTb[nb]]
        end
    end
    return HomogeneousProductTable(positions)
end

function generate_multiplication_tables(index_table, pos_table, order::Int)
    empty_table = HomogeneousProductTable(Matrix{Int}(undef, 0, 0))
    return [[order_a + order_b ≤ order ?
        _homogeneous_product_table(index_table, pos_table, order_a, order_b) :
        empty_table
        for order_b in 0:order] for order_a in 0:order]
end

function TaylorNSpace(order::Int, num_vars::Int, variable_names::Vector{String},
        variable_symbols::Vector{Symbol}, coeff_table::Vector{Vector{Vector{Int}}},
        index_table::Vector{Vector{Int}}, size_table::Vector{Int},
        pos_table::Vector{Dict{Int,Int}})
    mul_table = generate_multiplication_tables(index_table, pos_table, order)
    return TaylorNSpace(order, num_vars, variable_names, variable_symbols,
        coeff_table, index_table, size_table, pos_table, mul_table)
end

function TaylorNSpace(order::Int, variable_names::Vector{String})
    order ≥ 1 || error("Order must be at least 1")
    num_vars = length(variable_names)
    num_vars ≥ 1 || error("Number of variables must be at least 1")
    tables = generate_tables(num_vars, order)
    return TaylorNSpace(order, num_vars, copy(variable_names),
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


const coeff_table, index_table, size_table, pos_table =
    generate_tables(get_numvars(), get_order())

function _sync_legacy_tables!(space::TaylorNSpace)
    resize!(coeff_table, space.order+1)
    resize!(index_table, space.order+1)
    resize!(size_table, space.order+1)
    resize!(pos_table, space.order+1)

    coeff_table[:] = space.coeff_table
    index_table[:] = space.index_table
    size_table[:] = space.size_table
    pos_table[:] = space.pos_table
    return nothing
end

function set_default_space!(space::TaylorNSpace)
    active_space = if isassigned(default_space)
        dst = default_space[]
        dst.order = space.order
        dst.num_vars = space.num_vars
        dst.variable_names = space.variable_names
        dst.variable_symbols = space.variable_symbols
        dst.coeff_table = space.coeff_table
        dst.index_table = space.index_table
        dst.size_table = space.size_table
        dst.pos_table = space.pos_table
        dst.mul_table = space.mul_table
        dst
    else
        default_space[] = space
        space
    end
    _params_TaylorN_.order = active_space.order
    _params_TaylorN_.num_vars = active_space.num_vars
    _params_TaylorN_.variable_names = active_space.variable_names
    _params_TaylorN_.variable_symbols = active_space.variable_symbols
    _sync_legacy_tables!(active_space)
    GC.gc()
    return active_space
end

default_space[] = TaylorNSpace(get_order(), get_numvars(),
    copy(get_variable_names()), copy(get_variable_symbols()),
    coeff_table, index_table, size_table, pos_table)

# Garbage-collect here to free memory
GC.gc();


"""
    show_monomials(ord::Int) --> nothing

List the indices and corresponding of a `HomogeneousPolynomial`
of degree `ord`.
"""
function show_monomials(ord::Int)
    z = zeros(Int, TS.size_table[ord+1])
    for (index, value) in enumerate(TS.coeff_table[ord+1])
        z[index] = 1
        pol = HomogeneousPolynomial(z)
        println(" $index  -->  $(homogPol2str(pol)[4:end])")
        z[index] = 0
    end
    nothing
end
