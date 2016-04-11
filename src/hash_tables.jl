# This file is part of the TaylorSeries.jl Julia package, MIT license

# Hash tables for HomogeneousPolynomial and TaylorN

@doc """
    generate_tables(num_vars, order)

Return the hash tables `coeff_table`, `index_table`, `size_table` and
`pos_table`. See the help of each of them for specific information.
""" ->
function generate_tables(num_vars, order)
    coeff_table = [generate_index_vectors(num_vars, i) for i in 0:order]

    index_table = Vector{Int}[map(x->in_base(order, x), coeffs) for coeffs in coeff_table]

    pos_table = map(make_inverse_dict, index_table)
    size_table = map(length, index_table)

    coeff_table, index_table, size_table, pos_table
end

@doc """
    generate_index_vectors(num_vars, degree)

Return a vector of index vectors with `num_vars` (number of variables) and
degree.
""" ->
function generate_index_vectors(num_vars, degree)
    if num_vars == 1
        return Vector{Int}[ [degree] ]
    end

    indices = Vector{Int}[]

    for k in degree:-1:0

        new_indices = [ [k, x...] for x in generate_index_vectors(num_vars-1, degree-k) ]
        append!(indices, new_indices)

    end

    indices
end


function make_forward_dict(v::Vector)
    Dict([i=>x for (i,x) in enumerate(v)])
end

@doc """
    make_inverse_dict(v)

Return a Dict with the enumeration of `v`: the elements of `v` point to
the corresponding index.

It is used to construct `pos_table` from `index_table`.
""" ->
function make_inverse_dict(v::Vector)
    Dict([x=>i for (i,x) in enumerate(v)])
end

@doc """
    in_base(order, v)

Convert vector `v` of non-negative integers to base `order+1`.
""" ->
function in_base(order, v)
    order = order+1

    result = 0

    for i in v
        result = result*order + i
    end

    result
end


const coeff_table, index_table, size_table, pos_table =
    generate_tables(get_numvars(), get_order())
gc();

# Documentation of the hash tables
@doc """
    coeff_table::Array{Array{Array{Int64,1},1},1}

The `i+1`-th component contains a vector with vectors holding all the possible
degrees of the independent variables of the monomials of a
`HomogeneousPolynomial` of order `i`.
""" coeff_table

@doc """
    index_table ::Array{Array{Int64,1},1}

The `i+1`-th component contains a vector of (hashed) indices that represent
the distinct monomials of a `HomogeneousPolynomial` of order (degree) `i`;
see `coeff_table`.
""" index_table

@doc """
    size_table::Array{Int64,1}

The `i+1`-th component contains the number of distinct monomials of the
`HomogeneousPolynomial` of order `i`, equivalent to `length(coeff_table[i])`.
""" size_table

@doc """
    pos_table::Array{Dict{Int64,Int64},1}

The `i+1`-th component maps the hash index to the (lexicographic) position
of the corresponding monomial in `coeffs_table`.
""" pos_table
