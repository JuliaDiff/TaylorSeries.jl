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
