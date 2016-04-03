# This file is part of the TaylorSeries.jl Julia package, MIT license

# Hash tables for HomogeneousPolynomial and TaylorN

@doc """
  Generates the dictionaries `index_table`, `size_table` and `pos_table`:

  - `index_table`: vector that contains the dictionaries that link the
  lexicographic position of the monomial with the corresponding indexes of
  the powers that characterize the monomial of given degree.
  The vector entry `[k+1]` corresponds to the homogeneous polynomial of
  degree `k`.

  - `size_table`: vector containing the number of distinct monomials of the
  homogeneous polynomial, ordered by the degree of the polynomial.

  - `pos_table`: vector with the inverse of `index_table`, i.e., for a
  given degree `k` and a vector of indexes (hashed), it returns the
  (lexicographic) position of the corresponding monomial.
""" ->

function generate_tables(num_vars, order)

    coeff_table = [generate_index_vectors(num_vars, i) for i in 0:order]

    index_table = Vector{Int}[map(x->in_base(order, x), coeffs) for coeffs in coeff_table]

    pos_table = map(make_inverse_dict, index_table)
    size_table = map(length, index_table)

    coeff_table, index_table, size_table, pos_table

end

@doc "Make a list of index vectors with given number of variables and degree" ->
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

function make_inverse_dict(v::Vector)
    Dict([x=>i for (i,x) in enumerate(v)])
end

@doc "Convert vector v of non-negative integers to base order+1" ->
function in_base(order, v)
    order = order+1

    result = 0

    for i in v
        result = result*order + i
    end

    result
end


const coeff_table, index_table, size_table, pos_table = generate_tables(get_numvars(), get_order())
gc();
