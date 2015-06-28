# hashtables.jl: hash tables for HomogeneousPolynomial and TaylorN
#
# Last modification: 2015.06.28
#
# Luis Benet & David P. Sanders
# UNAM


@doc """
  Generates the dictionaries `index_table`, `size_table` and `pos_table`:

  - `index_table`: vector that contains the dictionaries that link the
  lexicographic position of the monomial with the corresponding indexes of
  the powers that characterize the monomial of given degree.
  The vector entry `[k+1]` corresponds to the homogeneous polynomial of
  degree `k`.
  - `size_table`: vector containing the number of distinct monomials of the
  homogeneous polynomial, ordered by the degree of the polynomial.
  - `pos_table`: vector with the inverse of `indicesTable`, i.e., for a
  given degree `k` and a vector of indexes (hashed), it returns the
  (lexicographic) position of the corresponding monomial.
""" ->

function generate_tables(num_vars=2, order=6)

    index_table = Vector{NTuple{num_vars, Int}}[generate_index_vectors(num_vars, i) for i in 0:order]
    pos_table = map(inverse_hash_map, index_table)
    size_table = map(length, index_table)

    index_table, size_table, pos_table

end

function generate_index_vectors(num_vars, degree)

    if num_vars == 1
        return [(degree,)]
    end

    indices = NTuple{num_vars, Int}[]

    for k in degree:-1:0

        new_indices = [tuple(k, x...) for x in generate_index_vectors(num_vars-1, degree-k)]
        append!(indices, new_indices)

    end

    indices

end

function inverse_hash_map(v::Vector)
    Dict([x=>i for (i,x) in enumerate(v)])
end




index_table, size_table, pos_table = generate_tables()
gc();

