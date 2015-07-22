# This file is part of TaylorSeries.jl
#
# Hash tables for HomogeneousPolynomial and TaylorN

@doc """
  Generates the array of dictionaries `index_table`, `size_table` and `pos_table`:

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
function generateTables()
    maxOrd = get_order()
    numVars = get_numvars()

    arrayInd  = Array(Dict{Int,Array{Int,1}},maxOrd+1)
    arraySize = Array(Int,maxOrd+1)
    arrayPos  = Array(Dict{UInt,Int},maxOrd+1)

    if numVars==1
        for k = 0:maxOrd
            dInd = Dict{Int, Array{Int,1}}()
            dPos = Dict{UInt,Int}()
            dInd[1] = [k]
            dPos[hash([k])] = 1
            arrayInd[k+1] = dInd
            arraySize[k+1] = 1
            arrayPos[k+1] = dPos
        end
        return arrayInd, arraySize, arrayPos
    end

    iindices = zeros(Int, numVars)
    idic = 0
    for kDeg = 0:maxOrd
        pos = 0
        dInd = Dict{Int, Array{Int,1}}()
        dPos = Dict{UInt, Int}()
        iV = numVars
        for iz = 0:kDeg
            @inbounds iindices[end] = iz
            pos, dInd[pos] = pos2indices!( iindices, pos, dInd, iV, kDeg )
        end
        nCoefH = length(dInd)
        for i=1:nCoefH
            kdic = hash(dInd[i])
            dPos[kdic] = i
        end
        idic += 1
        arrayInd[idic] = dInd
        arraySize[idic] = nCoefH
        arrayPos[idic] = dPos
    end

    return arrayInd, arraySize, arrayPos
end

function pos2indices!(iIndices::Array{Int,1}, pos::Int, dict::Dict{Int,Array{Int,1}},
    iV::Int, kDeg::Int )

    jVar = iV-1
    kDegNext = kDeg - iIndices[iV]
    if jVar > 1
        for jDeg = 0:kDegNext
            @inbounds iIndices[jVar] = jDeg
            pos, dict[pos] = pos2indices!( iIndices, pos, dict, jVar, kDegNext )
        end
    else
        iIndices[1] = kDegNext
        pos += 1
        @inbounds dict[pos] = iIndices
    end

    return pos, iIndices[1:end]
end

const index_table, size_table, pos_table = generateTables()
gc();

