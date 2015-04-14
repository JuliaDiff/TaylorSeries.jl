# hashtables.jl: hash tables for HomogPol and TaylorN
#
# Last modification: 2015.04.07
#
# Luis Benet & David P. Sanders
# UNAM
#

@doc """Type structure holding the current parameters for `TaylorN` and `HomogPol`.
    
    Fieldnames:

    - maxOrder: maximim degree of the polynomials
    - numVars : maximum number of variables

    These parameters can be changed using `set_ParamsTaylorN(order,numVars)`
    """ ->
type ParamsTaylorN
    maxOrder :: Int
    numVars  :: Int
end

global const _params = ParamsTaylorN(6,2)

@doc """Display the current parameters for `TaylorN` and `HomogPol`""" ->
function show_ParamsTaylorN()
    info( string("`TaylorN` and `HomogPol` parameters:\n    Maximum order       = ", _params.maxOrder,
        "\n    Number of variables = ", _params.numVars) )
    nothing
end

## Hash tables
@doc """
  Generates the array of dictionaries `indicesTable`, `sizeTable` and `posTable`:

  - `indicesTable`: vector that contains the dictionaries that link a (lexicographic) position
  of the monomial with the corresponding indexesfor a homogeneous polynomial of given degree. 
  The vector entry `[k+1]` corresponds to the homogeneous polynomial of degree `k`.
  - `sizeTable`: vector containing the number of distinct monomials of the homogeneous polynomial, 
  ordered by the degree of the polynomial.
  - `posTable`: vector with the inverse of `indicesTable`, i.e., for a given degree `k` and
  a vector of indexes, it returns the (lexicographic) position of the corresponding monomial.
""" ->
function generateTables()
    maxOrd = _params.maxOrder
    numVars = _params.numVars
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

const indicesTable, sizeTable, posTable = generateTables()
gc();


## Utilities to get/set the maximum order and number of variables; they reset the hash tables
get_maxOrder() = _params.maxOrder
set_maxOrder(n::Int) = set_ParamsTaylorN(n, _params.numVars)
get_numVars() = _params.numVars
set_numVars(n::Int) = set_ParamsTaylorN(_params.maxOrder, n)

function set_ParamsTaylorN(order::Int, numVars::Int)
    (order > 0 && numVars>=1) || error("Incorrect order or number of variables")
    order == _params.maxOrder && numVars == _params.numVars && return order, numVars
    oldOrder = _params.maxOrder
    oldVars = _params.numVars
    global _params = ParamsTaylorN(order, numVars)

    resize!(indicesTable,order+1)
    resize!(sizeTable,order+1)
    resize!(posTable,order+1)
    indicesTable[:], sizeTable[:], posTable[:] = generateTables()
    gc();
    show_ParamsTaylorN()

    return order, numVars
end
