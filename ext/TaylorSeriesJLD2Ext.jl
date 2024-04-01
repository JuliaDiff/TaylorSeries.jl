module TaylorSeriesJLD2Ext

import Base: convert
using TaylorSeries

if isdefined(Base, :get_extension)
    import JLD2: writeas
else
    import ..JLD2: writeas
end

@doc raw"""
    TaylorNSerialization{T}

Custom serialization struct to save a `TaylorN{T}` to a `.jld2` file.

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

# Tell JLD2 to save TaylorN{T} as TaylorNSerialization{T}
writeas(::Type{TaylorN{T}}) where {T} = TaylorNSerialization{T}

# Convert method to write .jld2 files
function convert(::Type{TaylorNSerialization{T}}, eph::TaylorN{T}) where {T}
    # Variables
    vars = TS.get_variable_names()
    # Number of variables
    n = length(vars)
    # TaylorN order
    varorder = eph.order
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
        TS.set_variables(T, vars, order = varorder)
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

end