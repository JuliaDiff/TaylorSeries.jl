# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# Handles Taylor series of arbitrary but finite order

"""
    TaylorSeries

A Julia package for Taylor expansions in one or more independent variables.

The basic constructors are [`Taylor1`](@ref) and [`TaylorN`](@ref);
see also [`HomogeneousPolynomial`](@ref).

"""
module TaylorSeries


using InteractiveUtils: subtypes
using SparseArrays: SparseMatrixCSC
using Markdown
using Requires

using LinearAlgebra: norm, mul!,
    lu, lu!, LinearAlgebra.lutype, LinearAlgebra.copy_oftype,
    LinearAlgebra.issuccess

import LinearAlgebra: norm, mul!, lu

import Base: ==, +, -, *, /, ^

import Base: iterate, size, eachindex, firstindex, lastindex,
    eltype, length, getindex, setindex!, axes, copyto!

import Base: zero, one, zeros, ones, isinf, isnan, iszero,
    convert, promote_rule, promote, show,
    real, imag, conj, adjoint,
    rem, mod, mod2pi, abs, abs2,
    sqrt, exp, log, sin, cos, sincos, tan,
    asin, acos, atan, sinh, cosh, tanh,
    power_by_squaring,
    rtoldefault, isfinite, isapprox, rad2deg, deg2rad

export Taylor1, TaylorN, HomogeneousPolynomial, AbstractSeries

export getcoeff, derivative, integrate, differentiate,
    evaluate, evaluate!, inverse, set_taylor1_varname,
    show_params_TaylorN, show_monomials, displayBigO, use_show_default,
    get_order, get_numvars,
    set_variables, get_variables,
    get_variable_names, get_variable_symbols,
    # jacobian, hessian, jacobian!, hessian!,
    ∇, taylor_expand, update!, constant_term, linear_polynomial,
    normalize_taylor

include("parameters.jl")
include("hash_tables.jl")
include("constructors.jl")
include("conversion.jl")
include("auxiliary.jl")
include("arithmetic.jl")
include("power.jl")
include("functions.jl")
include("other_functions.jl")
include("evaluate.jl")
include("calculus.jl")
include("dictmutfunct.jl")
include("broadcasting.jl")
include("printing.jl")

function __init__()
    @require IntervalArithmetic = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253" include("intervals.jl")
end

end # module
