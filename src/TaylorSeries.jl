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


using SparseArrays: SparseMatrixCSC
using Markdown

using LinearAlgebra: norm, mul!,
    lu, lu!, LinearAlgebra.lutype, LinearAlgebra.copy_oftype,
    LinearAlgebra.issuccess, NoPivot, RowMaximum

import LinearAlgebra: norm, mul!, lu

import Base: ==, +, -, *, /, ^

import Base: iterate, size, eachindex, firstindex, lastindex,
    eltype, length, getindex, setindex!, axes, copyto!

import Base: zero, one, zeros, ones, isinf, isnan, iszero, isless,
    convert, promote_rule, promote, show, sum!,
    real, imag, conj, adjoint,
    rem, mod, mod2pi, abs, abs2,
    sqrt, exp, expm1, log, log1p,
    sin, cos, sincos, sinpi, cospi, sincospi, tan,
    asin, acos, atan, sinh, cosh, tanh, atanh, asinh, acosh,
    power_by_squaring,
    rtoldefault, isfinite, isapprox, rad2deg, deg2rad

import Base.float

export Taylor1, TaylorN, HomogeneousPolynomial, AbstractSeries, TS

export getcoeff, derivative, integrate, differentiate,
    evaluate, evaluate!, inverse, inverse_map, set_taylor1_varname,
    show_params_TaylorN, show_monomials, displayBigO, use_show_default,
    get_order, get_numvars,
    set_variables, get_variables,
    get_variable_names, get_variable_symbols,
    # jacobian, hessian, jacobian!, hessian!,
    ∇, taylor_expand, update!,
    constant_term, linear_polynomial, nonlinear_polynomial,
    normalize_taylor, norm

const TS = TaylorSeries

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

end # module
