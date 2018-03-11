# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# Handles Taylor series of arbitrary but finite order

__precompile__(true)

"""
    TaylorSeries

A Julia package for Taylor expansions in one or more independent variables.

The basic constructors are [`Taylor1`](@ref) and [`TaylorN`](@ref);
see also [`HomogeneousPolynomial`](@ref).

"""
module TaylorSeries

using Compat

if VERSION <= v"0.7.0-DEV.2004"
    import Base: norm, gradient, A_mul_B!
    const mul! = A_mul_B!
    export mul!
else
    using InteractiveUtils: subtypes
    using SparseArrays: SparseMatrixCSC
    import LinearAlgebra: norm, gradient, mul!
    import Base: lastindex
    using Markdown
end

import Base: ==, +, -, *, /, ^

import Base: zero, one, zeros, ones, isinf, isnan, iszero,
    convert, promote_rule, promote, eltype, length, show,
    real, imag, conj, ctranspose,
    rem, mod, mod2pi, abs, abs2,
    sqrt, exp, log, sin, cos, tan,
    asin, acos, atan, sinh, cosh, tanh,
    power_by_squaring,
    getindex, setindex!, endof,
    rtoldefault, isfinite, isapprox, rad2deg, deg2rad

export Taylor1, TaylorN, HomogeneousPolynomial, AbstractSeries

export getcoeff, derivative, integrate, rectaylor1,
    evaluate, evaluate!, inverse,
    show_params_TaylorN, show_monomials, displayBigO, use_show_default,
    get_order, get_numvars,
    set_variables, get_variables,
    get_variable_names, get_variable_symbols,
    ∇, jacobian, jacobian!, hessian, hessian!,
    taylor_expand, update!, constant_term

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
include("printing.jl")

end # module
