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

if !isdefined(Base, :get_extension)
    using Requires
end

using LinearAlgebra: norm, mul!,
    lu, lu!, LinearAlgebra.lutype, LinearAlgebra.copy_oftype,
    LinearAlgebra.issuccess

if VERSION >= v"1.7.0-DEV.1188"
    using LinearAlgebra: NoPivot, RowMaximum
end

import LinearAlgebra: norm, mul!, lu

import Base: ==, +, -, *, /, ^

import Base: iterate, size, eachindex, firstindex, lastindex,
    eltype, length, getindex, setindex!, axes, copyto!

import Base: zero, one, zeros, ones, isinf, isnan, iszero, isless,
    convert, promote_rule, promote, show,
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
    evaluate, evaluate!, inverse, set_taylor1_varname,
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

function __init__()
    @static if !isdefined(Base, :get_extension)
        @require IntervalArithmetic = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253" begin
            include("../ext/TaylorSeriesIAExt.jl")
        end
        @require StaticArrays = "90137ffa-7385-5640-81b9-e52037218182" begin
            include("../ext/TaylorSeriesSAExt.jl")
        end
        @require JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819" begin
            include("../ext/TaylorSeriesJLD2Ext.jl")
        end
    end
end

end # module
