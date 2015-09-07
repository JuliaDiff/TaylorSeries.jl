# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Handles Taylor series of arbitrary but finite order

module TaylorSeries

if VERSION < v"0.4.0-dev"
    using Docile
end

using Compat

import Base: ==, +, -, *, /, ^
import Base: zero, one, zeros, ones, isinf, isnan,
    convert, promote_rule, promote, eltype, length, show,
    real, imag, conj, ctranspose,
    rem, mod, mod2pi, gradient,
    sqrt, exp, log, sin, cos, tan

export Taylor1, TaylorN, HomogeneousPolynomial
export taylor1_variable, taylorN_variable, get_coeff,
    diffTaylor, integTaylor, evaluate, deriv,
    show_params_TaylorN,
    get_order, get_numvars,
    set_variables, get_variables,
    ∇, jacobian, hessian

# one variable:
include("Taylor1.jl")

# several variables:
include("parameters.jl")
include("hash_tables.jl")
include("TaylorN.jl")

include("printing.jl")

end # module
