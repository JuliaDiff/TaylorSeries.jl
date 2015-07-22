# TaylorSeries.jl
#
# Julia module for handling Taylor series of arbitrary but finite order
#
# - utils_Taylor1.jl contains the constructors and methods for 1-variable expansions
#
# - utils_TaylorN.jl contains the constructors and methods for N-variable expansions
#
# Last modification: 2015.07.03
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

module TaylorSeries

## Documentation
if VERSION < v"0.4.0-dev"
    using Docile
end

## Compatibility v0.3 -> 0.4
using Compat
@compat sizehint!
@compat trunc
@compat eachindex
@compat round

import Base: zero, one, zeros, ones,
    convert, promote_rule, promote, eltype, length, show,
    real, imag, conj, ctranspose,
    rem, mod, mod2pi, gradient,
    sqrt, exp, log, sin, cos, tan


## Exported types and methods
export Taylor1, TaylorN, HomogeneousPolynomial
export taylor1_variable, taylorN_variable, get_coeff,
    diffTaylor, integTaylor, evaluate, deriv,
    show_params_TaylorN,
    get_order, get_numvars,
    set_variables,
    ∇, jacobian, hessian


# subscriptify is taken from ValidatedNumerics/src/nterval_definition.jl
# and is licensed under MIT "Expat".
# superscriptify is a small variation

const subscript_digits = [c for c in "₀₁₂₃₄₅₆₇₈₉"]
const superscript_digits = [c for c in "⁰¹²³⁴⁵⁶⁷⁸⁹"]

function subscriptify(n::Int)
    dig = reverse(digits(n))
    join([subscript_digits[i+1] for i in dig])
end
function superscriptify(n::Int)
    dig = reverse(digits(n))
    join([superscript_digits[i+1] for i in dig])
end


# one variable:
include("Taylor1.jl")

# several variables:
include("parameters.jl")
include("hash_tables.jl")
include("TaylorN.jl")

include("printing.jl")

end # module
