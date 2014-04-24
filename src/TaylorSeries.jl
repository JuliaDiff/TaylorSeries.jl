# TaylorSeries.jl
#
# Julia module for handling Taylor series of arbitrary but finite order
#
# utils_Taylor1.jl contains the constructors and methods for 1-variable expansions
#
# utils_TaylorN.jl contains the constructors and methods for N-variable expansions
#
# Last modification: 2014.04.12
#
# Luis Benet & David P. Sanders
# UNAM
#

module TaylorSeries

import Base: zero, one
import Base: convert, promote_rule, promote, eltype, length, showcompact
import Base: real, imag, conj, ctranspose
import Base: sqrt, exp, log, sin, cos, tan#, square

include("utils_Taylor1.jl")

include("utils_TaylorN.jl")

# Exports to Taylor
export Taylor, diffTaylor, integTaylor, evalTaylor, deriv

# Exports to TaylorN
export TaylorN
export set_maxOrder, get_maxOrder, set_numVars, get_numVars

end